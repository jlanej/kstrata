"""Figures for the k-mer strata resource. Reads the outputs of strata.py (classes.tsv, win*.tsv.gz,
summary.json) and the run json files (self histograms)."""
import glob, gzip, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "out/runs"; FIG = "fig"; os.makedirs(FIG, exist_ok=True)
SUBJ_ORDER = ["hg002", "chimp", "bonobo", "gorilla", "borang", "sorang", "siamang"]
SUBJ_LABEL = {"hg002": "HG002 (human)", "chimp": "chimpanzee", "bonobo": "bonobo", "gorilla": "gorilla",
              "borang": "Bornean orangutan", "sorang": "Sumatran orangutan", "siamang": "siamang"}
# approximate divergence times (Mya) used only for the x axis of the decay figure
DIV_MYA = {"hg002": 0.0, "chimp": 6.4, "bonobo": 6.4, "gorilla": 9.0, "borang": 15.5, "sorang": 15.5, "siamang": 19.5}
CLASS_LABEL = {"unique": "non-repeat, non-SD", "SD": "segmental dup.", "ct": "centromere transition", "mon": "monomeric alpha-sat.",
               "dhor": "diverged HOR", "hor": "inactive HOR", "hor_active": "active HOR", "bsat": "beta satellite", "gsat": "gamma satellite",
               "hsat1A": "HSat1A", "hsat1B": "HSat1B", "hsat2": "HSat2", "hsat3": "HSat3", "censat": "other censat", "rDNA": "rDNA"}

def read_classes(path):
    rows = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for l in f:
            v = l.rstrip("\n").split("\t"); rows[v[0]] = dict(zip(hdr[1:], map(int, v[1:])))
    return hdr, rows

def fig_class_heatmap(run):
    hdr, rows = read_classes(run + ".classes.tsv")
    subs = [c[len("shared_"):] for c in hdr if c.startswith("shared_")]
    subs = [s for s in SUBJ_ORDER if s in subs]
    classes = [c for c in CLASS_LABEL if c in rows and rows[c]["valid_positions"] > 0]
    M = np.array([[rows[c][f"shared_{s}"] / max(rows[c]["valid_positions"], 1) for s in subs] for c in classes])
    fig, ax = plt.subplots(figsize=(1.2 + 0.9 * len(subs), 0.42 * len(classes) + 1.2))
    im = ax.imshow(M, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(subs))); ax.set_xticklabels([SUBJ_LABEL[s] for s in subs], rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(classes))); ax.set_yticklabels([f"{CLASS_LABEL[c]} ({rows[c]['valid_positions']/1e6:.1f} Mb)" for c in classes], fontsize=8)
    for i in range(len(classes)):
        for j in range(len(subs)):
            ax.text(j, i, f"{100*M[i,j]:.0f}", ha="center", va="center", fontsize=7, color="white" if M[i, j] < 0.6 else "black")
    ax.set_title("Fraction of CHM13 31-mers present in each genome, by region class", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_class_sharing_k31.png", dpi=150); plt.close(fig)

def fig_multik(runs_by_k):
    ks = sorted(runs_by_k)
    subs_show = ["hg002", "chimp", "gorilla", "sorang", "siamang", "random"]
    classes = ["unique", "SD", "mon", "hor_active", "hor", "dhor", "hsat1B", "hsat2", "hsat3", "bsat", "censat", "rDNA"]
    tables = {k: read_classes(runs_by_k[k] + ".classes.tsv")[1] for k in ks}
    fig, axes = plt.subplots(3, 4, figsize=(13, 8.5), sharex=True, sharey=True)
    for ax, c in zip(axes.flat, classes):
        for s in subs_show:
            ys = []
            for k in ks:
                r = tables[k].get(c)
                ys.append(r[f"shared_{s}"] / r["valid_positions"] if r and r["valid_positions"] and f"shared_{s}" in r else np.nan)
            style = dict(marker="o", ms=3, label=SUBJ_LABEL.get(s, s))
            if s == "random": style = dict(ls="--", color="grey", lw=1, label="uniform random 3.1 Gb (chance)")
            ax.plot(ks, ys, **style)
        ax.set_xscale("log"); ax.set_title(CLASS_LABEL[c], fontsize=9); ax.grid(alpha=0.3)
        ax.set_ylim(0, 1.02)
        ax.set_xticks(ks); ax.set_xticklabels([str(k) for k in ks], fontsize=7, rotation=60); ax.minorticks_off()
    axes[0, 0].legend(fontsize=7)
    for ax in axes[-1]: ax.set_xlabel("k")
    for ax in axes[:, 0]: ax.set_ylabel("fraction of CHM13 k-mers shared")
    fig.suptitle("Exact k-mer sharing as a function of k, by region class and genome", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_multik_sharing.png", dpi=150); plt.close(fig)

def fig_decay(run):
    """Sharing vs divergence time at k=31 per class, with the apparent per-base divergence 1-f^(1/k)."""
    hdr, rows = read_classes(run + ".classes.tsv"); k = 31
    subs = [s for s in SUBJ_ORDER if f"shared_{s}" in hdr]
    classes = ["unique", "SD", "ct", "mon", "dhor", "hor", "hor_active", "hsat1B", "hsat2", "hsat3", "bsat", "censat"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    cmap = plt.get_cmap("tab20")
    for i, c in enumerate(classes):
        if c not in rows or rows[c]["valid_positions"] == 0: continue
        xs = [DIV_MYA[s] for s in subs]; ys = [rows[c][f"shared_{s}"] / rows[c]["valid_positions"] for s in subs]
        ax.plot(xs, ys, marker="o", ms=4, color=cmap(i), label=CLASS_LABEL[c])
    ax.set_xlabel("approximate divergence time from human (Mya)"); ax.set_ylabel("fraction of CHM13 31-mers present")
    ax.set_yscale("log"); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7, ncol=2)
    ax.set_title("Exact 31-mer survival across the apes, by region class", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_decay_k31.png", dpi=150); plt.close(fig)

def read_windows(path):
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in f]
    chrom = np.array([r[0] for r in rows]); start = np.array([int(r[1]) for r in rows]); end = np.array([int(r[2]) for r in rows])
    cols = {h: np.array([int(r[i]) for r in rows]) for i, h in enumerate(hdr) if i >= 3}
    return chrom, start, end, cols

def fig_centromeres(run, chroms=("chr1", "chr8", "chr11", "chr17", "chrX"), iv=None):
    chrom, start, end, cols = read_windows(run + ".win10000.tsv.gz")
    subs = ["hg002", "chimp", "gorilla", "sorang", "siamang"]
    fig, axes = plt.subplots(len(chroms), 1, figsize=(12, 2.6 * len(chroms)))
    for ax, ch in zip(np.atleast_1d(axes), chroms):
        m = chrom == ch
        from strata import CLASS_ID
        cen = [(s, e, c) for s, e, c in iv.get(ch, [])] if iv else []
        sat = [(s, e, c) for s, e, c in cen if c != CLASS_ID["SD"]]
        if sat:
            lo = max(min(s for s, e, c in sat) - 1_000_000, 0); hi = max(e for s, e, c in sat) + 1_000_000
        else:
            lo, hi = start[m].min(), end[m].max()
        sel = m & (start >= lo) & (end <= hi)
        x = (start[sel] + end[sel]) / 2 / 1e6; v = np.maximum(cols["valid"][sel], 1)
        for s in subs:
            ax.plot(x, cols[s][sel] / v, lw=0.9, label=SUBJ_LABEL[s])
        ax.set_xlim(lo / 1e6, hi / 1e6); ax.set_ylim(-0.12, 1.02); ax.set_ylabel("shared fraction"); ax.set_title(f"{ch}", fontsize=9, loc="left")
        # annotation bar
        from strata import CLASS_ORDER
        palette = {"ct": "#bbbbbb", "mon": "#f4a582", "dhor": "#d6604d", "hor": "#b2182b", "hor_active": "#67001f", "bsat": "#4393c3", "gsat": "#2166ac",
                   "hsat1A": "#a6dba0", "hsat1B": "#5aae61", "hsat2": "#1b7837", "hsat3": "#00441b", "censat": "#9970ab", "rDNA": "#e7d4e8", "SD": "#fee08b"}
        for s, e, c in cen:
            name = CLASS_ORDER[c]
            ax.add_patch(plt.Rectangle((s / 1e6, -0.11), (e - s) / 1e6, 0.09, color=palette.get(name, "#000000"), lw=0))
        ax.grid(alpha=0.2)
    np.atleast_1d(axes)[0].legend(fontsize=7, ncol=5, loc="upper right")
    np.atleast_1d(axes)[-1].set_xlabel("position (Mb); bar: censat class (grey ct, reds alpha-satellite, greens HSat, blues beta/gamma, purple other, yellow SD)")
    fig.suptitle("k-mer stratigraphy across centromeres: fraction of CHM13 31-mers (10 kb windows) present in each genome", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_centromere_strata_k31.png", dpi=150); plt.close(fig)

def fig_genome_heatmap(run, window=100000):
    chrom, start, end, cols = read_windows(run + f".win{window}.tsv.gz")
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    chroms = [c for c in chroms if (chrom == c).any()]
    maxlen = max(end[chrom == c].max() for c in chroms)
    nb = int(np.ceil(maxlen / window))
    layers = {"chimp": "shared with chimpanzee", "siamang": "shared with siamang", "human_specific": "human-specific (in HG002, in no ape)"}
    fig, axes = plt.subplots(len(layers), 1, figsize=(14, 3.2 * len(layers)))
    for ax, (key, title) in zip(axes, layers.items()):
        M = np.full((len(chroms), nb), np.nan)
        for i, c in enumerate(chroms):
            m = chrom == c
            v = cols["valid"][m]; ok = v > 0
            idx = (start[m] // window)[ok]
            M[i, idx] = cols[key][m][ok] / v[ok]
        im = ax.imshow(M, aspect="auto", cmap="magma", interpolation="nearest", vmin=0, vmax=np.nanpercentile(M, 99))
        ax.set_yticks(range(len(chroms))); ax.set_yticklabels(chroms, fontsize=7)
        ax.set_xticks(np.arange(0, nb, 250e6 / window / 5)); ax.set_xticklabels([f"{int(x*window/1e6)}" for x in ax.get_xticks()], fontsize=7)
        ax.set_title(f"fraction of CHM13 31-mers {title}, {window//1000} kb windows", fontsize=9, loc="left")
        fig.colorbar(im, ax=ax, fraction=0.015, pad=0.01)
    axes[-1].set_xlabel("position (Mb)")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_genome_strata_k31.png", dpi=130); plt.close(fig)

def fig_spectrum(run_jsons):
    """Self multiplicity spectrum of CHM13 across k: fraction of genome positions whose k-mer occurs m times."""
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    edges = [1, 2, 3, 5, 10, 30, 100, 1000, 10000, 10**9]
    labels = ["1", "2", "3-4", "5-9", "10-29", "30-99", "100-999", "1k-9,999", ">=10k"]
    ks = sorted(run_jsons)
    W = np.zeros((len(ks), len(labels)))
    for i, k in enumerate(ks):
        meta = json.load(open(run_jsons[k]))
        h = np.array(meta["self_histogram"], dtype=float)
        tot = (h[:, 0] * h[:, 1]).sum()
        for j in range(len(labels)):
            m = (h[:, 0] >= edges[j]) & (h[:, 0] < edges[j + 1])
            W[i, j] = (h[m, 0] * h[m, 1]).sum() / tot
    bottom = np.zeros(len(ks)); cmap = plt.get_cmap("magma_r")
    for j in range(len(labels)):
        ax.bar(range(len(ks)), W[:, j], bottom=bottom, color=cmap((j + 1) / (len(labels) + 1)), label=labels[j], width=0.8)
        bottom += W[:, j]
    ax.set_xticks(range(len(ks))); ax.set_xticklabels([str(k) for k in ks]); ax.set_xlabel("k"); ax.set_ylabel("fraction of genome positions")
    ax.set_ylim(0, 1); ax.legend(title="copies of the k-mer in CHM13", fontsize=7, title_fontsize=7, ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.set_title("Repeat census of T2T-CHM13 across scales: multiplicity of the k-mer at each position", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_repeat_census.png", dpi=150); plt.close(fig)

if __name__ == "__main__":
    sys.path.insert(0, "scripts")
    from strata import load_intervals
    flag = sys.argv[1] if len(sys.argv) > 1 else "all"
    main = f"{OUT}/chm13_k31_all"
    if flag in ("all", "k31") and os.path.exists(main + ".classes.tsv"):
        fig_class_heatmap(main); fig_decay(main); fig_centromeres(main, iv=load_intervals("chm13"))
        if os.path.exists(main + ".win100000.tsv.gz"): fig_genome_heatmap(main)
    runs = {}
    for p in glob.glob(f"{OUT}/chm13_k*_all.classes.tsv"):
        k = int(os.path.basename(p).split("_k")[1].split("_")[0]); runs[k] = p[:-len(".classes.tsv")]
    if flag in ("all", "multik") and len(runs) > 2: fig_multik(runs)
    js = {k: p + ".json" for k, p in runs.items()}
    if flag in ("all", "spectrum") and len(js) > 2: fig_spectrum(js)
    print("figures:", sorted(os.listdir(FIG)))
