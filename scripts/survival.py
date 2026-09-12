"""Survival length of exact k-mers: for every base of the spotlight regions and each genome, the largest k on the
ladder at which the k-mer starting there still exists in that genome. Writes docs/data/surv_<spot>.bin.gz
(one byte per base per genome, genome-major: the index into the ladder, 0 = not even 16; 255 = censored, no k-mer)
and the figure fig/fig_survival.png."""
import gzip, json, os, sys
import numpy as np
sys.path.insert(0, "scripts")
KS = [16, 21, 25, 31, 41, 51, 61, 81, 101, 151, 201, 301, 501, 701, 1001, 1501, 2001]
R = "out/runs"; OUT = "docs/data"
meta = json.load(open(f"{OUT}/meta.json")); ucs = json.load(open(f"{OUT}/ucs.json"))
subs = [s["name"] for s in json.load(open(f"{R}/regions_k16.json"))["subjects"]]
regions = [l.rstrip("\n").split("\t") for l in open("out/regions.tsv")]  # id chrom start end padstart padend
idx = {}
for l in open("out/prep/regions.idx"):
    n, s, ln = l.rstrip("\n").split("\t"); idx[n] = (int(s), int(ln))
pres = {k: np.fromfile(f"{R}/regions_k{k}.pres.u8", np.uint8) for k in KS}
# validity: the k-mer starting at a position must fit inside its (padded) region; CHM13v2.0 has no N gaps, so that is the only condition

def write_gz(path, data):
    if os.path.exists(path):
        try:
            with gzip.open(path, "rb") as f:
                if f.read() == data: return
        except OSError: pass
    with gzip.GzipFile(path, "wb", compresslevel=6, mtime=0) as f: f.write(data)

surv_all = {}
for sid, ch, s, e, ps, pe in regions:
    s, e, ps, pe = map(int, (s, e, ps, pe)); name = f"{ch}:{ps+1}-{pe}"; off, ln = idx[name]
    e = min(e, ps + ln); L = e - s; a = s - ps  # offset of the unpadded region inside the padded contig (clipped at the chromosome end)
    surv = np.zeros((len(subs), L), np.uint8)
    valid_any = np.zeros(L, bool)
    pos = a + np.arange(L)
    for ki, k in enumerate(KS):
        p = pres[k][off + a: off + a + L]; v = pos + k <= ln  # valid k-mer: fits inside the padded region
        if ki == 0: valid_any = v.copy()
        for j in range(len(subs)):
            hit = v & (((p >> j) & 1) == 1)
            surv[j][hit] = ki + 1           # monotone: the last hit along the ladder is the largest k
    # censored: the k-mer at the base could not be tested beyond the padded end (only the last 2 kb of a region)
    surv[:, ~valid_any] = 255
    write_gz(f"{OUT}/surv_{sid}.bin.gz", surv.tobytes())
    surv_all[sid] = (ch, s, e, surv)
json.dump({"ks": KS, "subjects": subs, "note": "surv_<spot>.bin.gz: genome-major bytes, value = 1-based index into ks of the largest k at which the k-mer starting at the base exists in that genome; 0 = not even k=16; 255 = untestable"}, open(f"{OUT}/survival_meta.json", "w"))
print("survival tracks written for", len(surv_all), "spotlights", file=sys.stderr)

# ---- figure: six regions, one strip per genome, colour = survival length ----
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib.colors import ListedColormap, BoundaryNorm
SUBJ_LAB = {"hg002": "HG002 (human)", "chimp": "chimpanzee", "bonobo": "bonobo", "gorilla": "gorilla", "borang": "B. orangutan", "sorang": "S. orangutan", "siamang": "siamang"}
panels = [
    ("chr20_unique", 20_100_000, 20_130_000, "ordinary non-repeat sequence (chr20)"),
    ("chr1_hor_start", 121_830_000, 121_860_000, "active HOR array, chr1 (kinetochore-bound alpha-satellite)"),
    ("chr8_layers", 43_800_000, 43_830_000, "monomeric alpha-satellite layer beside the chr8 array"),
    ("chr16_hsat2", 40_100_000, 40_130_000, "HSat2 block, chr16"),
    ("hsat1B_block", None, None, "HSat1B block, chrY"),
    ("chr9_hsat3", 50_100_000, 50_130_000, "HSat3 block, chr9"),
]
pax = [u for u in ucs if u["chrom"] == "chr11" and abs(u["start"] - 31_899_884) < 100]
if pax: panels.append((pax[0]["spot"], pax[0]["spot_start"] + 1000, pax[0]["spot_end"] - 1000, "ultraconserved element in the PAX6 region, chr11 (1.1 kb identical in all seven genomes)"))
edges = [0, 1, 2, 4, 6, 8, 10, 12, 14, 18]  # ladder-index bins: absent, 16-20, 21-30, 31-50, 51-80, 81-150, 151-300, 301-700, 701-2001
labels = ["absent at k=16", "16-20", "21-30", "31-50", "51-80", "81-150", "151-300", "301-700", "701-2001"]
cmap = ListedColormap(["#f0f0f0", "#fee8c8", "#fdbb84", "#fc8d59", "#e34a33", "#b30000", "#7a0177", "#3f007d", "#08306b"])
norm = BoundaryNorm(edges, cmap.N); cmap.set_bad('#ffffff')
fig, axes = plt.subplots(len(panels), 1, figsize=(13, 1.9 * len(panels)), gridspec_kw={"hspace": 0.9, "top": 0.9, "bottom": 0.08})
for ax, (sid, ws, we, title) in zip(axes, panels):
    ch, s, e, surv = surv_all[sid]
    if ws is None: ws, we = s + (e - s) // 2 - 15_000, s + (e - s) // 2 + 15_000
    a, b = ws - s, we - s; sub = surv[:, a:b].astype(float); sub[sub == 255] = np.nan
    # downsample by the median over 10-base bins for the raster
    n = (b - a) // 10 * 10; block = sub[:, :n].reshape(len(subs), n // 10, 10); med = np.nanmedian(block, axis=2)
    ax.imshow(med, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest", extent=[ws / 1e6, (ws + n) / 1e6, len(subs) - 0.5, -0.5])
    ax.set_yticks(range(len(subs))); ax.set_yticklabels([SUBJ_LAB[x] for x in subs], fontsize=8)
    ax.set_title(f"{title}: {ch}:{ws:,}-{we:,}", fontsize=9, loc="left"); ax.set_xlabel("position (Mb)", fontsize=8); ax.tick_params(axis="x", labelsize=8)
    ax.ticklabel_format(axis="x", useOffset=False, style="plain"); ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.4f}".rstrip("0").rstrip(".")))
    for y in range(len(subs)): ax.axhline(y + 0.5, color="white", lw=0.8)
handles = [plt.Rectangle((0, 0), 1, 1, color=cmap(i)) for i in range(cmap.N)]
fig.legend(handles, labels, title="longest exact k-mer starting at the base that exists in the genome (k)", ncol=9, fontsize=8, title_fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.005), frameon=False)
fig.suptitle("Survival length of exact k-mers: for each base (columns, 10 bp bins) and genome (rows), the longest k-mer starting there that still exists", fontsize=10, y=0.93)
fig.savefig("fig/fig_survival.png", dpi=150, bbox_inches="tight"); print("figure written", file=sys.stderr)
# summary numbers for the write-up: median survival per panel and genome
for sid, ws, we, title in panels:
    ch, s, e, surv = surv_all[sid]
    if ws is None: ws, we = s + (e - s) // 2 - 15_000, s + (e - s) // 2 + 15_000
    sub = surv[:, ws - s: we - s]
    med = [KS[int(np.median(sub[j][sub[j] != 255])) - 1] if np.median(sub[j][sub[j] != 255]) >= 1 else 0 for j in range(len(subs))]
    print(f"{title[:45]:45s} median survival: " + " ".join(f"{SUBJ_LAB[x].split()[0]}={m}" for x, m in zip(subs, med)))
