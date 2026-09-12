"""Render the README tables from the run outputs; text between <!-- table:NAME --> and <!-- /table --> is replaced."""
import glob, json, os, re, sys
import numpy as np
sys.path.insert(0, "scripts")
from ages import tables as age_tables

RUNS = "out/runs"; README = "README.md"
SUBS = ["hg002", "chimp", "bonobo", "gorilla", "borang", "sorang", "siamang"]
SUBL = {"hg002": "HG002", "chimp": "chimp", "bonobo": "bonobo", "gorilla": "gorilla", "borang": "B. orang", "sorang": "S. orang", "siamang": "siamang", "random": "random"}
CL = {"unique": "non-repeat, non-SD", "SD": "segmental duplication", "ct": "centromere transition", "mon": "monomeric alpha-sat", "dhor": "diverged HOR",
      "hor": "inactive HOR", "hor_active": "active HOR", "bsat": "beta satellite", "gsat": "gamma satellite", "hsat1A": "HSat1A", "hsat1B": "HSat1B",
      "hsat2": "HSat2", "hsat3": "HSat3", "censat": "other censat", "rDNA": "rDNA", "ALL": "**whole genome**"}

def classes(path):
    rows = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for l in f:
            v = l.rstrip("\n").split("\t"); rows[v[0]] = dict(zip(hdr[1:], map(int, v[1:])))
    return rows

def pct(a, b): return f"{100*a/b:.1f}" if b else "--"

def t_k31():
    r = classes(f"{RUNS}/chm13_k31_all.classes.tsv")
    out = ["| region class | Mb | " + " | ".join(SUBL[s] for s in SUBS) + " | any ape | human-specific | CHM13-only | single-copy in CHM13 |", "|" + "---|" * (len(SUBS) + 6)]
    for c, row in r.items():
        v = row["valid_positions"]
        if v == 0: continue
        out.append(f"| {CL[c]} | {v/1e6:.1f} | " + " | ".join(pct(row[f'shared_{s}'], v) for s in SUBS) + f" | {pct(row['ape_any'], v)} | {pct(row['human_specific'], v)} | {pct(row['chm13_only'], v)} | {pct(row['unique_in_query'], v)} |")
    return "\n".join(out)

def t_ages():
    out_t, labels = age_tables(f"{RUNS}/chm13_k31_all")
    cols = ["chm13_only", "human_specific", "pan", "gorilla", "pongo", "hylobatid"]
    names = {"chm13_only": "CHM13 only", "human_specific": "human only (in HG002)", "pan": "to Pan (>=6 My)", "gorilla": "to gorilla (>=9 My)", "pongo": "to orangutan (>=15 My)", "hylobatid": "to siamang (>=20 My)"}
    out = ["| region class | Mb | " + " | ".join(names[c] for c in cols) + " | gorilla not Pan | chimp not bonobo | bonobo not chimp |", "|" + "---|" * (len(cols) + 5)]
    for c, r in out_t.items():
        v = r["positions"]
        out.append(f"| {CL[c]} | {v/1e6:.1f} | " + " | ".join(pct(r[k], v) for k in cols) + f" | {pct(r['gorilla_not_pan'], v)} | {pct(r['chimp_not_bonobo'], v)} | {pct(r['bonobo_not_chimp'], v)} |")
    return "\n".join(out)

def ladder_runs():
    runs = {}
    for p in glob.glob(f"{RUNS}/chm13_k*_all.classes.tsv"):
        k = int(os.path.basename(p).split("_k")[1].split("_")[0]); runs[k] = p[:-len(".classes.tsv")]
    return dict(sorted(runs.items()))

def t_ladder():
    runs = ladder_runs(); ks = list(runs)
    show = [("unique", "chimp"), ("unique", "siamang"), ("unique", "random"), ("SD", "chimp"), ("hor_active", "chimp"), ("hor_active", "siamang"), ("hor", "chimp"), ("mon", "chimp"), ("hsat3", "chimp"), ("hsat3", "siamang"), ("hsat2", "chimp"), ("hsat1B", "chimp"), ("bsat", "chimp"), ("censat", "chimp")]
    out = ["| class, genome | " + " | ".join(f"k={k}" for k in ks) + " |", "|" + "---|" * (len(ks) + 1)]
    for c, s in show:
        cells = []
        for k in ks:
            r = classes(runs[k] + ".classes.tsv").get(c)
            cells.append(pct(r[f"shared_{s}"], r["valid_positions"]) if r and f"shared_{s}" in r else "--")
        out.append(f"| {CL[c]}, {SUBL[s]} | " + " | ".join(cells) + " |")
    return "\n".join(out)

def t_census():
    runs = ladder_runs()
    bins = [(1, 1), (2, 2), (3, 9), (10, 99), (100, 999), (1000, 10**12)]
    names = ["1 (single copy)", "2", "3-9", "10-99", "100-999", ">=1000"]
    out = ["| k | distinct k-mers (G) | " + " | ".join(names) + " |", "|" + "---|" * (len(names) + 2)]
    for k, p in runs.items():
        m = json.load(open(p + ".json")); h = np.array(m["self_histogram"], dtype=float); tot = (h[:, 0] * h[:, 1]).sum()
        cells = [pct(((h[:, 0] >= lo) & (h[:, 0] <= hi)).dot(h[:, 0] * h[:, 1]), tot) for lo, hi in bins]
        out.append(f"| {k} | {m['query_distinct']/1e9:.3f} | " + " | ".join(cells) + " |")
    return "\n".join(out)

def t_divergence():
    """Apparent per-base divergence from k-mer survival in non-repeat sequence: d = 1 - f^(1/k) at k=31, and the slope of ln f against k over the ladder."""
    runs = ladder_runs()
    lit = {"hg002": "~0.1 (heterozygosity)", "chimp": "1.23 (CSAC 2005)", "bonobo": "1.3 (Pruefer 2012)", "gorilla": "1.75 (Scally 2012)", "borang": "3.1 (Locke 2011)", "sorang": "3.1 (Locke 2011)", "siamang": "not tabulated here"}
    out = ["| genome | 31-mers shared, non-repeat (%) | d from k=31 alone (%) | d from the slope of ln f over k=31..201 (%) | published single-nucleotide divergence (%) |", "|---|---|---|---|---|"]
    r31 = classes(f"{RUNS}/chm13_k31_all.classes.tsv")["unique"]
    for s in SUBS:
        f31 = r31[f"shared_{s}"] / r31["valid_positions"]
        ks, ys = [], []
        for k, p in runs.items():
            if 31 <= k <= 201:
                r = classes(p + ".classes.tsv")["unique"]; ks.append(k); ys.append(np.log(r[f"shared_{s}"] / r["valid_positions"]))
        slope = np.polyfit(ks, ys, 1)[0] if len(ks) >= 2 else np.nan
        out.append(f"| {SUBL[s]} | {100*f31:.2f} | {100*(1-f31**(1/31)):.3f} | {100*(1-np.exp(slope)):.3f} | {lit.get(s, '')} |")
    return "\n".join(out)

def t_within():
    r = classes(f"{RUNS}/hg002mat_k31_patchm13.auto.classes.tsv")
    out = ["| region class (HG002 maternal autosomes) | Mb | in paternal haplotype | in CHM13 | in either | single-copy in maternal |", "|---|---|---|---|---|---|"]
    for c, row in r.items():
        v = row["valid_positions"]
        if v == 0: continue
        out.append(f"| {CL[c]} | {v/1e6:.1f} | {pct(row['shared_hg002_pat'], v)} | {pct(row['shared_chm13'], v)} | {pct(row['ape_any'], v)} | {pct(row['unique_in_query'], v)} |")
    return "\n".join(out)

def t_compute():
    out = ["| run | k | stride | query parts | wall (s) | peak footprint (GB) |", "|---|---|---|---|---|---|"]
    log = open("out/run_all.log").read() if os.path.exists("out/run_all.log") else ""
    peaks = {}
    for m in re.finditer(r"k=(\d+) stride=(\d+).*?done in ([\d.]+)s.*?(\d+)\s+peak memory footprint", log, re.S):
        peaks[int(m.group(1))] = (float(m.group(3)), int(m.group(4)) / 1e9)
    for k, p in ladder_runs().items():
        m = json.load(open(p + ".json"))
        wall, peak = peaks.get(k, (m["seconds"], float("nan")))
        out.append(f"| {os.path.basename(p)} | {k} | {m['stride']} | {m['query_parts']} | {wall:.0f} | {peak:.1f} |")
    return "\n".join(out)

def t_arrays():
    arr = json.load(open(f"{RUNS}/chm13_k31_all.arrays.json"))
    subs = ["hg002", "chimp", "bonobo", "gorilla", "borang", "sorang", "siamang"]
    out = ["| active HOR array | Mb | single-copy 31-mers (%) | " + " | ".join(SUBL[s] for s in subs) + " | highest ape |", "|" + "---|" * (len(subs) + 4)]
    for a in arr:
        out.append(f"| {a['chrom']}:{a['start']/1e6:.2f}-{a['end']/1e6:.2f} | {a['n']/1e6:.2f} | {a['single_copy_pct']:.1f} | " + " | ".join(f"{a['shared_pct'][s]:.1f}" for s in subs) + f" | {SUBL[a['closest_ape']]} |")
    return "\n".join(out)

def fill(text, name, body):
    pat = re.compile(rf"(<!-- table:{name} -->\n).*?(<!-- /table -->)", re.S)
    if not pat.search(text): raise SystemExit(f"marker {name} missing")
    return pat.sub(lambda m: m.group(1) + body + "\n" + m.group(2), text)

if __name__ == "__main__":
    text = open(README).read()
    for name, fn in [("k31", t_k31), ("ages", t_ages), ("ladder", t_ladder), ("census", t_census), ("divergence", t_divergence), ("within", t_within), ("compute", t_compute), ("arrays", t_arrays)]:
        try:
            text = fill(text, name, fn())
        except Exception as e:
            print(f"{name}: skipped ({e})")
    open(README, "w").write(text)
    print("README tables updated")
