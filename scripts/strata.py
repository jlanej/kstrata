"""Summaries of kstrata runs: per-class sharing tables, presence-pattern tables, windowed tracks.

A run is <prefix>.pres.u8 (one byte per query position/stride: bit i set when subject i carries
the k-mer), <prefix>.mult.u8 (multiplicity of the k-mer in the query, 0 = no valid k-mer at that
position) and <prefix>.json (metadata: k, stride, subjects in bit order).

Classes for CHM13v2.0: censat v2.1 class (prefix of the label; 'hor' split into hor_active for
live arrays, label ending in 'L)') painted over the SD annotation, everything else 'unique'.
"""
import gzip, json, os, re, sys
from collections import defaultdict
import numpy as np

CLASS_ORDER = ["unique", "SD", "ct", "mon", "dhor", "hor", "hor_active", "bsat", "gsat", "hsat1A", "hsat1B", "hsat2", "hsat3", "censat", "rDNA"]
CLASS_ID = {c: i for i, c in enumerate(CLASS_ORDER)}

def read_idx(prefix):
    contigs = []
    for l in open(prefix + ".idx"):
        n, s, ln = l.rstrip("\n").split("\t"); contigs.append((n, int(s), int(ln)))
    return contigs

def censat_class(label, genome):
    base = label.split("(")[0]
    fam = label[len(base):]
    if genome == "chm13":
        cls = base.split("_")[0]
        if cls == "hor":
            return "hor_active" if re.search(r"L\)$|L,", fam) else "hor"
        return cls
    # hg002 cenSat v2.0 labels: active_hor(...), hor(...), cenSat(...), HSat3, mixedAlpha, ...
    m = {"active_hor": "hor_active", "hor": "hor", "dhor": "dhor", "mon": "mon", "ct": "ct", "bSat": "bsat", "gSat": "gsat",
         "HSat1A": "hsat1A", "HSat1B": "hsat1B", "HSat2": "hsat2", "HSat3": "hsat3", "cenSat": "censat", "rDNA": "rDNA",
         "mixedAlpha": "hor", "GAP": None}
    return m.get(base.split("_")[0] if base.startswith("active") is False else "active_hor", None)

def load_intervals(genome, data="data"):
    """chrom -> list of (start, end, class_id); SD first (lower priority), censat after."""
    iv = defaultdict(list)
    if genome == "chm13":
        for l in open(f"{data}/chm13v2.0_SD.bed"):
            f = l.split("\t"); iv[f[0]].append((int(f[1]), int(f[2]), CLASS_ID["SD"]))
        bed = f"{data}/chm13v2.0_censat_v2.1.bed"
    else:
        bed = f"{data}/hg002v1.1.cenSatv2.0.bed"
    for l in open(bed):
        if l.startswith("track") or l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        c = censat_class(f[3], genome)
        if c is None or c not in CLASS_ID: continue
        iv[f[0]].append((int(f[1]), int(f[2]), CLASS_ID[c]))
    return iv

def paint(n, intervals):
    cls = np.zeros(n, np.uint8)
    for s, e, c in intervals: cls[max(s, 0):min(e, n)] = c
    return cls

def load_run(prefix):
    meta = json.load(open(prefix + ".json"))
    return meta, prefix + ".pres.u8", prefix + ".mult.u8"

def chrom_entries(path, start, length, stride):
    """Entries for genome positions [start, start+length) at the run's stride."""
    i0 = (start + stride - 1) // stride; i1 = (start + length + stride - 1) // stride
    return np.fromfile(path, np.uint8, count=i1 - i0, offset=i0), i0, i1

def bincount_u8(a, chunk=1 << 25):
    """np.bincount over a uint8 array without materializing an int64 copy of the whole array."""
    out = np.zeros(256, np.int64)
    for i in range(0, len(a), chunk): out += np.bincount(a[i:i + chunk], minlength=256)
    return out

def window_sums(vals, e):
    """Sum of a 0/1 array over entry ranges [e[j], e[j+1]) without index arrays."""
    cs = np.empty(len(vals) + 1, np.int32); cs[0] = 0
    np.cumsum(vals, dtype=np.int32, out=cs[1:])
    return cs[e[1:]] - cs[e[:-1]]

def summarize(run_prefix, genome, query_prefix, window=10000, data="data", chrom_filter=None, suffix=""):
    out_prefix = run_prefix + suffix
    meta, pres_path, mult_path = load_run(run_prefix)
    stride = meta["stride"]; subs = [s["name"] for s in meta["subjects"]]; ns = len(subs)
    iv = load_intervals(genome, data)
    nc = len(CLASS_ORDER)
    pat_counts = np.zeros((nc, 256), np.int64)         # class x presence pattern (valid entries)
    mult_counts = np.zeros((nc, 256), np.int64)        # class x multiplicity (0 = no k-mer)
    ape_mask = sum(1 << j for j, s in enumerate(subs) if s != "hg002")
    hb = (1 << subs.index("hg002")) if "hg002" in subs else 0
    win_rows = []; per_chrom = []; cols = None
    for chrom, start, length in read_idx(query_prefix):
        if chrom_filter and chrom not in chrom_filter: continue
        pres, i0, i1 = chrom_entries(pres_path, start, length, stride)
        mult, _, _ = chrom_entries(mult_path, start, length, stride)
        n = i1 - i0
        cls_full = paint(length, iv.get(chrom, []))
        cls = cls_full if stride == 1 else cls_full[np.arange(i0, i1, dtype=np.int64) * stride - start]
        del cls_full
        valid = mult > 0
        for c in range(nc):
            mc = cls == c
            if not mc.any(): continue
            pat_counts[c] += bincount_u8(pres[mc & valid])
            mult_counts[c] += bincount_u8(mult[mc])
            del mc
        nw = (length + window - 1) // window
        e = (start + np.arange(nw + 1, dtype=np.int64) * window + stride - 1) // stride - i0
        e = np.clip(e, 0, n)
        cols = {"valid": window_sums(valid, e), "uniq": window_sums(valid & (mult == 1), e)}
        for j, s in enumerate(subs):
            cols[s] = window_sums(valid & (((pres >> j) & 1) == 1), e)
        if hb:
            cols["human_specific"] = window_sums(valid & ((pres & hb) != 0) & ((pres & ape_mask) == 0), e)
            cols["chm13_only"] = window_sums(valid & (pres == 0), e)
        cols["ape_any"] = window_sums(valid & ((pres & ape_mask) != 0), e)
        for i in range(nw):
            win_rows.append((chrom, i * window, min((i + 1) * window, length)) + tuple(int(cols[c][i]) for c in cols))
        chrom_stats = {"chrom": chrom, "length": length, "valid": int(cols["valid"].sum())}
        for s in subs: chrom_stats[s] = int(cols[s].sum())
        per_chrom.append(chrom_stats)
        print(f"  {chrom}: {length} bp, valid {chrom_stats['valid']}", file=sys.stderr)
        del pres, mult, cls, valid
    col_names = ["chrom", "start", "end"] + list(cols.keys())
    with gzip.open(out_prefix + f".win{window}.tsv.gz", "wt") as f:
        f.write("\t".join(col_names) + "\n")
        for r in win_rows: f.write("\t".join(map(str, r)) + "\n")
    pat = pat_counts; patterns = np.arange(256)
    with open(out_prefix + ".classes.tsv", "w") as f:
        hdr = ["class", "valid_positions"] + [f"shared_{s}" for s in subs] + ["human_specific", "chm13_only", "ape_any", "unique_in_query", "multi_in_query"]
        f.write("\t".join(hdr) + "\n")
        rows = []
        for ci, c in enumerate(CLASS_ORDER):
            tot = int(pat[ci].sum())
            row = [c, tot]
            for j in range(ns): row.append(int(pat[ci][(patterns >> j) & 1 == 1].sum()))
            if hb:
                row.append(int(pat[ci][((patterns & hb) != 0) & ((patterns & ape_mask) == 0)].sum()))
            else:
                row.append(0)
            row.append(int(pat[ci][0]))
            row.append(int(pat[ci][(patterns & ape_mask) != 0].sum()))
            row.append(int(mult_counts[ci][1])); row.append(int(mult_counts[ci][2:].sum()))
            rows.append(row)
        tot = ["ALL", sum(r[1] for r in rows)] + [sum(r[i] for r in rows) for i in range(2, len(hdr))]
        for r in rows + [tot]: f.write("\t".join(map(str, r)) + "\n")
    np.save(out_prefix + ".patterns.npy", pat)
    np.save(out_prefix + ".multcounts.npy", mult_counts)
    json.dump({"subjects": subs, "k": meta["k"], "stride": stride, "per_chrom": per_chrom, "classes": CLASS_ORDER}, open(out_prefix + ".summary.json", "w"), indent=1)
    return rows

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("run_prefix"); ap.add_argument("--genome", default="chm13"); ap.add_argument("--query-prefix", default="out/prep/chm13")
    ap.add_argument("--window", type=int, default=10000); ap.add_argument("--data", default="data"); ap.add_argument("--chroms", default=None); ap.add_argument("--suffix", default="")
    a = ap.parse_args()
    rows = summarize(a.run_prefix, a.genome, a.query_prefix, a.window, a.data, set(a.chroms.split(",")) if a.chroms else None, a.suffix)
    for r in rows: print("\t".join(map(str, r)))
