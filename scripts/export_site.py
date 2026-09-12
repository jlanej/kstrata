"""Export the base-resolution k=31 strata into the static-site data files under docs/data/.

Levels: w1k_<chrom>.bin.gz and w10k_<chrom>.bin.gz (uint16, one row per window, 15 columns) and
w100k.bin.gz (uint32, all chromosomes concatenated; offsets in meta.json). Columns:
  0 valid, 1 single-copy, 2..8 present in [hg002, chimp, bonobo, gorilla, borang, sorang, siamang],
  9..14 age stratum [CHM13 only, human only, to Pan, to gorilla, to orangutan, to siamang].
Spotlights: spot_<id>.bin.gz = presence bytes then multiplicity bytes for a region (per base).
Annotation: annot_<chrom>.json = censat intervals [start, end, class] and merged SD intervals.
"""
import gzip, json, os, sys
import numpy as np
sys.path.insert(0, "scripts")
from strata import read_idx, load_intervals, CLASS_ORDER, CLASS_ID

RUN = "out/runs/chm13_k31_all"; PREP = "out/prep/chm13"; OUT = "docs/data"; os.makedirs(OUT, exist_ok=True)
meta = json.load(open(RUN + ".json")); subs = [s["name"] for s in meta["subjects"]]
NODES = [("hg002",), ("chimp", "bonobo"), ("gorilla",), ("borang", "sorang"), ("siamang",)]
lut = np.zeros(256, np.uint8)  # age stratum per presence pattern
for p in range(256):
    present = {s for j, s in enumerate(subs) if (p >> j) & 1}
    a = 0
    for i, members in enumerate(NODES):
        if present & set(members): a = i + 1
    lut[p] = a
contigs = [(c, s, l) for c, s, l in read_idx(PREP) if c != "chrM"]
iv = load_intervals("chm13")

def write_gz(path, data):
    """Write gzip(data) only when the decompressed content of the existing file differs (keeps git history quiet)."""
    if os.path.exists(path):
        try:
            with gzip.open(path, "rb") as f:
                if f.read() == data: return
        except OSError: pass
    with gzip.GzipFile(path, "wb", compresslevel=6, mtime=0) as f: f.write(data)

def gz(path, arr): write_gz(path, np.ascontiguousarray(arr).tobytes())

def level_sums(col, w):
    n = len(col); nw = (n + w - 1) // w
    pad = np.zeros(nw * w, np.uint8); pad[:n] = col
    return pad.reshape(nw, w).sum(axis=1, dtype=np.int64)

w100 = []; offsets = {}; total100 = 0; chrom_meta = []
for chrom, start, length in contigs:
    pres = np.fromfile(RUN + ".pres.u8", np.uint8, count=length, offset=start)
    mult = np.fromfile(RUN + ".mult.u8", np.uint8, count=length, offset=start)
    valid = mult > 0; age = lut[pres]
    cols = [valid, valid & (mult == 1)] + [valid & (((pres >> j) & 1) == 1) for j in range(7)] + [valid & (age == a) for a in range(6)]
    k1 = np.stack([level_sums(c.view(np.uint8), 1000) for c in cols], axis=1)      # (nwin, 15)
    del cols, pres, mult, valid, age
    n1 = k1.shape[0]
    gz(f"{OUT}/w1k_{chrom}.bin.gz", k1.astype(np.uint16))
    n10 = (n1 + 9) // 10; pad = np.zeros((n10 * 10, 15), np.int64); pad[:n1] = k1
    k10 = pad.reshape(n10, 10, 15).sum(axis=1); gz(f"{OUT}/w10k_{chrom}.bin.gz", k10.astype(np.uint16))
    n100 = (n10 + 9) // 10; pad = np.zeros((n100 * 10, 15), np.int64); pad[:n10] = k10
    k100 = pad.reshape(n100, 10, 15).sum(axis=1); w100.append(k100.astype(np.uint32)); offsets[chrom] = total100; total100 += n100
    # annotation
    cen = [[int(s), int(e), int(c)] for s, e, c in iv.get(chrom, []) if c != CLASS_ID["SD"]]
    sd = sorted((int(s), int(e)) for s, e, c in iv.get(chrom, []) if c == CLASS_ID["SD"]); merged = []
    for s, e in sd:
        if merged and s <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], e)
        else: merged.append([s, e])
    json.dump({"censat": cen, "sd": merged}, open(f"{OUT}/annot_{chrom}.json", "w"), separators=(",", ":"))
    chrom_meta.append({"name": chrom, "length": int(length), "n1k": int(n1), "n10k": int(n10), "n100k": int(n100), "off100k": offsets[chrom]})
    print(chrom, n1, file=sys.stderr)
gz(f"{OUT}/w100k.bin.gz", np.concatenate(w100))

# spotlights: fixed boundary regions plus the largest interval of several satellite classes
def largest(cls, min_len=0):
    best = None
    for chrom, start, length in contigs:
        for s, e, c in iv.get(chrom, []):
            if CLASS_ORDER[c] == cls and e - s > (best[2] - best[1] if best else min_len): best = (chrom, s, e)
    return best
spots = [
    ("chr1_hor_start", "chr1: monomeric and diverged layers into the active HOR array", "chr1", 121_500_000, 122_000_000),
    ("chr1_hor_end", "chr1: end of the active array and the satellite mosaic beyond it", "chr1", 126_200_000, 126_700_000),
    ("chr1_hsat2", "chr1: into the HSat2 block", "chr1", 128_800_000, 129_200_000),
    ("chr8_layers", "chr8: monomeric layer beside the active array", "chr8", 43_600_000, 44_500_000),
    ("chr17_layers", "chr17: HSat and monomeric layers 2 Mb from the active array", "chr17", 21_800_000, 22_400_000),
    ("chr17_hor_start", "chr17: start of the active array", "chr17", 23_700_000, 24_100_000),
    ("chrX_hor_start", "chrX: start of the active array", "chrX", 57_600_000, 58_000_000),
    ("chr9_hsat3", "chr9: inside the HSat3 block", "chr9", 50_000_000, 50_300_000),
    ("chr16_hsat2", "chr16: inside the HSat2 block", "chr16", 40_000_000, 40_300_000),
    ("chr20_unique", "chr20: ordinary non-repeat sequence (the alignment-validated region)", "chr20", 20_000_000, 20_300_000),
]
for cls, label in [("hsat1B", "largest HSat1B block"), ("hsat1A", "largest HSat1A block"), ("bsat", "largest beta-satellite block"), ("rDNA", "edge of the largest rDNA model"), ("gsat", "largest gamma-satellite block")]:
    b = largest(cls)
    if b:
        chrom, s, e = b
        if cls == "rDNA": spots.append((f"{cls}_edge", f"{chrom}: {label}", chrom, max(s - 150_000, 0), s + 150_000))
        else:
            mid = (s + e) // 2; half = min(150_000, (e - s) // 2 + 50_000)
            spots.append((f"{cls}_block", f"{chrom}: {label}", chrom, max(mid - half, 0), mid + half))
# largest merged SD
best = None
for chrom, start, length in contigs:
    sd = json.load(open(f"{OUT}/annot_{chrom}.json"))["sd"]
    for s, e in sd:
        if best is None or e - s > best[2] - best[1]: best = (chrom, s, e)
if best: spots.append(("sd_largest", f"{best[0]}: largest segmental duplication block", best[0], best[1], min(best[2], best[1] + 400_000)))
spot_meta = []
for sid, desc, chrom, s, e in spots:
    off = {c: st for c, st, l in contigs}[chrom]; L = {c: l for c, st, l in contigs}[chrom]; e = min(e, L)
    pres = np.fromfile(RUN + ".pres.u8", np.uint8, count=e - s, offset=off + s)
    mult = np.fromfile(RUN + ".mult.u8", np.uint8, count=e - s, offset=off + s)
    write_gz(f"{OUT}/spot_{sid}.bin.gz", pres.tobytes() + mult.tobytes())
    spot_meta.append({"id": sid, "desc": desc, "chrom": chrom, "start": int(s), "end": int(e)})
    print("spot", sid, chrom, s, e, file=sys.stderr)

# ultraconserved: identical segments >= 1001 bp with siamang or an orangutan, from the k=1001 run if present
ucs = []
K1 = "out/runs/chm13_k1001_all"
if os.path.exists(K1 + ".json"):
    m1 = json.load(open(K1 + ".json")); s1 = [s["name"] for s in m1["subjects"]]; st = m1["stride"]; k = m1["k"]
    want = {n: s1.index(n) for n in ("siamang", "borang", "sorang", "gorilla", "chimp") if n in s1}
    for chrom, start, length in contigs:
        i0 = (start + st - 1) // st; i1 = (start + length + st - 1) // st
        pres = np.fromfile(K1 + ".pres.u8", np.uint8, count=i1 - i0, offset=i0)
        pos = np.arange(i0, i1, dtype=np.int64) * st - start
        b = ((pres >> want["siamang"]) & 1).astype(np.int8) if "siamang" in want else None
        if b is None: break
        d = np.diff(np.concatenate(([0], b, [0]))); starts = np.nonzero(d == 1)[0]; ends = np.nonzero(d == -1)[0]
        for a, e in zip(starts, ends):
            ps, pe = int(pos[a]), int(pos[e - 1]) + k
            bits = {n: bool(((pres[a:e] >> j) & 1).all()) for n, j in want.items()}
            ucs.append({"chrom": chrom, "start": ps, "end": pe, "len": pe - ps, "in": [n for n, v in bits.items() if v]})
    ucs.sort(key=lambda u: -u["len"])
    for i, u in enumerate(ucs):
        off = {c: st_ for c, st_, l in contigs}[u["chrom"]]; s = max(u["start"] - 3000, 0); e = u["end"] + 3000
        pres = np.fromfile(RUN + ".pres.u8", np.uint8, count=e - s, offset=off + s); mult = np.fromfile(RUN + ".mult.u8", np.uint8, count=e - s, offset=off + s)
        sid = f"ucs{i:02d}"
        write_gz(f"{OUT}/spot_{sid}.bin.gz", pres.tobytes() + mult.tobytes())
        u["spot"] = sid; u["spot_start"] = s; u["spot_end"] = e
json.dump(ucs, open(f"{OUT}/ucs.json", "w"), separators=(",", ":"))

def tsv(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path)]; return {r[0]: dict(zip(rows[0][1:], map(int, r[1:]))) for r in rows[1:]}
ladder = {}
for f in sorted(os.listdir("tables")):
    if f.startswith("chm13_k") and f.endswith("_all.classes.tsv"):
        k = int(f.split("_k")[1].split("_")[0]); ladder[k] = tsv("tables/" + f)
json.dump({
    "k": 31, "subjects": subs, "ages": ["CHM13 only", "human only", "to Pan", "to gorilla", "to orangutan", "to siamang"],
    "classes": CLASS_ORDER, "chroms": chrom_meta, "spots": spot_meta,
    "k31": tsv("tables/chm13_k31_all.classes.tsv"), "ladder": {str(k): v for k, v in sorted(ladder.items())},
    "arrays": json.load(open("tables/chm13_k31_all.arrays.json")), "ages_table": json.load(open("tables/chm13_k31_all.ages.json")),
    "within": tsv("tables/hg002mat_k31_patchm13.auto.classes.tsv"),
}, open(f"{OUT}/meta.json", "w"), separators=(",", ":"))
print("exported", file=sys.stderr)
