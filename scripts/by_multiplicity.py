"""Do the strata depend on how many copies the k-mer has in CHM13? Counts of (class, multiplicity bin, presence
pattern) over every base of the k=31 run, then sharing per genome and age composition per multiplicity bin."""
import json, sys
import numpy as np
sys.path.insert(0, "scripts")
from strata import read_idx, load_intervals, paint, CLASS_ORDER, bincount_u8
RUN = "out/runs/chm13_k31_all"; PREP = "out/prep/chm13"
meta = json.load(open(RUN + ".json")); subs = [s["name"] for s in meta["subjects"]]
BINS = [(1, 1, "1"), (2, 2, "2"), (3, 9, "3-9"), (10, 99, "10-99"), (100, 254, "100-254"), (255, 255, ">=255")]
mbin = np.zeros(256, np.int64); mbin[0] = -1
for i, (lo, hi, _) in enumerate(BINS): mbin[lo:hi + 1] = i
iv = load_intervals("chm13"); nc = len(CLASS_ORDER)
counts = np.zeros((nc, len(BINS), 256), np.int64)
for chrom, start, length in read_idx(PREP):
    if chrom == "chrM": continue
    pres = np.fromfile(RUN + ".pres.u8", np.uint8, count=length, offset=start)
    mult = np.fromfile(RUN + ".mult.u8", np.uint8, count=length, offset=start)
    cls = paint(length, iv.get(chrom, []))
    for c in range(nc):
        mc = cls == c
        if not mc.any(): continue
        for b in range(len(BINS)):
            sel = mc & (mbin[mult] == b)
            if sel.any(): counts[c, b] += bincount_u8(pres[sel])
    print(chrom, file=sys.stderr)
NODES = [("hg002",), ("chimp", "bonobo"), ("gorilla",), ("borang", "sorang"), ("siamang",)]
age = np.zeros(256, np.int64)
for p in range(256):
    present = {s for j, s in enumerate(subs) if (p >> j) & 1}; a = 0
    for i, m in enumerate(NODES):
        if present & set(m): a = i + 1
    age[p] = a
pats = np.arange(256)
out = {"bins": [b[2] for b in BINS], "subjects": subs, "classes": {}}
def summarize(cnt):  # cnt: (bins, 256)
    rows = []
    for b in range(len(BINS)):
        tot = int(cnt[b].sum()); row = {"positions": tot}
        for j, s in enumerate(subs): row["shared_" + s] = int(cnt[b][((pats >> j) & 1) == 1].sum())
        for a in range(6): row["age%d" % a] = int(cnt[b][age == a].sum())
        rows.append(row)
    return rows
out["classes"]["ALL"] = summarize(counts.sum(axis=0))
for c in range(nc):
    if counts[c].sum(): out["classes"][CLASS_ORDER[c]] = summarize(counts[c])
json.dump(out, open("tables/chm13_k31_by_multiplicity.json", "w"), indent=1)
AGES = ["CHM13 only", "human only", "to Pan", "to gorilla", "to orangutan", "to siamang"]
for cname in ["ALL", "unique", "SD", "hor_active", "hsat3"]:
    print(f"\n{cname}: multiplicity bin | share of class | %chimp %gorilla %sorang %siamang | age: " + " ".join(AGES))
    rows = out["classes"][cname]; tot = sum(r["positions"] for r in rows)
    for (lo, hi, lab), r in zip(BINS, rows):
        n = r["positions"]
        if not n: continue
        print(f"  {lab:8s} {100*n/tot:5.1f}% | " + " ".join(f"{100*r['shared_'+s]/n:5.1f}" for s in ["chimp", "gorilla", "sorang", "siamang"]) + " | " + " ".join(f"{100*r['age%d'%a]/n:5.1f}" for a in range(6)))
