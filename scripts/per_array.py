"""Per-centromere table: for every active HOR array of CHM13 (censat v2.1 labels ending in L),
the fraction of its 31-mers present in each genome, computed exactly from the base-resolution arrays."""
import json, sys
import numpy as np
sys.path.insert(0, "scripts")
from strata import read_idx, load_intervals, CLASS_ID

run = sys.argv[1] if len(sys.argv) > 1 else "out/runs/chm13_k31_all"
meta = json.load(open(run + ".json")); subs = [s["name"] for s in meta["subjects"]]
off = {c: s for c, s, l in read_idx("out/prep/chm13")}
iv = load_intervals("chm13")
order = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
print("chrom\tarray (Mb)\tsize (Mb)\tsingle-copy %\t" + "\t".join(f"{s} %" for s in subs) + "\tclosest ape")
rows = []
for ch in order:
    act = sorted((s, e) for s, e, c in iv.get(ch, []) if c == CLASS_ID["hor_active"])
    if not act: continue
    # merge adjacent active intervals (gaps < 500 kb) into arrays
    arrays = []
    for s, e in act:
        if arrays and s - arrays[-1][1] < 500_000: arrays[-1][1] = max(arrays[-1][1], e)
        else: arrays.append([s, e])
    for s, e in arrays:
        pres = np.concatenate([np.fromfile(run + ".pres.u8", np.uint8, count=b - a, offset=off[ch] + a) for a, b in act if a >= s and b <= e])
        mult = np.concatenate([np.fromfile(run + ".mult.u8", np.uint8, count=b - a, offset=off[ch] + a) for a, b in act if a >= s and b <= e])
        v = mult > 0; n = v.sum()
        fr = [100 * ((pres[v] >> j) & 1).mean() for j in range(len(subs))]
        apes = {sub: f for sub, f in zip(subs, fr) if sub != "hg002"}
        best = max(apes, key=apes.get)
        rows.append((ch, s, e, n, 100 * (mult[v] == 1).mean(), fr, best))
        print(f"{ch}\t{s/1e6:.2f}-{e/1e6:.2f}\t{n/1e6:.2f}\t{100*(mult[v]==1).mean():.1f}\t" + "\t".join(f"{f:.1f}" for f in fr) + f"\t{best}")
json.dump([{"chrom": r[0], "start": int(r[1]), "end": int(r[2]), "n": int(r[3]), "single_copy_pct": r[4], "shared_pct": dict(zip(subs, r[5])), "closest_ape": r[6]} for r in rows], open(run + ".arrays.json", "w"), indent=1)
