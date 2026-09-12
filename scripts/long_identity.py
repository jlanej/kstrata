"""From a long-k run (stride 8): count the identical segments of length >= k between CHM13 and each
subject (runs of consecutive sampled positions whose k-mer is present), and their total span, by class."""
import json, sys
import numpy as np
sys.path.insert(0, "scripts")
from strata import read_idx, load_intervals, paint, CLASS_ORDER

run = sys.argv[1]
meta = json.load(open(run + ".json")); subs = [s["name"] for s in meta["subjects"]]; k = meta["k"]; stride = meta["stride"]
iv = load_intervals("chm13")
out = {s: {"segments": 0, "positions": 0, "by_class": {}} for s in subs}
for chrom, start, length in read_idx("out/prep/chm13"):
    i0 = (start + stride - 1) // stride; i1 = (start + length + stride - 1) // stride
    pres = np.fromfile(run + ".pres.u8", np.uint8, count=i1 - i0, offset=i0)
    pos = np.arange(i0, i1, dtype=np.int64) * stride - start
    cls = paint(length, iv.get(chrom, []))[pos]
    for j, s in enumerate(subs):
        b = ((pres >> j) & 1).astype(np.int8)
        starts = np.nonzero(np.diff(np.concatenate(([0], b))) == 1)[0]
        out[s]["segments"] += len(starts); out[s]["positions"] += int(b.sum())
        for c in np.unique(cls[starts]):
            out[s]["by_class"][CLASS_ORDER[c]] = out[s]["by_class"].get(CLASS_ORDER[c], 0) + int((cls[starts] == c).sum())
print(f"k={k}: segments of >= {k} identical bases (runs of sampled start positions, stride {stride}) and sampled start positions")
for s in subs:
    o = out[s]
    top = sorted(o["by_class"].items(), key=lambda x: -x[1])[:4]
    print(f"{s:8s} segments={o['segments']:8d} start_positions~={o['positions']*stride:10d} classes: " + ", ".join(f"{c}={n}" for c, n in top))
json.dump(out, open(run + ".segments.json", "w"), indent=1)
