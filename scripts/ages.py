"""Age strata from the presence patterns of a run (class x 256-pattern counts saved by strata.py).

Age of a k-mer = the deepest lineage in which it is present (ignoring absences in shallower
lineages, which are losses or incomplete lineage sorting). Bits are the subject order in the run.
"""
import json, sys
import numpy as np

NODES = [  # (age label, approx Mya, set of subjects that define it)
    ("chm13_only", 0.0, set()),
    ("human_specific", 0.0, {"hg002"}),
    ("pan", 6.4, {"chimp", "bonobo"}),
    ("gorilla", 9.0, {"gorilla"}),
    ("pongo", 15.5, {"borang", "sorang"}),
    ("hylobatid", 19.5, {"siamang"}),
]

def age_of(pattern, subs):
    present = {s for j, s in enumerate(subs) if (pattern >> j) & 1}
    if pattern == 0: return "chm13_only"
    age = "human_specific" if "hg002" in present else "chm13_only"
    for label, _, members in NODES[2:]:
        if present & members: age = label
    if age == "chm13_only" and present: age = "not_in_hg002_but_in_ape"
    return age

def tables(run_prefix):
    S = json.load(open(run_prefix + ".summary.json")); subs = S["subjects"]; classes = S["classes"]
    pat = np.load(run_prefix + ".patterns.npy")  # class x 256
    pats = np.arange(256)
    ages = [age_of(p, subs) for p in pats]
    labels = [n[0] for n in NODES] + ["not_in_hg002_but_in_ape"]
    out = {}
    for ci, c in enumerate(classes):
        tot = pat[ci].sum()
        if tot == 0: continue
        row = {"positions": int(tot)}
        for lab in labels:
            row[lab] = int(pat[ci][[i for i, a in enumerate(ages) if a == lab]].sum())
        # lineage-sorting patterns among k-mers present in at least one of Pan/gorilla and absent from orangutans and siamang
        def has(p, names): return any((p >> subs.index(n)) & 1 for n in names if n in subs)
        pan_not_gor = sum(int(pat[ci][p]) for p in pats if has(p, ["chimp", "bonobo"]) and not has(p, ["gorilla"]) and not has(p, ["borang", "sorang", "siamang"]))
        gor_not_pan = sum(int(pat[ci][p]) for p in pats if has(p, ["gorilla"]) and not has(p, ["chimp", "bonobo"]) and not has(p, ["borang", "sorang", "siamang"]))
        both = sum(int(pat[ci][p]) for p in pats if has(p, ["gorilla"]) and has(p, ["chimp", "bonobo"]) and not has(p, ["borang", "sorang", "siamang"]))
        chimp_not_bonobo = sum(int(pat[ci][p]) for p in pats if has(p, ["chimp"]) and not has(p, ["bonobo"]))
        bonobo_not_chimp = sum(int(pat[ci][p]) for p in pats if has(p, ["bonobo"]) and not has(p, ["chimp"]))
        row.update({"pan_not_gorilla": pan_not_gor, "gorilla_not_pan": gor_not_pan, "pan_and_gorilla_only": both,
                    "chimp_not_bonobo": chimp_not_bonobo, "bonobo_not_chimp": bonobo_not_chimp})
        out[c] = row
    tot = {k: sum(r[k] for r in out.values()) for k in next(iter(out.values()))}
    out["ALL"] = tot
    return out, labels

if __name__ == "__main__":
    out, labels = tables(sys.argv[1])
    cols = ["positions"] + labels + ["pan_not_gorilla", "gorilla_not_pan", "pan_and_gorilla_only", "chimp_not_bonobo", "bonobo_not_chimp"]
    print("class\t" + "\t".join(cols))
    for c, r in out.items(): print(c + "\t" + "\t".join(str(r[k]) for k in cols))
    json.dump(out, open(sys.argv[1] + ".ages.json", "w"), indent=1)
