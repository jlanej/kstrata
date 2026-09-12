"""Independent check of the chimpanzee presence bits on 5 Mb of CHM13 chr20 (20-25 Mb).

minimap2 (asm20, --eqx) aligns the region to the chimp assembly. A query position p is
'exact by alignment' when some alignment covers [p, p+31) with '=' operations only and no indel
boundary inside. Every such 31-mer exists in chimp, so kstrata must report it present
(soundness). Positions whose best alignment carries a mismatch/indel inside the window may
still be present elsewhere in chimp (repeats), so their presence rate is reported, not asserted.
"""
import json, re, sys
import numpy as np

K = 31
REGION_CHROM, REGION_START, REGION_LEN = "chr20", 20_000_000, 5_000_000  # 0-based start
paf = sys.argv[1]; run = sys.argv[2]; idx = sys.argv[3]

exact_any = np.zeros(REGION_LEN, bool)     # base is '=' in some alignment
brk_any = np.ones(REGION_LEN, bool)        # provisional: position p can start an exact 31-mer in some alignment
can_start = np.zeros(REGION_LEN, bool)
prim_mismatch = np.zeros(REGION_LEN, bool) # in the primary alignment, the 31-mer window at p contains a mismatch/indel
prim_cov = np.zeros(REGION_LEN, bool)
n_aln = 0
for l in open(paf):
    f = l.rstrip("\n").split("\t")
    qs, qe, strand = int(f[2]), int(f[3]), f[4]
    tags = {t[:4]: t[5:] for t in f[12:]}
    cg = tags.get("cg:Z"); tp = tags.get("tp:A", "P")
    if cg is None: continue
    n_aln += 1
    ops = re.findall(r"(\d+)([=XIDM])", cg)
    # walk in query coordinates; for '-' strand minimap2 reports the CIGAR along the query in the reversed orientation,
    # so we build per-alignment arrays in target-walk order and then flip to query orientation.
    L = qe - qs
    eq = np.zeros(L, bool); brk = np.zeros(L, bool)  # brk[i]: an indel sits between query base i and i+1
    i = 0
    for n, op in ops:
        n = int(n)
        if op == "=": eq[i:i + n] = True; i += n
        elif op in "XM": i += n
        elif op == "I": i += n              # inserted query bases: not '='
        elif op == "D":
            if 0 < i <= L: brk[i - 1] = True  # deletion between query base i-1 and i
    assert i == L, (i, L)
    if strand == "-": eq = eq[::-1]; brk = np.concatenate(([False], brk[:-1][::-1])) if L > 1 else brk
    # a 31-mer at query position p (region coords qs+p) is exact when eq[p..p+30] all True and no brk in [p, p+29]
    cs_eq = np.concatenate(([0], np.cumsum(~eq))); cs_brk = np.concatenate(([0], np.cumsum(brk)))
    p = np.arange(0, L - K + 1)
    ok = (cs_eq[p + K] - cs_eq[p] == 0) & (cs_brk[p + K - 1] - cs_brk[p] == 0)
    can_start[qs + p[ok]] = True
    if tp == "P":
        prim_cov[qs:qe] = True
        prim_mismatch[qs + p[~ok]] = True
        prim_mismatch[qs + p[ok]] = False

meta = json.load(open(run + ".json")); subs = [s["name"] for s in meta["subjects"]]; bit = 1 << subs.index("chimp")
off = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(idx)}[REGION_CHROM] + REGION_START
pres = np.fromfile(run + ".pres.u8", np.uint8, count=REGION_LEN, offset=off)
mult = np.fromfile(run + ".mult.u8", np.uint8, count=REGION_LEN, offset=off)
present = (pres & bit) != 0
valid = mult > 0
print(f"alignments: {n_aln}; region positions: {REGION_LEN}; with a valid 31-mer: {valid.sum()}")
a = can_start & valid
print(f"exact-by-alignment 31-mers: {a.sum()} ({100*a.sum()/valid.sum():.2f}% of valid); kstrata present among them: {present[a].sum()} ({100*present[a].mean():.4f}%)")
b = prim_cov & prim_mismatch & ~can_start & valid
print(f"primary alignment has a mismatch/indel inside the window and no other alignment is exact: {b.sum()}; kstrata present: {present[b].sum()} ({100*present[b].mean():.2f}%)")
bu = b & (mult == 1)
print(f"   of those with a single-copy 31-mer in CHM13: {bu.sum()}; present in chimp: {present[bu].sum()} ({100*present[bu].mean():.2f}%)")
c = ~prim_cov & valid
print(f"not covered by a primary alignment: {c.sum()}; kstrata present: {present[c].sum()} ({100*present[c].mean():.2f}%)")
print(f"overall: alignment-exact {100*a.sum()/valid.sum():.2f}% vs kstrata-present {100*present[valid].mean():.2f}%")
