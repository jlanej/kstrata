"""Brute-force validation of kstrata on synthetic genomes.

Builds a query genome with planted repeats (forward and reverse-complement copies, N runs,
several contigs) and two subject genomes that share planted segments with it, runs kstrata for
several k (both the exact k<=32 path and the hash k>32 path, with stride 1 and stride 4, with
one and several partitions), and compares every output byte with Python sets of canonical
k-mers.
"""
import json, os, random, subprocess, sys, tempfile
import numpy as np

BIN = sys.argv[1]
rng = random.Random(7)
COMP = str.maketrans("ACGT", "TGCA")

def rc(s): return s.translate(COMP)[::-1]
def rand(n): return "".join(rng.choice("ACGT") for _ in range(n))

def make_query():
    a = rand(20000)
    rep = rand(400)
    b = a[:5000] + rep + a[5000:9000] + rc(rep) + "N" * 37 + a[9000:14000] + rep[:200] + a[14000:]
    c = rand(3000) + rep + rand(500)          # third contig, another copy
    pal = rand(60); pal = pal + rc(pal)       # perfect palindrome
    d = rand(1000) + pal + rand(1000)
    return {"c1": b, "c2": c, "c3": d}

def make_subjects(q):
    s1 = {"s1a": rand(4000) + q["c1"][1000:3000] + rand(2000) + rc(q["c1"][10000:12000]) + rand(1000),
          "s1b": rand(500) + q["c2"][3000:3400] + rand(500)}
    s2 = {"s2a": rand(6000) + q["c3"][900:1300] + rand(100) + q["c1"][6000:6020] + rand(3000)}
    return s1, s2

def write_fa(path, d):
    with open(path, "w") as f:
        for n, s in d.items():
            f.write(f">{n} desc\n")
            for i in range(0, len(s), 60): f.write(s[i:i+60] + "\n")

def kmers_of(d, k):
    out = {}
    for s in d.values():
        for i in range(len(s) - k + 1):
            w = s[i:i+k]
            if "N" in w: continue
            c = min(w, rc(w))
            out[c] = out.get(c, 0) + 1
    return out

def expected(q, subs, k, stride):
    qk = kmers_of(q, k)
    sk = [set(kmers_of(s, k)) for s in subs]
    total = sum(len(s) for s in q.values())
    n = (total + stride - 1) // stride
    pres = np.zeros(n, np.uint8); mult = np.zeros(n, np.uint8)
    off = 0
    for s in q.values():
        for i in range(len(s)):
            g = off + i
            if g % stride: continue
            if i + k > len(s): continue
            w = s[i:i+k]
            if "N" in w: continue
            c = min(w, rc(w))
            mult[g // stride] = min(qk[c], 255)
            bits = 0
            for j, ss in enumerate(sk):
                if c in ss: bits |= 1 << j
            pres[g // stride] = bits
        off += len(s)
    return pres, mult

def main():
    q = make_query(); s1, s2 = make_subjects(q)
    with tempfile.TemporaryDirectory() as td:
        for name, d in [("q", q), ("s1", s1), ("s2", s2)]:
            write_fa(f"{td}/{name}.fa", d)
            subprocess.run([BIN, "prep", f"{td}/{name}.fa", f"{td}/{name}"], check=True)
        fails = 0
        for k, stride, parts in [(21, 1, 1), (21, 1, 4), (32, 1, 2), (16, 4, 2), (33, 1, 1), (51, 1, 4), (51, 4, 2), (201, 1, 2), (201, 8, 4), (25, 1, "small"), (101, 2, "small")]:
            out = f"{td}/o_{k}_{stride}_{parts}"
            extra = ["--small-query"] if parts == "small" else ["--query-parts", str(parts)]
            subprocess.run([BIN, "run", "--query", f"{td}/q", "--query-name", "q", "--subject", f"s1={td}/s1", "--subject", f"s2={td}/s2",
                            "-k", str(k), "--stride", str(stride), "--budget-gb", "0.0001", "--chunk-mb", "0.005", "--out", out] + extra,
                           check=True, stderr=subprocess.DEVNULL)
            pres = np.fromfile(out + ".pres.u8", np.uint8); mult = np.fromfile(out + ".mult.u8", np.uint8)
            ep, em = expected(q, [s1, s2], k, stride)
            meta = json.load(open(out + ".json"))
            ok = np.array_equal(pres, ep) and np.array_equal(mult, em)
            qk = kmers_of(q, k)
            ok2 = meta["query_distinct"] == len(qk)
            hist = dict((int(a), int(b)) for a, b in meta["self_histogram"])
            exp_hist = {}
            for c in qk.values(): exp_hist[c] = exp_hist.get(c, 0) + 1
            ok3 = hist == exp_hist
            print(f"k={k:4d} stride={stride} parts={str(parts):5s} mode={meta['mode']:7s} arrays={'OK' if ok else 'FAIL'} distinct={'OK' if ok2 else 'FAIL'} hist={'OK' if ok3 else 'FAIL'} "
                  f"shared={[ (s['name'], s['shared_entries']) for s in meta['subjects']]} nonzero_pres={int((pres>0).sum())} mult>1={int((mult>1).sum())}")
            if not (ok and ok2 and ok3):
                fails += 1
                bad = np.nonzero(pres != ep)[0][:10]; print("  pres mismatches at", bad, pres[bad], ep[bad])
                bad = np.nonzero(mult != em)[0][:10]; print("  mult mismatches at", bad, mult[bad], em[bad])
        print("ALL OK" if fails == 0 else f"{fails} FAILURES"); sys.exit(fails)

main()
