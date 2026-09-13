"""Figure: sharing and age strata of CHM13 31-mers by their multiplicity in CHM13 (from tables/chm13_k31_by_multiplicity.json)."""
import json, sys
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
D = json.load(open("tables/chm13_k31_by_multiplicity.json")); bins = D["bins"]; subs = D["subjects"]
AGE_COL = ["#bbbbbb", "#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4"]; AGES = ["CHM13 only", "human only", "to Pan", "to gorilla", "to orangutan", "to siamang"]
SUBJ_COL = {"hg002": "#1f77b4", "chimp": "#ff7f0e", "gorilla": "#2ca02c", "sorang": "#d62728", "siamang": "#9467bd"}
SUBJ_LAB = {"hg002": "HG002 (human)", "chimp": "chimpanzee", "gorilla": "gorilla", "sorang": "Sumatran orangutan", "siamang": "siamang"}
classes = [("ALL", "whole genome"), ("unique", "non-repeat, non-SD"), ("SD", "segmental duplication"), ("hor_active", "active HOR"), ("hsat3", "HSat3")]
fig, axes = plt.subplots(2, len(classes), figsize=(3.1 * len(classes), 6.4), sharey="row")
for ci, (c, lab) in enumerate(classes):
    rows = D["classes"][c]; tot = sum(r["positions"] for r in rows); x = np.arange(len(bins))
    ax = axes[0, ci]
    for s in ["hg002", "chimp", "gorilla", "sorang", "siamang"]:
        ax.plot(x, [100 * r["shared_" + s] / r["positions"] if r["positions"] else np.nan for r in rows], marker="o", ms=4, color=SUBJ_COL[s], label=SUBJ_LAB[s])
    ax.set_title(lab, fontsize=9); ax.set_xticks(x); ax.set_xticklabels(bins, fontsize=8); ax.grid(alpha=0.3); ax.set_ylim(0, 102)
    if ci == 0: ax.set_ylabel("% of 31-mers present in the genome", fontsize=9); ax.legend(fontsize=7)
    ax2 = axes[1, ci]; bottom = np.zeros(len(bins))
    for a in range(6):
        vals = np.array([100 * r["age%d" % a] / r["positions"] if r["positions"] else 0 for r in rows]); ax2.bar(x, vals, bottom=bottom, color=AGE_COL[a], label=AGES[a], width=0.8); bottom += vals
    share = [100 * r["positions"] / tot for r in rows]
    for i, sh in enumerate(share): ax2.text(i, 101, f"{sh:.0f}%", ha="center", fontsize=7, color="#555")
    ax2.set_xticks(x); ax2.set_xticklabels(bins, fontsize=8); ax2.set_ylim(0, 108); ax2.set_xlabel("copies of the 31-mer in CHM13", fontsize=9)
    if ci == 0: ax2.set_ylabel("age strata (% of the bin's 31-mers)", fontsize=9); ax2.legend(fontsize=7, loc="lower left")
fig.suptitle("Sharing and age of CHM13 31-mers by their copy number in CHM13 (numbers above the bars: share of the class's positions in each bin)", fontsize=10)
fig.tight_layout(); fig.savefig("fig/fig_multiplicity.png", dpi=150); print("figure written")
