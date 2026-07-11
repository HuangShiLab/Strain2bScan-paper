#!/usr/bin/env python3
"""Fig 4 (revised): impact of reference-genome completeness on strain identification across all 15
species. Precision (A) and recall (B) vs reference completeness; one thin line per species + bold
median. Data: sim/refqual15/refqual_15species.tsv."""
import os, numpy as np
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
PAPER = "/Users/macstudio/Downloads/Strain2bScan-paper"
rows = [l.rstrip("\n").split("\t") for l in open(f"{WORK}/sim/refqual15/refqual_15species.tsv")][1:]
P = defaultdict(dict); R = defaultdict(dict)
for r in rows:
    sp, comp, prec, rec = r[0], float(r[3]), float(r[7]), float(r[8])
    P[sp][comp] = prec; R[sp][comp] = rec
species = sorted(P)
comps = sorted({float(r[3]) for r in rows}, reverse=True)   # 1.0 -> 0.5
xt = [c * 100 for c in comps]

def med(D):
    return [np.median([D[sp][c] for sp in species if c in D[sp]]) for c in comps]

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5.2), sharex=True)
for ax, D, lab in [(a1, P, "precision"), (a2, R, "recall")]:
    for sp in species:
        xs = [c * 100 for c in comps if c in D[sp]]
        ys = [D[sp][c] for c in comps if c in D[sp]]
        ax.plot(xs, ys, "-", color="#bbbbbb", lw=0.9, alpha=0.8, zorder=1)
    ax.plot(xt, med(D), "-o", color="#1f77b4" if lab == "precision" else "#2ca02c", lw=2.5, ms=6,
            label=f"median {lab}", zorder=3)
    ax.set_xlabel("reference-genome completeness (%)"); ax.set_ylabel(lab)
    ax.set_ylim(0, 1.05); ax.set_xticks(xt); ax.invert_xaxis()
    ax.grid(alpha=0.25); ax.legend(loc="lower left", fontsize=10)
    ax.set_title(f"{lab.capitalize()} vs reference completeness\n(15 species; grey = per species, bold = median)", fontsize=11)
# annotate contamination/fragmentation co-vary
a1.annotate("100% = complete;\nlower % also adds\ncontamination + fragmentation",
            (0.97, 0.06), xycoords="axes fraction", ha="right", fontsize=7.5, color="#555")
fig.suptitle("Fig 4. Reference-genome completeness controls strain-identification accuracy (all 15 species)", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(f"{PAPER}/figures/refqual_figure.png", dpi=150)
fig.savefig(f"{PAPER}/figures/refqual_figure.pdf")
# summary numbers
print("completeness ladder:", xt)
print("median precision:", [round(x, 2) for x in med(P)])
print("median recall:   ", [round(x, 2) for x in med(R)])
print("wrote figures/refqual_figure.png/pdf")
