#!/usr/bin/env python3
"""Fig 4 (problem -> solution): reference-genome completeness vs strain-identification accuracy across
15 species, and how the --containment clustering mode addresses it. Default Jaccard clustering
(dashed) loses precision/recall as references degrade because incomplete genomes split from their
complete relatives; max-containment clustering (solid) keeps them together, restoring accuracy down to
~80% completeness. Median over the 14 resolvable species (near-clonal M. tuberculosis reported inset)."""
import os, numpy as np
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

PAPER = "/Users/macstudio/Downloads/Strain2bScan-paper"
def load(p):
    P = defaultdict(dict); R = defaultdict(dict)
    for l in open(p):
        f = l.rstrip("\n").split("\t")
        if l.startswith("#") or f[0] == "species": continue
        P[f[0]][float(f[3])] = float(f[7]); R[f[0]][float(f[3])] = float(f[8])
    return P, R
Pj, Rj = load(f"{PAPER}/results/refqual_15species.tsv")
Pc, Rc = load(f"{PAPER}/results/refqual_15species_containment.tsv")
comps = sorted(set(Pj[next(iter(Pj))]), reverse=True)
xt = [c * 100 for c in comps]
OUT = "Mycobacterium_tuberculosis"
res = [s for s in Pj if s != OUT]
def med(D, sps): return [np.median([D[s][c] for s in sps if c in D[s]]) for c in comps]

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 5.4), sharex=True)
for ax, Dj, Dc, lab, col in [(a1, Pj, Pc, "precision", "#1f77b4"), (a2, Rj, Rc, "recall", "#2ca02c")]:
    for s in res:  # light per-species containment lines (spread)
        ax.plot(xt, [Dc[s][c] for c in comps], "-", color=col, lw=0.7, alpha=0.22, zorder=1)
    ax.plot(xt, med(Dj, res), "--o", color="#888", lw=2.4, ms=6, zorder=4, label="default (Jaccard)")
    ax.plot(xt, med(Dc, res), "-o", color=col, lw=3.0, ms=7, zorder=5, label="--containment")
    ax.fill_between(xt, med(Dj, res), med(Dc, res), color=col, alpha=0.12, zorder=0)
    ax.set_xlabel("reference-genome completeness (%)"); ax.set_ylabel(f"median {lab} (14 species)")
    ax.set_ylim(0, 1.05); ax.set_xticks(xt); ax.grid(alpha=0.25)
    ax.legend(loc="lower right", fontsize=10)
    ax.set_title(f"{lab.capitalize()}: containment recovers accuracy\nunder incomplete references", fontsize=11)
# inset: the M. tuberculosis near-clonal artifact fix (recall)
axi = a2.inset_axes([0.10, 0.12, 0.42, 0.42])
axi.plot(xt, [Rj[OUT][c] for c in comps], "--o", color="#888", lw=1.5, ms=3)
axi.plot(xt, [Rc[OUT][c] for c in comps], "-o", color="#d62728", lw=2, ms=3)
axi.set_ylim(0, 1.08); axi.set_xticks([50, 70, 90, 100]); axi.tick_params(labelsize=6.5)
axi.set_title("M. tuberculosis recall\n(near-clonal; artifact fixed)", fontsize=6.8)
axi.grid(alpha=0.25)
fig.suptitle("Fig 4. Reference completeness limits strain ID under Jaccard clustering; the --containment "
             "mode restores it (15 species)", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(f"{PAPER}/figures/refqual_figure.png", dpi=150)
fig.savefig(f"{PAPER}/figures/refqual_figure.pdf")
print("Jaccard  precision:", [round(x,2) for x in med(Pj,res)], "recall:", [round(x,2) for x in med(Rj,res)])
print("Contain  precision:", [round(x,2) for x in med(Pc,res)], "recall:", [round(x,2) for x in med(Rc,res)])
print("wrote figures/refqual_figure.png/pdf")
