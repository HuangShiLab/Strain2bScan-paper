#!/usr/bin/env python3
"""Fig 4 (revised, labeled): impact of reference-genome completeness on strain identification across
all 15 species. Precision (A) and recall (B) vs completeness; each species a labeled coloured line;
bold median over the RESOLVABLE species. The near-clonal outlier (single 0.95 cluster) is drawn
separately (dashed) and excluded from the median, because once its references degrade the single
cluster splits into spurious singletons — a cluster-relabeling artifact, not a completeness effect."""
import os, numpy as np
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
PAPER = "/Users/macstudio/Downloads/Strain2bScan-paper"
rows = [l.rstrip("\n").split("\t") for l in open(f"{WORK}/sim/refqual15/refqual_15species.tsv")][1:]
# clusters per species (to flag near-clonal outliers = 1 cluster in the pool)
nclust = {}
for l in open(f"{PAPER}/figure_raw_data/sim_pool_summary.tsv").read().splitlines()[1:]:
    sp, ng, nc = l.split("\t"); nclust[sp] = int(nc)

P = defaultdict(dict); R = defaultdict(dict)
for r in rows:
    sp, comp = r[0], float(r[3]); P[sp][comp] = float(r[7]); R[sp][comp] = float(r[8])
species = sorted(P)
comps = sorted({float(r[3]) for r in rows}, reverse=True)
xt = [c * 100 for c in comps]

# outlier = near-clonal (single pool cluster) AND recall flat-low after any degradation
def is_outlier(sp):
    deg = [R[sp][c] for c in comps if c < 1.0]
    return nclust.get(sp, 9) <= 1 and np.mean(deg) < 0.2
outliers = [sp for sp in species if is_outlier(sp)]
resolvable = [sp for sp in species if sp not in outliers]

def med(D, sps):
    return [np.median([D[sp][c] for sp in sps if c in D[sp]]) for c in comps]

def short(sp):
    g, s = sp.split("_", 1); return f"{g[0]}. {s.replace('_',' ')}"

cmap = plt.get_cmap("tab20")
colors = {sp: cmap(i % 20) for i, sp in enumerate(species)}

fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 5.8), sharex=True)
for ax, D, lab, mc in [(a1, P, "precision", "#111"), (a2, R, "recall", "#111")]:
    for sp in resolvable:
        xs = [c * 100 for c in comps]; ys = [D[sp][c] for c in comps]
        ax.plot(xs, ys, "-", color=colors[sp], lw=1.4, alpha=0.9, label=short(sp))
    for sp in outliers:
        xs = [c * 100 for c in comps]; ys = [D[sp][c] for c in comps]
        ax.plot(xs, ys, "--", color="#333", lw=1.6, alpha=0.9, label=short(sp) + " (near-clonal)")
    ax.plot(xt, med(D, resolvable), "-o", color=mc, lw=3.2, ms=7, zorder=5,
            label=f"median ({len(resolvable)} resolvable)")
    ax.set_xlabel("reference-genome completeness (%)"); ax.set_ylabel(lab)
    ax.set_ylim(0, 1.05); ax.set_xticks(xt)
    ax.grid(alpha=0.25)
    ax.set_title(f"{lab.capitalize()} vs reference completeness", fontsize=11)
# single shared legend (species labels) to the right
h, l = a1.get_legend_handles_labels()
fig.legend(h, l, loc="center left", bbox_to_anchor=(0.995, 0.5), fontsize=8.2, frameon=True)
fig.suptitle("Fig 4. Reference-genome completeness controls strain-identification accuracy (15 species; "
             f"median over {len(resolvable)} resolvable)", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 0.83, 0.95])
fig.savefig(f"{PAPER}/figures/refqual_figure.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{PAPER}/figures/refqual_figure.pdf", bbox_inches="tight")
print("outliers (reported separately):", outliers)
print("median precision (resolvable):", [round(x, 2) for x in med(P, resolvable)])
print("median recall (resolvable):   ", [round(x, 2) for x in med(R, resolvable)])
print("wrote figures/refqual_figure.png/pdf")
