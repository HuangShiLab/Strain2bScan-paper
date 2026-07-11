#!/usr/bin/env python3
"""Supplementary: per-species RANK-RANK scatter of between-strain distance — whole-genome (x) vs
2bRAD and 16S (y). Since the metric is Spearman (rank) correlation, ranks are plotted directly: a
perfect monotonic relationship lies on the diagonal. 2bRAD hugs the diagonal in every species;
16S is a diffuse off-diagonal cloud (cannot order strains) — most extreme in Phocaeicola dorei and
M. tuberculosis, where 16S carries essentially no strain signal."""
import os, glob, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
PAPER = os.path.expanduser("~/Downloads/Strain2bScan-paper")
PD = f"{WORK}/results/pairdist"

# Spearman values (+n) from the recomputed table
stat = {}
for l in open(f"{PAPER}/results/mash_2brad_vs_16s_recomputed.tsv"):
    if l.startswith("#") or l.startswith("species"): continue
    f = l.split("\t"); stat[f[0]] = (int(f[1]), float(f[3]), float(f[6]))

def norm(v):
    """average-rank -> [0,1] percentile (ties get mean rank), matching the Spearman computation."""
    v = np.asarray(v, float); n = len(v)
    order = np.argsort(v, kind="mergesort"); r = np.empty(n); i = 0
    while i < n:
        j = i
        while j + 1 < n and v[order[j + 1]] == v[order[i]]: j += 1
        r[order[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r / (n - 1) if n > 1 else r * 0

files = sorted(glob.glob(f"{PD}/*.tsv"))
order = sorted(files, key=lambda p: stat.get(os.path.basename(p)[:-4], (0, 0, 1))[2])  # worst 16S first
ncol, nrow = 3, 5
fig, axes = plt.subplots(nrow, ncol, figsize=(12.5, 16))
axes = axes.ravel()
for k, p in enumerate(order):
    sp = os.path.basename(p)[:-4]
    rows = [l.rstrip("\n").split("\t") for l in open(p)][1:]
    if not rows: axes[k].axis("off"); continue
    wgs = norm([float(r[2]) for r in rows])
    b2  = norm([float(r[3]) for r in rows])
    s16 = norm([float(r[4]) for r in rows])
    ax = axes[k]
    ax.scatter(wgs, b2, s=9, color="#1f77b4", alpha=0.35, edgecolors="none", label="2bRAD")
    ax.scatter(wgs, s16, s=9, color="#d62728", alpha=0.35, edgecolors="none", label="16S")
    ax.plot([0, 1], [0, 1], color="#999", lw=0.8, ls="--")
    n, r2, r16 = stat.get(sp, (len(rows), float("nan"), float("nan")))
    g, s = sp.split("_", 1)
    ax.set_title(f"$\\it{{{g}}}$ $\\it{{{s}}}$  (n={n})\n"
                 f"ρ 2bRAD={r2:.2f}, 16S={r16:.2f}", fontsize=9.5)
    ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
    ax.tick_params(labelsize=7.5); ax.grid(alpha=0.2)
    if k == 0: ax.legend(fontsize=8, loc="upper left", markerscale=1.5, framealpha=0.9)
for k in range(len(order), len(axes)): axes[k].axis("off")
fig.supxlabel("whole-genome (Mash) between-strain distance  —  rank percentile", fontsize=12, y=0.005)
fig.supylabel("marker-space between-strain distance  —  rank percentile", fontsize=12, x=0.005)
fig.suptitle("Rank–rank between-strain distance vs whole-genome: 2bRAD tracks the diagonal, 16S does not\n"
             "(15 species, complete/near-complete genomes; sorted worst→best 16S; ρ = Spearman)", fontsize=13, y=0.998)
fig.tight_layout(rect=[0.02, 0.01, 1, 0.985])
fig.savefig(f"{PAPER}/figures/mash_2brad_vs_16s_scatter.png", dpi=140)
fig.savefig(f"{PAPER}/figures/mash_2brad_vs_16s_scatter.pdf")
print("wrote figures/mash_2brad_vs_16s_scatter.png + .pdf")
