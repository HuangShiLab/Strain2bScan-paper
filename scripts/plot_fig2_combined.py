#!/usr/bin/env python3
"""Fig 2 (combined): 2bRAD captures strain signal 16S cannot.
LEFT panel  (A) — per-species Spearman of between-strain distance vs whole-genome, 2bRAD vs 16S,
                  complete/near-complete genomes, 95% CI error bars, genome count per species.
RIGHT panel (B) — 3x5 matrix of per-species RANK-RANK scatters (whole-genome vs 2bRAD/16S); 2bRAD
                  hugs the diagonal, 16S is flat/off-diagonal. Merges the former Fig 2 + Fig S1."""
import os, glob, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

PAPER = os.path.expanduser("~/Downloads/Strain2bScan-paper")
PD = f"{PAPER}/results/pairdist"
TSV = f"{PAPER}/results/mash_2brad_vs_16s_recomputed.tsv"

# ---- load summary stats ----
rows = [l.rstrip("\n").split("\t") for l in open(TSV) if l.strip()
        and not l.startswith("#") and not l.startswith("species")]
stat = {}   # species -> (n, sp2, lo2, hi2, sp16, lo16, hi16)
for r in rows:
    stat[r[0]] = (int(r[1]), float(r[3]), float(r[4]), float(r[5]), float(r[6]), float(r[7]), float(r[8]))

def rankpct(v):
    v = np.asarray(v, float); n = len(v)
    order = np.argsort(v, kind="mergesort"); rk = np.empty(n); i = 0
    while i < n:
        j = i
        while j + 1 < n and v[order[j + 1]] == v[order[i]]: j += 1
        rk[order[i:j + 1]] = (i + j) / 2.0; i = j + 1
    return rk / (n - 1) if n > 1 else rk * 0

def gs(sp):  # "Genus species" (spaces)
    g, s = sp.split("_", 1); return g, s.replace("_", " ")

fig = plt.figure(figsize=(19, 8.8))
outer = fig.add_gridspec(1, 2, width_ratios=[0.33, 0.67], wspace=0.13)

# ================= LEFT: summary bar chart =================
axb = fig.add_subplot(outer[0, 0])
d = [(sp, *stat[sp]) for sp in stat]
d.sort(key=lambda x: x[2] - x[5])                      # by 2bRAD advantage
labels = [f"{gs(x[0])[0]}\n{gs(x[0])[1]}" for x in d]
ns = [x[1] for x in d]
s2 = np.array([x[2] for x in d]); e2 = np.array([[x[2]-x[3] for x in d], [x[4]-x[2] for x in d]])
s16 = np.array([x[5] for x in d]); e16 = np.array([[x[5]-x[6] for x in d], [x[7]-x[5] for x in d]])
y = np.arange(len(d)); h = 0.38
axb.barh(y + h/2, s2, h, color="#1f77b4", label="2bRAD", edgecolor="k", linewidth=0.3,
         xerr=e2, error_kw=dict(ecolor="#0b3d66", elinewidth=1.0, capsize=2))
axb.barh(y - h/2, s16, h, color="#d62728", label="16S rRNA", edgecolor="k", linewidth=0.3,
         xerr=e16, error_kw=dict(ecolor="#7a1416", elinewidth=1.0, capsize=2))
axb.set_yticks(y); axb.set_yticklabels(labels, fontsize=8, style="italic")
axb.set_xlabel("Spearman ρ (between-strain vs whole-genome distance)", fontsize=10)
axb.set_xlim(-0.25, 1.06); axb.set_ylim(-0.7, len(d)-0.3); axb.axvline(0, color="#888", lw=0.8)
med2, med16 = np.median(s2), np.median(s16)
axb.axvline(med2, color="#1f77b4", ls="--", lw=1.1, alpha=0.8)
axb.axvline(med16, color="#d62728", ls="--", lw=1.1, alpha=0.8)
for yi, n in zip(y, ns): axb.text(1.05, yi, f"n={n}", fontsize=7, va="center", ha="right", color="#444")
axb.text(med2, len(d)-0.35, f"2bRAD med {med2:.2f}", color="#1f77b4", fontsize=8, ha="center", va="bottom")
axb.text(med16, -0.62, f"16S med {med16:.2f}", color="#d62728", fontsize=8, ha="center", va="bottom")
axb.legend(loc="lower left", fontsize=9, framealpha=0.95); axb.grid(alpha=0.25, axis="x")
axb.set_title("Complete/near-complete genomes; error bars 95% CI", fontsize=9)
axb.text(-0.02, 1.02, "A", transform=axb.transAxes, fontsize=16, fontweight="bold", ha="right")

# ================= RIGHT: 3x5 rank-rank scatter matrix =================
inner = outer[0, 1].subgridspec(3, 5, wspace=0.32, hspace=0.55)
files = sorted(glob.glob(f"{PD}/*.tsv"))
order = sorted(files, key=lambda p: stat.get(os.path.basename(p)[:-4], (0,0,0,0,1))[4])  # worst 16S first
for k, p in enumerate(order):
    sp = os.path.basename(p)[:-4]
    ax = fig.add_subplot(inner[k // 5, k % 5])
    pr = [l.rstrip("\n").split("\t") for l in open(p)][1:]
    if pr:
        wgs = rankpct([float(r[2]) for r in pr])
        b2  = rankpct([float(r[3]) for r in pr]); s16r = rankpct([float(r[4]) for r in pr])
        ax.scatter(wgs, b2, s=5, color="#1f77b4", alpha=0.30, edgecolors="none")
        ax.scatter(wgs, s16r, s=5, color="#d62728", alpha=0.30, edgecolors="none")
        ax.plot([0,1],[0,1], color="#999", lw=0.7, ls="--")
    n, r2, _, _, r16, _, _ = stat.get(sp, (0,0,0,0,0,0,0))
    g, s = gs(sp)
    ax.set_title(f"$\\it{{{g[0]}.\\ {s}}}$\nρ {r2:.2f} / {r16:.2f}", fontsize=7.5)
    ax.set_xlim(-0.04, 1.04); ax.set_ylim(-0.04, 1.04)
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1]); ax.tick_params(labelsize=6.5); ax.grid(alpha=0.18)
# shared axis labels + legend + panel letter
fig.text(0.66, 0.045, "whole-genome distance (rank percentile)", ha="center", fontsize=10)
fig.text(0.352, 0.5, "2bRAD / 16S distance (rank percentile)", va="center", rotation=90, fontsize=10)
axb2 = fig.add_axes([0.355, 0.93, 0.001, 0.001]); axb2.axis("off")
fig.text(0.362, 0.955, "B", fontsize=16, fontweight="bold")
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([0],[0], marker="o", color="w", markerfacecolor="#1f77b4", label="2bRAD", markersize=7),
                    Line2D([0],[0], marker="o", color="w", markerfacecolor="#d62728", label="16S rRNA", markersize=7)],
           loc="upper right", bbox_to_anchor=(0.995, 0.965), fontsize=9, ncol=2, frameon=True)
fig.suptitle("2bRAD captures strain-level signal that 16S cannot (15 species): summary (A) and per-species rank–rank scatter (B)",
             fontsize=13, y=0.995)
fig.subplots_adjust(left=0.075, right=0.985, top=0.9, bottom=0.09)
fig.savefig(f"{PAPER}/figures/mash_2brad_vs_16s.png", dpi=150)
fig.savefig(f"{PAPER}/figures/mash_2brad_vs_16s.pdf")
print(f"2bRAD median {med2:.3f} vs 16S {med16:.3f}; wrote combined Fig 2 (A bars + B 3x5 scatter)")
