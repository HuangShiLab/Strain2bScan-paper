#!/usr/bin/env python3
"""Fig 4A schematic: why reference incompleteness breaks Jaccard clustering and how --containment
fixes it. An incomplete genome's tags are a subset of a complete relative's; Jaccard drops below the
0.95 cut and splits them (shared tags demoted to non-discriminating, detection corrupted), whereas
max-containment stays ~1 and merges them (cluster marker set = union of members, discriminating tags
preserved)."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

PAPER = "/Users/macstudio/Downloads/Strain2bScan-paper"
BLUE, GREEN, GREY, RED = "#1f77b4", "#2ca02c", "#bbbbbb", "#d62728"

fig, ax = plt.subplots(figsize=(12, 4.6)); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, fc="#f4f4f4", ec="#555", lw=1.2, r=0.02):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.2,rounding_size={r*100}",
                                fc=fc, ec=ec, lw=lw, mutation_aspect=0.5))

def tags(x, y, n_shared, n_extra, extra_col=GREEN, s=42):
    for i in range(n_shared):
        ax.scatter(x + i * 2.3, y, s=s, marker="s", color=BLUE, edgecolors="k", linewidths=0.3, zorder=3)
    for i in range(n_extra):
        ax.scatter(x + (n_shared + i) * 2.3, y, s=s, marker="s", color=extra_col, edgecolors="k", linewidths=0.3, zorder=3)

def arrow(x1, y1, x2, y2, col="#333"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14, color=col, lw=1.6))

# ---------- headers ----------
ax.text(25, 97, "Default (Jaccard)", ha="center", fontsize=13, fontweight="bold", color="#333")
ax.text(75, 97, "--containment (max-containment)", ha="center", fontsize=13, fontweight="bold", color=GREEN)
ax.plot([50, 50], [2, 92], color="#ccc", lw=1, ls=":")

for x0 in (0, 50):
    # input genomes (same in both columns)
    box(x0 + 5, 78, 40, 12, fc="#eef4fb")
    ax.text(x0 + 7, 86.5, "complete genome", fontsize=8.5, color="#333")
    tags(x0 + 9, 81.5, 6, 4)
    box(x0 + 5, 62, 40, 12, fc="#fbeeee")
    ax.text(x0 + 7, 70.5, "incomplete genome (~70%)", fontsize=8.5, color="#333")
    tags(x0 + 9, 65.5, 6, 0)
    ax.text(x0 + 9 + 6 * 2.3 + 2, 65.5, "= subset of the complete one", fontsize=7.5, va="center", color="#777", style="italic")

# ---------- left: Jaccard splits ----------
ax.text(25, 56, "Jaccard = |A∩B| / |A∪B| = 6/10 = 0.60  <  0.95", ha="center", fontsize=9.5, color=RED, fontweight="bold")
arrow(15, 61, 12, 44); arrow(35, 61, 38, 44)
box(3, 30, 20, 13, fc="#fff"); ax.text(13, 40.5, "cluster 1", ha="center", fontsize=8); tags(6, 35, 6, 4)
box(27, 30, 20, 13, fc="#fff"); ax.text(37, 40.5, "cluster 2", ha="center", fontsize=8); tags(30, 35, 6, 0)
ax.text(25, 25, "→ split. The 6 shared tags now sit in TWO clusters →", ha="center", fontsize=8, color="#333")
ax.text(25, 21, "demoted to $\\it{SharedPartial}$ (not discriminating).", ha="center", fontsize=8, color=RED)
ax.text(25, 13, "Reads (from the complete strain) match both →\nfalse positives + missed calls (precision & recall ↓)",
        ha="center", fontsize=8, color=RED)

# ---------- right: containment merges ----------
ax.text(75, 56, "max-containment = |A∩B| / min(|A|,|B|) = 6/6 = 1.00  ≥  0.95", ha="center", fontsize=9.5, color=GREEN, fontweight="bold")
arrow(65, 61, 75, 46)
box(58, 30, 34, 14, fc="#f0faf0", ec=GREEN, lw=1.6)
ax.text(75, 41.5, "one cluster  (marker set = UNION of members)", ha="center", fontsize=8)
tags(66, 35.5, 6, 4)
ax.text(75, 25, "→ merge. The complete member supplies the missing tags;", ha="center", fontsize=8, color="#333")
ax.text(75, 21, "shared+extra tags stay $\\it{ClusterSpecific}$ (discriminating).", ha="center", fontsize=8, color=GREEN)
ax.text(75, 13, "Reads match the one correct cluster →\ncorrect detection (accuracy restored)", ha="center", fontsize=8, color=GREEN)

# legend for tag colors
ax.scatter(30, 6.5, s=42, marker="s", color=BLUE, edgecolors="k", linewidths=0.3); ax.text(31.5, 6.5, "shared tag", va="center", fontsize=7.5)
ax.scatter(45, 6.5, s=42, marker="s", color=GREEN, edgecolors="k", linewidths=0.3); ax.text(46.5, 6.5, "tag only in the complete genome", va="center", fontsize=7.5)

fig.suptitle("Reference incompleteness fragments Jaccard clusters (wrong markers); max-containment keeps them intact",
             fontsize=12, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f"{PAPER}/figures/containment_mechanism.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{PAPER}/figures/containment_mechanism.pdf", bbox_inches="tight")
print("wrote figures/containment_mechanism.png/pdf")
