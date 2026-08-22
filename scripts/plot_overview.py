#!/usr/bin/env python3
"""Fig 1 - Strain2bScan algorithm overview.

Hand-drawn schematic (no data inputs), kept as a script so the figure stays reproducible and
consistent with the repo's matplotlib PDF+PNG convention. Panels: (A) database construction,
(B) per-sample profiling with the two input modes, (C) what actually decides a call.

Every number shown is the shipped default, read from Strain2bScan src/identify.rs and src/main.rs.
"""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

# palette: blue = pipeline/reference ; green = 2bRAD data modes + interoperability ; amber = decisions
BLUE_F, BLUE_E, DB_F = "#e7f0f9", "#1f77b4", "#cfe0f2"
GREEN_F, GREEN_E = "#e4f1e8", "#2e8b57"
AMB_F, AMB_E = "#fdf0dc", "#c8791a"
TXT, ARR, MUT = "#12354f", "#5a6b7a", "#6b7a88"


def box(ax, cx, cy, w, h, text, fc=BLUE_F, ec=BLUE_E, tc=TXT, fs=7.6, weight="normal"):
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                 boxstyle="round,pad=0.15,rounding_size=1.2",
                 fc=fc, ec=ec, lw=1.3, zorder=3))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs, color=tc,
            weight=weight, zorder=5, linespacing=1.45)


def arrow(ax, x1, y1, x2, y2, label=None, color=ARR, rad=0.0, ls="-", lw=1.5, fs=6.5):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls,
                                connectionstyle=f"arc3,rad={rad}", shrinkA=2, shrinkB=2), zorder=2)
    if label:
        ax.text((x1 + x2) / 2, (y1 + y2) / 2, label, ha="center", va="center", fontsize=fs,
                color=color, style="italic", zorder=6,
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none", alpha=0.94))


fig = plt.figure(figsize=(13.2, 7.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

ax.text(50, 97.6, "Strain2bScan: strain-level profiling on a 1–2 % genomic subsample",
        ha="center", va="center", fontsize=13, weight="bold", color=TXT)

# ============================== A: database construction ==============================
AX = 16
ax.text(AX, 92.6, "A   Database construction  (per species)", ha="center", fontsize=9.6,
        weight="bold", color=BLUE_E)

box(ax, AX, 86, 27, 6.4, "Reference genomes  $G_1 \\dots G_n$", weight="bold")
box(ax, AX, 72, 27, 9.5,
    "Single-copy 2bRAD tags\n25–33 bp  ·  ~1–2 % of the genome\ncanonical tag → 64-bit marker")
box(ax, AX, 57, 27, 7.0, "Within-species clusters\n(finest unit short reads resolve)", fc=DB_F)
box(ax, AX, 42.5, 27, 8.5,
    "Marker classes\nspecies-core / shared-partial /\ncluster- and strain-specific")
box(ax, AX, 27, 27, 8.0, "cluster × marker database\n+ inverted degree index", fc=DB_F, weight="bold")

arrow(ax, AX, 82.6, AX, 77.0, "digest · 16 type-IIB enzymes")
arrow(ax, AX, 67.0, AX, 60.7, "single-linkage · Jaccard 0.95\nMinHash above 96 genomes")
arrow(ax, AX, 53.3, AX, 47.0, "within-species incidence")
arrow(ax, AX, 38.1, AX, 31.2, "unique = absent from every\nother cluster at any copy number")

# ============================== B: sample profiling ==============================
BX = 50
ax.text(BX, 92.6, "B   Sample profiling", ha="center", fontsize=9.6, weight="bold", color=BLUE_E)

box(ax, BX - 8.5, 86, 15, 6.4, "Shotgun / WMS reads", fc=GREEN_F, ec=GREEN_E)
box(ax, BX + 8.5, 86, 15, 6.4, "Native BcgI 2bRAD\nreads (already tags)", fc=GREEN_F, ec=GREEN_E)
box(ax, BX - 8.5, 77.5, 15, 4.6, "in-silico digest", fs=7.2)
box(ax, BX, 68.5, 30, 5.6, "canonical marker counts  $c_m$", weight="bold")
ax.text(BX, 63.9, "the sample is digested ONCE, then matched against every species database",
        ha="center", fontsize=6.6, color=MUT, style="italic")

box(ax, BX, 57.5, 30, 7.4,
    "Layer-1 · species gate\npresent ≥ max(⌈200·r⌉, 10),  r = max(1−e$^{-λ}$, 0.25)",
    fc=AMB_F, ec=AMB_E)
box(ax, BX, 45.5, 30, 6.6,
    "cross-species marker restriction\nkeep only tags specific to this species panel-wide",
    fc=AMB_F, ec=AMB_E)
box(ax, BX, 33, 30, 8.2,
    "Layer-2 · detection\nsupport ≥ 8  ·  coverage ≥ 0.1\ncoverage / (1−e$^{-depth}$) ≥ 0.5",
    fc=AMB_F, ec=AMB_E)
box(ax, BX, 20, 30, 7.4,
    "depth$_j$ = mean over the WHOLE panel\n(zeros included, top 1 % winsorized)")
box(ax, BX, 8, 30, 7.6,
    "abundance  ·  global_abundance  ·  sample_fraction\nwithin-species  ·  cells  ·  DNA",
    fc=DB_F, weight="bold")

arrow(ax, BX - 8.5, 82.6, BX - 8.5, 80.0)
arrow(ax, BX - 8.5, 75.1, BX - 4, 71.5)
arrow(ax, BX + 8.5, 82.6, BX + 4, 71.5)
arrow(ax, BX, 65.6, BX, 61.4)
arrow(ax, BX, 53.7, BX, 49.0)
arrow(ax, BX, 42.1, BX, 37.3)
arrow(ax, BX, 28.8, BX, 23.9)
arrow(ax, BX, 16.2, BX, 12.0)

arrow(ax, 30.2, 27, 34.4, 33, "reference DB", color=BLUE_E, rad=-0.16, ls=(0, (4, 2)))

# ============================== C: what decides a call ==============================
CX = 84
ax.text(CX, 92.6, "C   What separates a real strain from an artefact", ha="center",
        fontsize=9.6, weight="bold", color=AMB_E)

# C1 - depth-breadth consistency, drawn as the actual curve
axc = fig.add_axes([0.715, 0.545, 0.225, 0.245])
lam = np.linspace(0.02, 9, 400)
axc.plot(lam, 1 - np.exp(-lam), color=BLUE_E, lw=2)
axc.text(4.4, 0.90, "expected breadth  $1-e^{-\\lambda}$", fontsize=6.4, color=BLUE_E, ha="center")

# Both points are measured, from the shadow experiment: near-identical breadth, 17x apart in depth.
axc.scatter([0.44], [0.392], s=54, color=GREEN_E, zorder=5, edgecolor="white", linewidth=1.1)
axc.annotate("genuine rare strain\n0.44x, breadth 0.39", (0.44, 0.392), (1.5, 0.17),
             fontsize=6.2, color=GREEN_E, ha="left",
             arrowprops=dict(arrowstyle="-", color=GREEN_E, lw=0.9))
axc.scatter([7.68], [0.350], s=60, color="#c0392b", zorder=5, marker="X",
            edgecolor="white", linewidth=1.1)
axc.annotate("shadow\n7.7x, breadth 0.35", (7.68, 0.350), (5.0, 0.40),
             fontsize=6.2, color="#c0392b", ha="left",
             arrowprops=dict(arrowstyle="-", color="#c0392b", lw=0.9))
# the gap the test reads
axc.plot([7.68, 7.68], [0.350, 1 - np.exp(-7.68)], color="#c0392b", lw=0.9, ls=(0, (2, 2)))
axc.text(7.55, 0.68, "rejected", fontsize=6.0, color="#c0392b", ha="right", rotation=90, va="center")

axc.set_xlabel("estimated depth $\\lambda$ (reads / tag)", fontsize=6.6)
axc.set_ylabel("observed breadth", fontsize=6.6)
axc.set_ylim(0, 1.08); axc.set_xlim(0, 9.2)
axc.tick_params(labelsize=6)
axc.set_title("near-identical breadth, 17x apart in depth", fontsize=7, color=TXT, pad=3)
for sp in ("top", "right"):
    axc.spines[sp].set_visible(False)

box(ax, CX, 41, 30, 9.0,
    "A cluster is called only if its breadth matches\nthe depth it claims. A strain carrying part of\n"
    "another cluster's loci fires it at full depth\nacross a fraction of the panel — and is rejected.",
    fc=AMB_F, ec=AMB_E, fs=7.0)

box(ax, CX, 26.5, 30, 9.6,
    "Depth-adaptive evidence\n$c\\geq2$ above 3 reads/tag (errors dominate singletons)\n"
    "$c\\geq1$ below it — at λ=0.5, 78 % of detected\nmarkers are singletons and a fixed rule discards them",
    fc=AMB_F, ec=AMB_E, fs=7.0)

box(ax, CX, 12, 30, 8.4,
    "Zero-inclusive depth\nAveraging only DETECTED markers pins a rare\n"
    "cluster near 1 read/tag however rare it is,\nflattening the whole composition.",
    fc=AMB_F, ec=AMB_E, fs=7.0)

# ============================== interoperability band ==============================
box(ax, 50, 2.2, 96, 3.6,
    "Tags are identical to Fast2bRAD-M / 2bRADExtraction.pl  →  species layer (Fast2bRAD-M) ∘ strain layer (Strain2bScan) share one tag space   ·   "
    "Rust, no third-party dependencies   ·   streaming gzip I/O, peak RSS ~11 MB on 4 M reads",
    fc=GREEN_F, ec=GREEN_E, tc="#1d5637", fs=7.2, weight="bold")

os.makedirs("figures", exist_ok=True)
fig.savefig("figures/overview.png", dpi=200, bbox_inches="tight")
fig.savefig("figures/overview.pdf", bbox_inches="tight")
print("wrote figures/overview.png + .pdf")
