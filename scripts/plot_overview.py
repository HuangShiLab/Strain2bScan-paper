#!/usr/bin/env python3
"""Fig 1 - overview schematic: Strain2bScan pipeline, two data modes, Fast2bRAD-M interoperability.

Hand-drawn schematic (no data inputs). Kept as a script so the figure is reproducible/editable
and consistent with the repo's matplotlib PDF+PNG convention.
"""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# --- palette (2 ramps): blue = pipeline/data ; green = 2bRAD data modes + Fast2bRAD-M interop ---
BLUE_F, BLUE_E, DB_F = "#e7f0f9", "#1f77b4", "#cfe0f2"
GREEN_F, GREEN_E = "#e4f1e8", "#2e8b57"
TXT, ARR = "#12354f", "#5a6b7a"

def box(ax, cx, cy, w, h, text, fc=BLUE_F, ec=BLUE_E, tc=TXT, fs=8.5, weight="normal"):
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                 boxstyle="round,pad=0.15,rounding_size=1.4",
                 fc=fc, ec=ec, lw=1.4, zorder=3))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs, color=tc, weight=weight, zorder=5)

def arrow(ax, x1, y1, x2, y2, label=None, color=ARR, rad=0.0, ls="-", lw=1.7):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls,
                                connectionstyle=f"arc3,rad={rad}", shrinkA=2, shrinkB=2), zorder=2)
    if label:
        ax.text((x1 + x2) / 2, (y1 + y2) / 2, label, ha="center", va="center", fontsize=7,
                color=color, style="italic", zorder=6,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.92))

fig, ax = plt.subplots(figsize=(11, 6.7))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

# title + column headers + legend
ax.text(50, 97.5, "Strain2bScan: strain-level profiling on 2bRAD-reduced markers",
        ha="center", va="center", fontsize=13.5, weight="bold", color=TXT)
ax.text(50, 92.7,
        "blue = reference build / pipeline      green = 2bRAD data modes  &  Fast2bRAD-M interoperability",
        ha="center", va="center", fontsize=8, color="#6b7a88")
ax.text(25, 88.5, "A   Database construction  (per species)", ha="center", fontsize=10.5,
        weight="bold", color=BLUE_E)
ax.text(74, 88.5, "B   Sample profiling", ha="center", fontsize=10.5, weight="bold", color=BLUE_E)

# ---------------- A: database construction (left column, cx=25) ----------------
box(ax, 25, 82, 36, 8,  "Reference genomes\n$G_1\\ G_2\\ \\dots\\ G_n$", weight="bold")
box(ax, 25, 64, 36, 9,  "Single-copy 2bRAD tags\n32–38 bp  ·  ~1–2% of genome  ·  50–100× sparser than full $k$-mers")
box(ax, 25, 47.5, 36, 7, "Within-species clusters", fc=DB_F)
box(ax, 25, 31, 36, 8,  "Marker classes\nspecies-core / cluster-specific / strain-specific")
box(ax, 25, 13, 36, 8,  "cluster × marker\ndatabase", fc=DB_F, weight="bold")

arrow(ax, 25, 78, 25, 68.7, "2bRAD digest · 16 type-IIB enzymes (BcgI + …)")
arrow(ax, 25, 59.4, 25, 51.2, "cluster · 0.95 Jaccard · MinHash")
arrow(ax, 25, 44, 25, 35.2, "classify markers")
arrow(ax, 25, 27, 25, 17.2)

# ---------------- B: sample profiling (right column, cx=74) ----------------
box(ax, 65, 82, 17.5, 8, "Shotgun\nmetagenome reads", fc=GREEN_F, ec=GREEN_E, fs=8)
box(ax, 84, 82, 15.5, 8, "BcgI 2bRAD\nreads (= tags)", fc=GREEN_F, ec=GREEN_E, fs=8)
box(ax, 65, 71.5, 17.5, 5, "in-silico 2bRAD digest", fs=8)
box(ax, 74, 61, 35, 6, "canonical marker counts", weight="bold")
box(ax, 74, 46, 35, 7, "species gate  (Fast2bRAD-M markers)", fc=GREEN_F, ec=GREEN_E)
box(ax, 74, 30, 35, 8, "strain detection + abundance\nunique markers · depth / regression")
box(ax, 74, 13, 35, 8, "strain-level profile\nclusters present + relative abundance", fc=DB_F, weight="bold")

arrow(ax, 65, 78, 65, 74.2)                    # shotgun -> digest
arrow(ax, 65, 69, 70.5, 64.2)                  # digest -> counts
arrow(ax, 84, 78, 79, 64.2)                    # 2bRAD -> counts (direct)
arrow(ax, 74, 58, 74, 49.7, "Layer-1")
arrow(ax, 74, 42.5, 74, 34.2, "Layer-2")
arrow(ax, 74, 26, 74, 17.2)

# cross-link: reference DB feeds the matching step
arrow(ax, 43.2, 14, 56.2, 29, "reference DB", color=BLUE_E, rad=-0.18, ls=(0, (4, 2)), lw=1.6)

# ---------------- interoperability band ----------------
box(ax, 49.25, 4, 88, 5.4,
    "Tags identical to Fast2bRAD-M / 2bRADExtraction.pl  →  "
    "species layer (Fast2bRAD-M) ∘ strain layer (Strain2bScan) share one tag space",
    fc=GREEN_F, ec=GREEN_E, tc="#1d5637", fs=8.2, weight="bold")

fig.tight_layout()
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/overview.png", dpi=150, bbox_inches="tight")
fig.savefig("figures/overview.pdf", bbox_inches="tight")
print("wrote figures/overview.png + .pdf")
