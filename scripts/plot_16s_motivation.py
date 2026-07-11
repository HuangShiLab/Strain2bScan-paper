#!/usr/bin/env python3
"""Motivation figure: 2bRAD captures strain-level signal that 16S cannot. Per-species Spearman
correlation of between-strain distances with whole-genome distance — 2bRAD (high) vs 16S (low),
across the 15 slide-deck species."""
import os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

PAPER = os.path.expanduser("~/Downloads/Strain2bScan-paper")
rows = [l.rstrip("\n").split("\t") for l in open(f"{PAPER}/results/mash_2brad_vs_16s_recomputed.tsv")
        if l.strip() and not l.startswith("#") and not l.startswith("species")]
d = [(r[0].replace("_", " "), float(r[3]), float(r[4])) for r in rows]
d.sort(key=lambda x: x[1] - x[2])          # by 2bRAD advantage
names = [x[0] for x in d]; b2 = [x[1] for x in d]; b16 = [x[2] for x in d]
y = np.arange(len(names)); h = 0.38

fig, ax = plt.subplots(figsize=(9, 7.5))
ax.barh(y + h/2, b2, h, color="#1f77b4", label="2bRAD", edgecolor="k", linewidth=0.3)
ax.barh(y - h/2, b16, h, color="#d62728", label="16S rRNA", edgecolor="k", linewidth=0.3)
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=9, style="italic")
ax.set_xlabel("Spearman correlation of between-strain distance with whole-genome distance")
ax.set_xlim(0, 1.02); ax.set_ylim(-0.7, len(names) - 0.3)
med2, med16 = np.median(b2), np.median(b16)
ax.axvline(med2, color="#1f77b4", ls="--", lw=1.2, alpha=0.8)
ax.axvline(med16, color="#d62728", ls="--", lw=1.2, alpha=0.8)
ax.text(med2 - 0.01, -0.62, f"2bRAD median {med2:.2f}", color="#1f77b4", fontsize=8.5, ha="right", va="bottom")
ax.text(med16 + 0.01, -0.62, f"16S median {med16:.2f}", color="#d62728", fontsize=8.5, ha="left", va="bottom")
ax.set_title("2bRAD captures strain-level signal that 16S cannot\n(between-strain distance vs whole-genome, 15 species)", fontsize=12)
ax.legend(loc="lower right", fontsize=10); ax.grid(alpha=0.25, axis="x")
fig.tight_layout()
fig.savefig(f"{PAPER}/figures/mash_2brad_vs_16s.png", dpi=150)
fig.savefig(f"{PAPER}/figures/mash_2brad_vs_16s.pdf")
print(f"wrote figures/mash_2brad_vs_16s.png + .pdf")
print(f"2bRAD median Spearman {med2:.3f}  vs  16S median {med16:.3f}")
