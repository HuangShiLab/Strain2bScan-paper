#!/usr/bin/env python3
"""Motivation figure v2: high-quality (complete/near-complete) genomes only, per-species genome
count in the labels, bootstrap 95% CI error bars, genus and species on separate lines."""
import os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

PAPER = os.path.expanduser("~/Downloads/Strain2bScan-paper")
rows = [l.rstrip("\n").split("\t") for l in open(f"{PAPER}/results/mash_2brad_vs_16s_recomputed.tsv")
        if l.strip() and not l.startswith("#") and not l.startswith("species")]
# species n_hq n_pairs sp_2brad lo hi sp_16s lo hi n_16s
d = []
for r in rows:
    sp = r[0]; n = int(r[1])
    s2, lo2, hi2 = float(r[3]), float(r[4]), float(r[5])
    s16, lo16, hi16 = float(r[6]), float(r[7]), float(r[8])
    d.append((sp, n, s2, lo2, hi2, s16, lo16, hi16))
d.sort(key=lambda x: x[2] - x[5])                 # by 2bRAD advantage

def two_line(sp):
    g, s = sp.split("_", 1)
    return f"{g}\n{s.replace('_', ' ')}"

labels = [two_line(x[0]) for x in d]; ns = [x[1] for x in d]
s2 = np.array([x[2] for x in d]); e2 = np.array([[x[2]-x[3] for x in d], [x[4]-x[2] for x in d]])
s16 = np.array([x[5] for x in d]); e16 = np.array([[x[5]-x[6] for x in d], [x[7]-x[5] for x in d]])
y = np.arange(len(labels)); h = 0.38

fig, ax = plt.subplots(figsize=(9.5, 8.6))
ax.barh(y + h/2, s2, h, color="#1f77b4", label="2bRAD", edgecolor="k", linewidth=0.3,
        xerr=e2, error_kw=dict(ecolor="#0b3d66", elinewidth=1.1, capsize=2.5))
ax.barh(y - h/2, s16, h, color="#d62728", label="16S rRNA", edgecolor="k", linewidth=0.3,
        xerr=e16, error_kw=dict(ecolor="#7a1416", elinewidth=1.1, capsize=2.5))
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8.5, style="italic")
ax.set_xlabel("Spearman correlation of between-strain distance with whole-genome distance")
ax.set_xlim(-0.25, 1.06); ax.set_ylim(-0.7, len(labels) - 0.3)
ax.axvline(0, color="#888", lw=0.8)
med2, med16 = np.median(s2), np.median(s16)
ax.axvline(med2, color="#1f77b4", ls="--", lw=1.2, alpha=0.8)
ax.axvline(med16, color="#d62728", ls="--", lw=1.2, alpha=0.8)
# per-row genome count at the right margin
for yi, n in zip(y, ns):
    ax.text(1.045, yi, f"n={n}", fontsize=7.5, va="center", ha="right", color="#444")
ax.text(1.045, len(labels)-0.5, "genomes", fontsize=7.5, va="center", ha="right", color="#444", style="italic")
ax.text(med2, len(labels) - 0.35, f"2bRAD median {med2:.2f}", color="#1f77b4", fontsize=8.5, ha="center", va="bottom")
ax.text(med16, -0.62, f"16S median {med16:.2f}", color="#d62728", fontsize=8.5, ha="center", va="bottom")
ax.set_title("2bRAD captures strain-level signal that 16S cannot\n"
             "(between-strain vs whole-genome distance; complete/near-complete genomes;\n"
             "error bars = 95% CI over genome subsamples)", fontsize=11.5)
ax.legend(loc="lower left", fontsize=10, framealpha=0.95); ax.grid(alpha=0.25, axis="x")
fig.tight_layout()
fig.savefig(f"{PAPER}/figures/mash_2brad_vs_16s.png", dpi=150)
fig.savefig(f"{PAPER}/figures/mash_2brad_vs_16s.pdf")
print(f"wrote figures/mash_2brad_vs_16s.png + .pdf")
print(f"2bRAD median {med2:.3f} vs 16S median {med16:.3f}; genomes/species {min(ns)}-{max(ns)}")
