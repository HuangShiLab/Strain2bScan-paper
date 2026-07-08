#!/usr/bin/env python3
"""P. copri panel-size test: does growing the reference panel fix the recombination-driven
false-uniqueness precision problem? Prediction was that it would. It did NOT -- precision (and
recall) DECLINE as the panel grows 40 -> 80 -> 111 (StrainScan's whole P. copri panel), on the
identical 5 samples. Strictly-nested panels, Strain2bScan only. results/panelsize_prevotella.tsv."""
import csv, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

rows = list(csv.DictReader(open("results/panelsize_prevotella.tsv"), delimiter="\t"))
n = [int(r["panel_size"]) for r in rows]
P = [float(r["precision"]) for r in rows]
R = [float(r["recall"]) for r in rows]

fig, ax = plt.subplots(figsize=(6.6, 4.4))
ax.plot(n, P, "o-", color="#d62728", lw=2, label="precision")
ax.plot(n, R, "s--", color="#1f77b4", lw=2, label="recall")
for x, y in zip(n, P):
    ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points", xytext=(0, -14), ha="center", fontsize=8, color="#d62728")
for x, y in zip(n, R):
    ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8, color="#1f77b4")
ax.set_xlabel("reference panel size (genomes; all singleton clusters)")
ax.set_ylabel("accuracy")
ax.set_ylim(0, 1.0)
ax.set_xticks(n)
ax.legend(loc="upper right")
ax.set_title("P. copri: growing the panel does NOT fix false uniqueness\n(prediction refuted — precision & recall both decline)")
ax.annotate("more genomes -> more singleton clusters whose\n'unique' markers a mosaic strain spuriously matches\n(FP up), while the true strain's own uniques get\nreclassified as shared (TP down)",
            (0.5, 0.40), xycoords="axes fraction", ha="center", va="center", fontsize=7.5,
            color="#555555", bbox=dict(boxstyle="round", fc="#f5f5f5", ec="#cccccc"))
fig.tight_layout(); os.makedirs("figures", exist_ok=True)
fig.savefig("figures/panelsize_prevotella.png", dpi=150)
fig.savefig("figures/panelsize_prevotella.pdf")
print("wrote figures/panelsize_prevotella.png + .pdf")
