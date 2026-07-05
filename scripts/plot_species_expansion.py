#!/usr/bin/env python3
"""Species-expansion figure: Strain2bScan vs real StrainScan (StrainScan's own pre-built DBs,
profiling-only -- no dashing needed) on 3 additional StrainScan-native species (Akkermansia
muciniphila, Prevotella copri, M. tuberculosis), each on a 40-genome subset drawn FROM
StrainScan's own DB accession list (guarantees genome overlap between the two tools' DBs).
StrainScan did not complete on M. tuberculosis within a practical time/memory budget (DNF)."""
import csv, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import numpy as np

rows = list(csv.DictReader(open("results/species_expansion.tsv"), delimiter="\t"))
labels = [r["species"].replace("_", " ").replace("Mycobacterium", "M.").replace("Akkermansia", "A.").replace("Prevotella", "P.") for r in rows]
x = np.arange(len(rows)); w = 0.35

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

# Panel A: precision/recall
ax = axes[0]
s2b_p = [float(r["s2b_precision"]) for r in rows]; s2b_r = [float(r["s2b_recall"]) for r in rows]
ss_p = [float(r["ss_precision"]) if r["ss_precision"] != "DNF" else np.nan for r in rows]
ss_r = [float(r["ss_recall"]) if r["ss_recall"] != "DNF" else np.nan for r in rows]
ax.bar(x - 1.5*w/2, s2b_p, w/2, label="Strain2bScan P", color="#1f77b4")
ax.bar(x - 0.5*w/2, s2b_r, w/2, label="Strain2bScan R", color="#7fb3d9")
ax.bar(x + 0.5*w/2, ss_p, w/2, label="StrainScan P", color="#8c8c8c")
ax.bar(x + 1.5*w/2, ss_r, w/2, label="StrainScan R", color="#c9c9c9")
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=8)
ax.set_ylabel("accuracy"); ax.set_ylim(0, 1.15); ax.legend(fontsize=7, ncol=2)
ax.set_title("Precision / recall per species")
for i, r in enumerate(rows):
    if r["ss_precision"] == "DNF":
        ax.annotate("DNF", (i + 0.5*w/2, 0.05), ha="center", fontsize=9, color="#d62728", fontweight="bold")

# Panel B: time (log)
ax = axes[1]
s2b_t = [float(r["s2b_profile_s"]) for r in rows]
ss_t = [float(r["ss_profile_s"]) if "DNF" not in r["ss_profile_s"] else None for r in rows]
ax.bar(x - w/2, s2b_t, w, label="Strain2bScan", color="#1f77b4")
ax.bar([xi + w/2 for xi, t in zip(x, ss_t) if t is not None], [t for t in ss_t if t is not None], w, label="StrainScan", color="#8c8c8c")
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=8)
ax.set_ylabel("profile time / sample (s, log)"); ax.legend(fontsize=8)
ax.set_title("Profiling time per sample")
ymax_t = max(max(s2b_t), max(t for t in ss_t if t is not None))
ax.set_ylim(top=ymax_t * 30)
for i, r in enumerate(rows):
    if "DNF" in r["ss_profile_s"]:
        ax.annotate("StrainScan\nDNF\n(>3.3 h,\n1 sample)", (i, ymax_t * 5), ha="center", va="bottom",
                    fontsize=7.5, color="#d62728", fontweight="bold", clip_on=False, zorder=10)

# Panel C: memory (log)
ax = axes[2]
s2b_m = [float(r["s2b_rss_mb"]) for r in rows]
ss_m = [float(r["ss_rss_mb"]) if "DNF" not in r["ss_rss_mb"] else None for r in rows]
ax.bar(x - w/2, s2b_m, w, label="Strain2bScan", color="#1f77b4")
ax.bar([xi + w/2 for xi, m in zip(x, ss_m) if m is not None], [m for m in ss_m if m is not None], w, label="StrainScan", color="#8c8c8c")
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=8)
ax.set_ylabel("peak RSS (MB, log)"); ax.legend(fontsize=8)
ax.set_title("Peak memory")
ymax_m = max(max(s2b_m), max(m for m in ss_m if m is not None))
ax.set_ylim(top=ymax_m * 30)
for i, r in enumerate(rows):
    if "DNF" in r["ss_rss_mb"]:
        ax.annotate("StrainScan\nDNF\n(>25.5 GB,\nstill rising)", (i, ymax_m * 5), ha="center", va="bottom",
                    fontsize=7.5, color="#d62728", fontweight="bold", clip_on=False, zorder=10)

fig.suptitle("Species expansion: Strain2bScan vs StrainScan's own pre-built DBs (3 new species)", y=1.03)
fig.tight_layout(); os.makedirs("figures", exist_ok=True)
fig.savefig("figures/species_expansion.png", dpi=150, bbox_inches="tight")
fig.savefig("figures/species_expansion.pdf", bbox_inches="tight")
print("wrote figures/species_expansion.png + .pdf")
