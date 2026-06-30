#!/usr/bin/env python3
"""Plot the reference-quality degradation figure.

Reads results/refqual_degradation.tsv (Strain2bScan arm) and, if present, the StrainScan arm
(refqual_ss/refqual_degradation_strainscan.tsv or results/refqual_degradation_strainscan.tsv),
overlaying both. Usage: python scripts/plot_refqual.py  (run from the repo root or a workdir).
"""
import csv, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    comp = [float(r["completeness"]) * 100 for r in rows]
    return (comp,
            [float(r["precision"]) for r in rows],
            [float(r["recall"]) for r in rows],
            [float(r["bray_curtis"]) for r in rows],
            [float(r["contamination"]) * 100 for r in rows])


s2b = next((p for p in ["results/refqual_degradation.tsv", "refqual/refqual_degradation.tsv"]
            if os.path.exists(p)), None)
if not s2b:
    sys.exit("refqual_degradation.tsv not found")
comp, prec, rec, bc, contam = load(s2b)

ss_path = next((p for p in ["results/refqual_degradation_strainscan.tsv",
                            "refqual_ss/refqual_degradation_strainscan.tsv"]
                if os.path.exists(p)), None)

fig, ax = plt.subplots(figsize=(6.4, 4.3))
ax.plot(comp, prec, "o-", color="#1f77b4", label="Strain2bScan precision")
ax.plot(comp, rec, "s-", color="#2ca02c", label="Strain2bScan recall")
ax.plot(comp, bc, "^--", color="#d62728", label="Strain2bScan Bray-Curtis")
for x, y, c in zip(comp, bc, contam):
    ax.annotate(f"{c:.0f}%", (x, y), textcoords="offset points", xytext=(0, 6),
                fontsize=7, color="#d62728")
if ss_path:
    c2, p2, r2, b2, _ = load(ss_path)
    ax.plot(c2, p2, "o:", color="#1f77b4", alpha=0.5, label="StrainScan precision")
    ax.plot(c2, r2, "s:", color="#2ca02c", alpha=0.5, label="StrainScan recall")
    ax.plot(c2, b2, "^:", color="#d62728", alpha=0.5, label="StrainScan Bray-Curtis")

ax.set_xlabel("Reference completeness (%)  [contamination 0->10%, contigs 1->400 co-vary]")
ax.set_ylabel("metric")
ax.set_title("Reference-genome quality vs strain-profiling accuracy\n"
             "(C. acnes, 14 truth strains degraded in the DB, 5 fixed samples)")
ax.invert_xaxis()
ax.set_ylim(0, 1)
ax.grid(alpha=0.3)
ax.legend(loc="center left", fontsize=8)
fig.tight_layout()
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/refqual_figure.png", dpi=150)
fig.savefig("figures/refqual_figure.pdf")
print("wrote figures/refqual_figure.png + .pdf")
