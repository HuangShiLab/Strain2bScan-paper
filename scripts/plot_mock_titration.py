#!/usr/bin/env python3
"""Fig 9 - real ATCC MSA-1002 DNA mock community (native BcgI 2bRAD): Strain2bScan holds
precision 1.0 and full 20/20 species recall down to 0.1 ng input, on real error-containing reads."""
import os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, numpy as np


def rows(path):
    return [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]


r = rows("results/mock_msa1002_titration.tsv")
h = r[0]
d = [dict(zip(h, x)) for x in r[1:]]
inp = [float(x["dna_input_ng"]) for x in d]
recall = [float(x["species_recall"]) for x in d]
prec = [float(x["species_precision"]) for x in d]
reads = [int(x["reads"]) for x in d]
ndet = [int(x["n_detected"]) for x in d]

fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.plot(inp, prec, "-s", color="#2ca02c", lw=2, ms=8, label="species precision")
ax.plot(inp, recall, "-o", color="#1f77b4", lw=2, ms=8, label="species recall (detected / 20)")
ax.set_xscale("log")
ax.set_xlabel("input DNA (ng)  —  low-biomass axis"); ax.set_ylabel("score")
ax.set_ylim(0, 1.06); ax.set_xticks(inp); ax.set_xticklabels(["0.001", "0.01", "0.1", "1"])
ax.set_title("Real ATCC MSA-1002 mock (native BcgI 2bRAD):\nprecision 1.0 throughout, full recall down to 0.1 ng", fontsize=11)
ax.grid(alpha=0.25, which="both")
# annotate reads + detected count at each point
for x, y, n, rd in zip(inp, recall, ndet, reads):
    ax.annotate(f"{n}/20\n({rd/1000:.0f}k reads)", (x, y), textcoords="offset points",
                xytext=(0, -28 if y > 0.5 else 10), ha="center", fontsize=7.5, color="#1f77b4")
ax.axhspan(0.99, 1.06, color="#2ca02c", alpha=0.06)
ax.legend(loc="lower right", fontsize=9)
fig.tight_layout(); os.makedirs("figures", exist_ok=True)
fig.savefig("figures/mock_msa1002_titration.png", dpi=150)
fig.savefig("figures/mock_msa1002_titration.pdf")
print("wrote figures/mock_msa1002_titration.png + .pdf")
