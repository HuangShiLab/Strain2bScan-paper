#!/usr/bin/env python3
"""P. copri panel-size test, BEFORE vs AFTER the strand-invariance digestion fix.

The original purpose was to test whether growing the reference panel fixes P. copri's apparent
over-detection. Diagnosing *why* it didn't uncovered the real cause: tag extraction was not
reverse-complement invariant (genomes digested forward-only, reads from both strands), so a
genome's reverse-strand tags collided with other genomes' forward tags. After fixing digestion
to scan both strands, P. copri profiles at precision 1.0 across all panel sizes -- the earlier
"recombination-driven false uniqueness" was a strand artifact. results/panelsize_prevotella.tsv."""
import csv, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import numpy as np

rows = list(csv.DictReader((l for l in open("results/panelsize_prevotella.tsv") if not l.startswith("#")), delimiter="\t"))
sizes = sorted({int(r["panel_size"]) for r in rows})
def series(cond, metric):
    d = {int(r["panel_size"]): float(r[metric]) for r in rows if r["condition"] == cond}
    return [d[s] for s in sizes]

fig, ax = plt.subplots(figsize=(7, 4.6))
x = np.arange(len(sizes))
ax.plot(x, series("before", "precision"), "o--", color="#c0392b", alpha=0.6, lw=1.8, label="precision (strand bug)")
ax.plot(x, series("before", "recall"), "s--", color="#2980b9", alpha=0.6, lw=1.8, label="recall (strand bug)")
ax.plot(x, series("after", "precision"), "o-", color="#c0392b", lw=2.5, label="precision (fixed)")
ax.plot(x, series("after", "recall"), "s-", color="#2980b9", lw=2.5, label="recall (fixed)")
ax.set_xticks(x); ax.set_xticklabels([f"{s}" for s in sizes])
ax.set_xlabel("reference panel size (genomes)"); ax.set_ylabel("accuracy"); ax.set_ylim(0, 1.08)
ax.legend(loc="center right", fontsize=8)
ax.set_title("P. copri: a strand-invariance digestion fix, not a bigger panel\n(precision 0.19 → 1.0; the 'recombination false uniqueness' was an artifact)")
ax.annotate("both-strand digestion:\nnear-identical genomes now cluster correctly\n(40→24, 80→44, 112→53 clusters) and\nreverse-strand tags no longer cross-detect",
            (0.03, 0.52), xycoords="axes fraction", ha="left", va="center", fontsize=7.5,
            color="#333333", bbox=dict(boxstyle="round", fc="#eef7ee", ec="#bcd"))
fig.tight_layout(); os.makedirs("figures", exist_ok=True)
fig.savefig("figures/panelsize_prevotella.png", dpi=150)
fig.savefig("figures/panelsize_prevotella.pdf")
print("wrote figures/panelsize_prevotella.png + .pdf")
