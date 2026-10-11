#!/usr/bin/env python3
"""Plot directional saliva concordance from the committed summary table.

The shotgun input was an R1 prefix subsample. The summary table provides call counts and
medians, so panel B reports group medians rather than reconstructing an unavailable raw distribution.
"""
from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PAPER = Path(__file__).resolve().parents[1]

def main():
    with (PAPER / "results/saliva_concordance.tsv").open() as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if not r["sample"].startswith("#")]
    totals = None
    p = None
    with (PAPER / "results/saliva_concordance.tsv").open() as fh:
        for line in fh:
            if line.startswith("#wms_confirmed"):
                totals = line.rstrip("\n").split("\t")
            elif line.startswith("#mannwhitney_p"):
                p = float(line.rstrip("\n").split("\t")[1])
    if totals is None or p is None:
        raise SystemExit("archived totals or p missing")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    ax = axes[0]
    x = np.arange(len(rows))
    shared = [int(r["wms_confirmed_by_2brad"]) for r in rows]
    only = [int(r["2brad_only"]) for r in rows]
    ax.bar(x, shared, color="#1f77b4", label="Shared with shotgun")
    ax.bar(x, only, bottom=shared, color="#ff7f0e", label="Native-2bRAD candidate only")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r['sample']}\n(subject {r['sample'].split('-')[1]})" for r in rows], fontsize=9)
    ax.set_ylabel("candidate strain-cluster calls")
    ax.set_title("A  Paired-sample call sets")
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.25, axis="y")

    ax = axes[1]
    med_shared = [float(r["med_relab_shared"]) for r in rows]
    med_only = [float(r["med_relab_2bradonly"]) for r in rows]
    positions = np.arange(len(rows))
    ax.bar(positions - 0.19, med_shared, width=0.36, color="#1f77b4", label="Shared median")
    ax.bar(positions + 0.19, med_only, width=0.36, color="#ff7f0e", label="Candidate-only median")
    ax.set_yscale("log")
    ax.set_xticks(positions)
    ax.set_xticklabels([r["sample"] for r in rows])
    ax.set_ylabel("community relative abundance (median)")
    ax.set_title(f"B  Candidate-only calls are lower abundance\np = {p:.1e}")
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.25, axis="y", which="both")

    fig.suptitle("Directional paired saliva comparison in three usable subjects", fontsize=12,
                 fontweight="bold")
    fig.text(0.5, 0.005,
             "Shotgun was an R1 prefix subsample; candidate-only calls are not independent strain validation.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=[0, 0.03, 1, 0.94])
    out = PAPER / "figures/saliva_fig08"
    fig.savefig(out.with_suffix(".png"), dpi=180, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(f"wrote {out}.png/.pdf; totals={totals[1]}")

if __name__ == "__main__":
    main()
