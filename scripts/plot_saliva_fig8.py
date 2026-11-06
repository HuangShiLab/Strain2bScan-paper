#!/usr/bin/env python3
"""Build the full-depth paired saliva comparison figure from the 32-pair rebuild."""
from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PAPER = Path(__file__).resolve().parents[1]


def read_tsv(path):
    with (PAPER / path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    pairs = read_tsv("results/wms32_rebuild/per_sample_summary.tsv")
    medians = read_tsv("results/wms32_rebuild/per_sample_shared_native_median.tsv")
    abundance = {r["group"]: r for r in
                 read_tsv("results/wms32_rebuild/shared_native_abundance_summary.tsv")}
    if len(pairs) != 32 or len(medians) != 32:
        raise SystemExit("expected the complete 32-pair full-depth rebuild")

    def sort_key(row):
        prefix, subject = row["sample"].split("-", 1)
        return int(subject), int(prefix.removeprefix("S"))
    pairs = sorted(pairs, key=sort_key)
    medians = sorted(medians, key=sort_key)
    labels = [row["sample"] for row in pairs]
    shared = np.array([int(row["shared_calls"]) for row in pairs])
    wms_only = np.array([int(row["wms_only_calls"]) for row in pairs])
    native_only = np.array([int(row["native_only_calls"]) for row in pairs])
    shared_median = np.array([float(row["shared_median_sample_fraction"]) for row in medians])
    native_median = np.array([float(row["native_only_median_sample_fraction"]) for row in medians])
    p_value = float(abundance["paired_per_sample_median_wilcoxon_p"]["median_sample_fraction"])

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.7), gridspec_kw={"width_ratios": [2.15, 1]})
    ax = axes[0]
    x = np.arange(len(pairs))
    ax.bar(x, shared, width=0.72, color="#2166AC", label="Shared")
    ax.bar(x, wms_only, bottom=shared, width=0.72, color="#BDD7E7", label="WMS only")
    ax.bar(x, native_only, bottom=shared + wms_only, width=0.72,
           color="#F4A582", label="Native 2bRAD only")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, fontsize=6.0)
    ax.set_ylabel("strain-cluster calls")
    ax.set_ylim(0, max(shared + wms_only + native_only) * 1.16)
    ax.set_title("A  Gate-matched full-depth paired calls (n = 32)",
                 loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper center")
    ax.grid(axis="y", color="0.90", linewidth=0.6)
    ax.set_axisbelow(True)

    ax = axes[1]
    for i in range(len(medians)):
        ax.plot([0, 1], [shared_median[i], native_median[i]], color="0.70",
                linewidth=0.55, alpha=0.65, zorder=1)
    rng = np.random.default_rng(20261106)
    ax.scatter(rng.normal(0, 0.045, len(shared_median)), shared_median, s=17,
               color="#2166AC", alpha=0.55, linewidths=0, zorder=2)
    ax.scatter(rng.normal(1, 0.045, len(native_median)), native_median, s=17,
               color="#D7301F", alpha=0.55, linewidths=0, zorder=2)
    ax.boxplot([shared_median, native_median], positions=[0, 1], widths=0.42,
               patch_artist=True, showfliers=False, medianprops={"color": "black", "linewidth": 1.0},
               boxprops={"facecolor": "none", "edgecolor": "0.45", "linewidth": 0.8},
               whiskerprops={"color": "0.45", "linewidth": 0.8},
               capprops={"color": "0.45", "linewidth": 0.8}, zorder=3)
    ax.set_yscale("log")
    ax.set_xlim(-0.5, 1.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Shared", "Native 2bRAD only"])
    ax.set_ylabel("per-sample median sample fraction")
    ax.set_title("B  Native-only calls are lower abundance",
                 loc="left", fontweight="bold")
    ax.text(0.97, 0.94, f"paired Wilcoxon p = {p_value:.1g}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="0.75", linewidth=0.5))
    fig.suptitle("Full-depth paired saliva 2bRAD-M versus shotgun WMS: 8 subjects × 4 timepoints",
                 x=0.012, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0.035, 1, 0.95])
    fig.text(0.012, 0.012,
             "Same 19-species oral panel and gates; differences in depth, representation and host content prevent an accuracy interpretation.",
             ha="left", fontsize=8.5, color="0.25")
    out = PAPER / "figures/source_panels/saliva_fig08"
    fig.savefig(out.with_suffix(".png"), dpi=180, bbox_inches="tight", facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    print(f"wrote {out}.png/.pdf; n={len(pairs)}")

if __name__ == "__main__":
    main()
