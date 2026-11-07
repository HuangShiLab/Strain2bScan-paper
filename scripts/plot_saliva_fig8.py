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
    class_abundance = read_tsv("results/wms32_rebuild/call_class_abundance_summary.tsv")
    if len(pairs) != 32:
        raise SystemExit("expected the complete 32-pair full-depth rebuild")

    def sort_key(row):
        prefix, subject = row["sample"].split("-", 1)
        return int(subject), int(prefix.removeprefix("S"))
    pairs = sorted(pairs, key=sort_key)
    labels = [row["sample"] for row in pairs]
    shared = np.array([int(row["shared_calls"]) for row in pairs])
    wms_only = np.array([int(row["wms_only_calls"]) for row in pairs])
    native_only = np.array([int(row["native_only_calls"]) for row in pairs])
    class_lookup = {(r["modality"], r["group"]): r for r in class_abundance}
    abundance_specs = [
        ("WMS", "WMS only", float(class_lookup[("WMS", "shared")]["fraction_of_modality_detected_abundance"]),
         float(class_lookup[("WMS", "WMS_only")]["fraction_of_modality_detected_abundance"]),
         int(class_lookup[("WMS", "shared")]["n_calls"]), int(class_lookup[("WMS", "WMS_only")]["n_calls"])),
        ("Native 2bRAD", "Native only",
         float(class_lookup[("native_2bRAD", "shared")]["fraction_of_modality_detected_abundance"]),
         float(class_lookup[("native_2bRAD", "native_only")]["fraction_of_modality_detected_abundance"]),
         int(class_lookup[("native_2bRAD", "shared")]["n_calls"]),
         int(class_lookup[("native_2bRAD", "native_only")]["n_calls"])),
    ]

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
    y = np.arange(len(abundance_specs))
    shared_pct = np.array([row[2] * 100 for row in abundance_specs])
    unique_pct = np.array([row[3] * 100 for row in abundance_specs])
    shared_calls = np.array([row[4] for row in abundance_specs])
    unique_calls = np.array([row[5] for row in abundance_specs])
    unique_names = [row[1] for row in abundance_specs]
    ax.barh(y, shared_pct, height=0.52, color="#2166AC", label="Shared")
    ax.barh(y, unique_pct, left=shared_pct, height=0.52, color="#F4A582",
            label="Modality-only calls")
    for i in range(len(abundance_specs)):
        ax.text(shared_pct[i] / 2, y[i], f"{shared_pct[i]:.1f}%\n{shared_calls[i]:,} calls",
                ha="center", va="center", color="white", fontsize=7.6)
        ax.text(shared_pct[i] + unique_pct[i] / 2, y[i],
                f"{unique_pct[i]:.1f}%\n{unique_calls[i]:,} calls",
                ha="center", va="center", color="black", fontsize=6.7)
    ax.set_yticks(y)
    ax.set_yticklabels([row[0] for row in abundance_specs])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0", "25", "50", "75", "100%"])
    ax.set_xlabel("Percent of modality-detected sample fraction")
    ax.set_title("B  Unique calls are abundance-poor",
                 loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.6, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
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
