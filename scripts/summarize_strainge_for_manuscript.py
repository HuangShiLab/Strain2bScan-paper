#!/usr/bin/env python3
"""Create the source-derived StrainGST summary used by manuscript Table 9."""
from __future__ import annotations

import csv
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/strainge_rerun"


def read(path: str) -> list[dict[str, str]]:
    with (ROOT / path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def median(rows: list[dict[str, str]], field: str) -> float:
    values = [float(row[field]) for row in rows if math.isfinite(float(row[field]))]
    return statistics.median(values) if values else math.nan


def main() -> int:
    accuracy = read("results/strainge_rerun/simulation_persample.tsv")
    efficiency = {
        row["sample"]: row
        for row in read("results/strainge_rerun/simulation_run_efficiency.tsv")
    }
    rows = []
    for kind, label in (("single", "sim_single"), ("multi", "sim_multi")):
        selected = [row for row in accuracy if row["kind"] == kind]
        wall = median([efficiency[row["sample"]] for row in selected], "total_wall_s")
        rss = median([efficiency[row["sample"]] for row in selected], "maxrss_bytes")
        rows.append({
            "dataset": label,
            "n": len(selected),
            "precision": median(selected, "precision"),
            "recall": median(selected, "recall"),
            "f1": median(selected, "f1"),
            "median_wall_s": wall,
            "median_peak_rss_gb": rss / 1e9,
        })
    accuracy = read("results/strainge_rerun/mock_persample.tsv")
    efficiency = {
        row["sample"]: row
        for row in read("results/strainge_rerun/mock_run_efficiency.tsv")
    }
    for row in accuracy:
        eff = efficiency[row["sample"]]
        rows.append({
            "dataset": row["sample"],
            "n": 1,
            "precision": float(row["precision"]),
            "recall": float(row["recall"]),
            "f1": float(row["f1"]),
            "median_wall_s": float(eff["total_wall_s"]),
            "median_peak_rss_gb": int(eff["maxrss_bytes"]) / 1e9,
        })
    fields = ["dataset", "n", "precision", "recall", "f1", "median_wall_s", "median_peak_rss_gb"]
    with (OUT / "manuscript_summary.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    # Runtime context uses the same 225 matched simulations and 12 communities
    # as the StrainGST rerun. Accuracy is not pooled because cluster spaces differ.
    s2b_single = [
        row for row in read("figure_raw_data/sim_headtohead/strain2bscan_single_native_persample.tsv")
        if "_rep1_" in row["sample"]
        and (row["strategy"] == "diff" or row["species"] == "Mycobacterium_tuberculosis")
    ]
    s2b_multi = [
        row
        for row in read("figure_raw_data/sim_headtohead/strain2bscan_multi_native_persample.tsv")
        if row["sample"] in {f"sample0{index}" for index in range(1, 5)}
    ]
    comparisons = []
    for label, selected, strainge in (
        ("sim_single", s2b_single, rows[0]),
        ("sim_multi", s2b_multi, rows[1]),
    ):
        s2b_wall = statistics.median(float(row["wall_s"]) for row in selected)
        s2b_rss = statistics.median(float(row["maxrss_mb"]) for row in selected) / 1024
        comparisons.append({
            "dataset": label,
            "n": len(selected),
            "s2b_median_wall_s": s2b_wall,
            "strainge_median_wall_s": strainge["median_wall_s"],
            "strainge_s2b_wall_ratio": strainge["median_wall_s"] / s2b_wall,
            "s2b_median_peak_rss_gb": s2b_rss,
            "strainge_median_peak_rss_gb": strainge["median_peak_rss_gb"],
            "strainge_s2b_rss_ratio": strainge["median_peak_rss_gb"] / s2b_rss,
        })
    with (OUT / "runtime_comparison_summary.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(comparisons[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(comparisons)
    print(f"wrote {OUT / 'runtime_comparison_summary.tsv'}")
    print(f"wrote {OUT / 'manuscript_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
