#!/usr/bin/env python3
"""Aggregate the current Strain2bScan regression against the reproducible StrainScan rerun."""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

PAPER = Path(__file__).resolve().parents[1]
OUT = PAPER / "results/strainscan_rerun"
RNG = np.random.default_rng(20261006)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def median(rows: list[dict[str, str]], field: str) -> float:
    values = sorted(float(row[field]) for row in rows)
    if not values:
        return math.nan
    mid = len(values) // 2
    return values[mid] if len(values) % 2 else (values[mid - 1] + values[mid]) / 2


def metric_group(rows: list[dict[str, str]], prefix: str) -> dict[str, float]:
    return {f"{prefix}_{metric}": median(rows, f"{prefix}_{metric}") for metric in ("precision", "recall", "f1")}


def main() -> int:
    current = read_tsv(OUT / "current_single_persample.tsv")
    scored = read_tsv(PAPER / "work/strainscan_headtohead_rerun/scored_runs.tsv")
    strainscan = {
        row["sample"]: row for row in scored
        if row["kind"] == "single" and row["status"] == "completed"
    }
    paired = []
    for row in current:
        comparator = strainscan.get(row["sample"])
        if comparator is None:
            continue
        paired.append({
            "sample": row["sample"], "species": row["species"], "depth": row["depth"],
            "s2b_precision": row["precision"], "s2b_recall": row["recall"], "s2b_f1": row["f1"],
            "ss_precision": comparator["precision"], "ss_recall": comparator["recall"],
            "ss_f1": comparator["f1"],
            "precision_difference": float(row["precision"]) - float(comparator["precision"]),
            "recall_difference": float(row["recall"]) - float(comparator["recall"]),
            "f1_difference": float(row["f1"]) - float(comparator["f1"]),
        })
    paired_fields = ["sample", "species", "depth", "s2b_precision", "s2b_recall", "s2b_f1",
                     "ss_precision", "ss_recall", "ss_f1", "precision_difference",
                     "recall_difference", "f1_difference"]
    write_tsv(OUT / "current_paired_persample.tsv", paired_fields, paired)

    by_species = defaultdict(list)
    by_depth = defaultdict(list)
    for row in paired:
        by_species[row["species"]].append(row)
        by_depth[row["depth"]].append(row)
    aggregate_fields = ["group", "n", "s2b_precision", "s2b_recall", "s2b_f1",
                        "ss_precision", "ss_recall", "ss_f1",
                        "precision_difference", "recall_difference", "f1_difference"]
    species_rows = []
    for species in sorted(by_species):
        rows = by_species[species]
        species_rows.append({"group": species, "n": len(rows), **metric_group(rows, "s2b"),
                             **metric_group(rows, "ss"),
                             "precision_difference": median(rows, "precision_difference"),
                             "recall_difference": median(rows, "recall_difference"),
                             "f1_difference": median(rows, "f1_difference")})
    write_tsv(OUT / "current_single_by_species.tsv", aggregate_fields, species_rows)
    depth_rows = []
    for depth in ("0.5", "1", "3", "5", "10"):
        rows = by_depth[depth]
        depth_rows.append({"group": depth, "n": len(rows), **metric_group(rows, "s2b"),
                           **metric_group(rows, "ss"),
                           "precision_difference": median(rows, "precision_difference"),
                           "recall_difference": median(rows, "recall_difference"),
                           "f1_difference": median(rows, "f1_difference")})
    write_tsv(OUT / "current_single_by_depth.tsv", aggregate_fields, depth_rows)

    by_species_pairs = defaultdict(list)
    for row in paired:
        by_species_pairs[row["species"]].append(row)
    species = sorted(by_species_pairs)
    stats = []
    for _ in range(20000):
        chosen = RNG.choice(species, len(species), replace=True)
        selected = [row for name in chosen for row in by_species_pairs[name]]
        stats.append((np.mean([row["recall_difference"] for row in selected]),
                      np.mean([row["f1_difference"] for row in selected])))
    point_recall = float(np.mean([row["recall_difference"] for row in paired]))
    point_f1 = float(np.mean([row["f1_difference"] for row in paired]))
    ci = np.quantile(stats, [0.025, 0.975], axis=0)
    uncertainty = [{
        "analysis": "current_completed_single_species:recall_difference",
        "point": point_recall, "ci_low": ci[0, 0], "ci_high": ci[1, 0],
    }, {
        "analysis": "current_completed_single_species:f1_difference",
        "point": point_f1, "ci_low": ci[0, 1], "ci_high": ci[1, 1],
    }]
    # Failure-as-absence sensitivity: recall is defined because absent predictions have zero recall;
    # precision and F1 are not defined for execution/no-report failures and are not reported here.
    completed_samples = set(strainscan)
    failure_rows = []
    for row in current:
        comparator = strainscan.get(row["sample"])
        ss_r = float(comparator["recall"]) if comparator else 0.0
        failure_rows.append({
            "sample": row["sample"], "species": row["species"], "status": "completed" if comparator else "no_call",
            "s2b_recall": row["recall"], "ss_recall": ss_r,
            "recall_difference": float(row["recall"]) - ss_r,
        })
    failure_recall = float(np.mean([row["recall_difference"] for row in failure_rows]))
    failure_species = sorted({row["species"] for row in failure_rows})
    failure_by_species = defaultdict(list)
    for row in failure_rows:
        failure_by_species[row["species"]].append(row)
    failure_stats = []
    for _ in range(20000):
        chosen = RNG.choice(failure_species, len(failure_species), replace=True)
        selected = [row for name in chosen for row in failure_by_species[name]]
        failure_stats.append(np.mean([row["recall_difference"] for row in selected]))
    failure_ci = np.quantile(failure_stats, [0.025, 0.975])
    uncertainty.extend([
        {"analysis": "current_failure_as_absence:recall_difference", "point": failure_recall,
         "ci_low": failure_ci[0], "ci_high": failure_ci[1]},
        {"analysis": "current_failure_as_absence:no_call_rate", "point": 21/225,
         "ci_low": math.nan, "ci_high": math.nan},
    ])
    write_tsv(OUT / "current_uncertainty_summary.tsv",
              ["analysis", "point", "ci_low", "ci_high"], uncertainty)
    write_tsv(OUT / "current_failure_as_absence.tsv",
              ["sample", "species", "status", "s2b_recall", "ss_recall", "recall_difference"], failure_rows)
    print(f"paired={len(paired)} recall={point_recall:.6f} f1={point_f1:.6f} "
          f"failure_recall={failure_recall:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
