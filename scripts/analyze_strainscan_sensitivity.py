#!/usr/bin/env python3
"""Failure-aware sensitivity and database-resolution audit for the StrainScan rerun.

The primary rerun reports medians over the 204 completed StrainScan pairs. This script retains
all 225 matched published-run entries. For the 21 StrainScan no-call runs it assigns precision,
recall and F1 of zero; the corresponding Strain2bScan run is unchanged. It also reports the
number of reference genomes collapsed to clusters by each tool's own 0.95 database.

Outputs are written under results/strainscan_rerun/:
  single_all225_persample.tsv
  single_all225_by_species.tsv
  single_all225_by_depth.tsv
  sensitivity_failure_as_absence.tsv
  database_cluster_resolution.tsv
"""
from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

PAPER = Path(__file__).resolve().parents[1]
OUT = PAPER / "results/strainscan_rerun"
RNG = np.random.default_rng(20261006)


def read_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def median(values):
    values = sorted(values)
    if not values:
        return math.nan
    mid = len(values) // 2
    return values[mid] if len(values) % 2 else (values[mid - 1] + values[mid]) / 2


def metric_group(rows):
    return {
        "n": len(rows),
        "no_call": sum(r["status"] == "no_call" for r in rows),
        "s2b_precision": median(float(r["s2b_precision"]) for r in rows),
        "s2b_recall": median(float(r["s2b_recall"]) for r in rows),
        "s2b_f1": median(float(r["s2b_f1"]) for r in rows),
        "ss_precision": median(float(r["ss_precision"]) for r in rows),
        "ss_recall": median(float(r["ss_recall"]) for r in rows),
        "ss_f1": median(float(r["ss_f1"]) for r in rows),
        "recall_difference": median(float(r["recall_difference"]) for r in rows),
        "f1_difference": median(float(r["f1_difference"]) for r in rows),
    }


def main():
    native = {
        r["sample"]: r
        for r in read_tsv(PAPER / "figure_raw_data/sim_headtohead/strain2bscan_single_native_persample.tsv")
    }
    scored = read_tsv(PAPER / "work/strainscan_headtohead_rerun/scored_runs.tsv")
    selected = [r for r in scored if r["kind"] == "single" and "_rep1_" in r["sample"]]
    if len(selected) != 225:
        raise SystemExit(f"expected 225 selected single runs, got {len(selected)}")

    rows = []
    for scored_row in sorted(selected, key=lambda r: r["sample"]):
        sample = scored_row["sample"]
        native_row = native.get(sample)
        if native_row is None:
            raise SystemExit(f"missing native Strain2bScan metrics for {sample}")
        no_call = scored_row["status"] == "no_call"
        ss_precision = 0.0 if no_call else fnum(scored_row["precision"])
        ss_recall = 0.0 if no_call else fnum(scored_row["recall"])
        ss_f1 = 0.0 if no_call else fnum(scored_row["f1"])
        s2b_precision = fnum(native_row["precision"])
        s2b_recall = fnum(native_row["recall"])
        s2b_f1 = fnum(native_row["f1"])
        rows.append({
            "sample": sample,
            "species": scored_row["species"],
            "depth": scored_row["depth"],
            "status": scored_row["status"],
            "n_truth_s2b_clusters": int(native_row["n_truth_cl"]),
            "n_truth_strainscan_clusters": int(scored_row["n_truth"]),
            "s2b_precision": s2b_precision,
            "s2b_recall": s2b_recall,
            "s2b_f1": s2b_f1,
            "ss_precision": ss_precision,
            "ss_recall": ss_recall,
            "ss_f1": ss_f1,
            "recall_difference": s2b_recall - ss_recall,
            "f1_difference": s2b_f1 - ss_f1,
        })

    fields = list(rows[0])
    write_tsv(OUT / "single_all225_persample.tsv", fields, rows)

    by_species = defaultdict(list)
    by_depth = defaultdict(list)
    for row in rows:
        by_species[row["species"]].append(row)
        by_depth[row["depth"]].append(row)
    group_fields = ["group", "n", "no_call", "s2b_precision", "s2b_recall", "s2b_f1",
                    "ss_precision", "ss_recall", "ss_f1", "recall_difference", "f1_difference"]
    species_rows = [dict(group=species, **metric_group(by_species[species])) for species in sorted(by_species)]
    write_tsv(OUT / "single_all225_by_species.tsv", group_fields, species_rows)
    depth_rows = [dict(group=depth, **metric_group(by_depth[depth])) for depth in ("0.5", "1", "3", "5", "10")]
    write_tsv(OUT / "single_all225_by_depth.tsv", group_fields, depth_rows)

    by_species_pairs = defaultdict(list)
    for row in rows:
        by_species_pairs[row["species"]].append(row)
    species_names = sorted(by_species_pairs)
    point_recall = float(np.mean([r["recall_difference"] for r in rows]))
    point_f1 = float(np.mean([r["f1_difference"] for r in rows]))
    stats = []
    for _ in range(20000):
        chosen = RNG.choice(species_names, len(species_names), replace=True)
        sampled = [r for species in chosen for r in by_species_pairs[species]]
        stats.append((np.mean([r["recall_difference"] for r in sampled]),
                      np.mean([r["f1_difference"] for r in sampled])))
    ci_recall = np.quantile([x[0] for x in stats], [0.025, 0.975])
    ci_f1 = np.quantile([x[1] for x in stats], [0.025, 0.975])
    summary = [{
        "analysis": "all_225_failure_as_absence",
        "n_entries": len(rows),
        "no_call_entries": sum(r["status"] == "no_call" for r in rows),
        "recall_difference_point": point_recall,
        "recall_difference_ci_low": ci_recall[0],
        "recall_difference_ci_high": ci_recall[1],
        "f1_difference_point": point_f1,
        "f1_difference_ci_low": ci_f1[0],
        "f1_difference_ci_high": ci_f1[1],
        "bootstrap": 20000,
        "resampling": "species clusters with replacement",
    }]
    write_tsv(OUT / "sensitivity_failure_as_absence.tsv",
              list(summary[0]), summary)

    pool = defaultdict(list)
    for row in read_tsv(PAPER / "figure_raw_data/sim_pool_manifest.tsv"):
        pool[row["species"]].append(row)
    db = PAPER / "work/strainscan_headtohead_rerun/db"
    resolution = []
    for species in sorted(pool):
        s2b_clusters = len({row["cluster"] for row in pool[species]})
        map_path = db / species / "Cluster_Result/hclsMap_95_recls.txt"
        strainscan_clusters = 0
        if map_path.exists():
            with map_path.open() as fh:
                for line in fh:
                    if line.strip():
                        strainscan_clusters += 1
        resolution.append({
            "species": species,
            "pool_genomes": len(pool[species]),
            "s2b_0.95_clusters": s2b_clusters,
            "strainscan_0.95_clusters": strainscan_clusters,
        })
    write_tsv(OUT / "database_cluster_resolution.tsv",
              ["species", "pool_genomes", "s2b_0.95_clusters", "strainscan_0.95_clusters"], resolution)
    print(f"wrote failure-aware all-225 sensitivity ({len(rows)} entries; "
          f"{sum(r['status'] == 'no_call' for r in rows)} scored as no-call)")
    print(f"failure-aware mean differences: recall={point_recall:.4f}, F1={point_f1:.4f}")


if __name__ == "__main__":
    main()
