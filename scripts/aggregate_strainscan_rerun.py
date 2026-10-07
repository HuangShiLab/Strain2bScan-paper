#!/usr/bin/env python3
"""Build reviewer-facing aggregate tables from the reproducible StrainScan rerun."""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
WORK = PAPER / "work/strainscan_headtohead_rerun"
OUT = PAPER / "results/strainscan_rerun"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def finite(value: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number


def mean(rows: list[dict[str, str]], field: str) -> float:
    values = [finite(row[field]) for row in rows]
    values = [value for value in values if math.isfinite(value)]
    return sum(values) / len(values) if values else math.nan


def median(rows: list[dict[str, str]], field: str) -> float:
    values = sorted(finite(row[field]) for row in rows)
    values = [value for value in values if math.isfinite(value)]
    if not values:
        return math.nan
    middle = len(values) // 2
    if len(values) % 2:
        return values[middle]
    return (values[middle - 1] + values[middle]) / 2


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    manifest = read_tsv(WORK / "run_manifest.tsv")
    native_rows = read_tsv(
        PAPER / "figure_raw_data/sim_headtohead/strain2bscan_single_native_persample.tsv"
    )
    native_multi_rows = read_tsv(
        PAPER / "figure_raw_data/sim_headtohead/strain2bscan_multi_native_persample.tsv"
    )
    scored = read_tsv(WORK / "scored_runs.tsv")

    native_by_sample = {row["sample"]: row for row in native_rows}
    rerun_completed = {
        row["sample"]: row
        for row in scored
        if row["kind"] == "single" and row["status"] == "completed"
    }
    rerun_multi = {
        row["sample"]: row for row in scored if row["kind"] == "multi" and row["status"] == "completed"
    }
    native_multi_by_sample = {
        f'{row["depth"]}_{row["sample"]}': row for row in native_multi_rows
    }

    paired = []
    attempted_single = sum(row["kind"] == "single" for row in manifest)
    selected_single = sum(
        row["kind"] == "single" and "_rep1_" in row["sample"] for row in manifest
    )
    for sample, rerun in sorted(rerun_completed.items()):
        native = native_by_sample.get(sample)
        if native is None:
            continue
        paired.append(
            {
                "sample": sample,
                "species": rerun["species"],
                "depth": rerun["depth"],
                "s2b_precision": native["precision"],
                "s2b_recall": native["recall"],
                "s2b_f1": native["f1"],
                "ss_precision": rerun["precision"],
                "ss_recall": rerun["recall"],
                "ss_f1": rerun["f1"],
                "recall_difference": finite(native["recall"]) - finite(rerun["recall"]),
                "f1_difference": finite(native["f1"]) - finite(rerun["f1"]),
            }
        )
    paired_fields = [
        "sample", "species", "depth", "s2b_precision", "s2b_recall", "s2b_f1",
        "ss_precision", "ss_recall", "ss_f1", "recall_difference", "f1_difference",
    ]
    write_tsv(OUT / "single_paired_persample.tsv", paired_fields, paired)

    grouped_fields = [
        "group", "n", "s2b_precision", "s2b_recall", "s2b_f1",
        "ss_precision", "ss_recall", "ss_f1",
        "recall_difference", "f1_difference",
    ]

    def grouped(group: str, rows: list[dict[str, str]]) -> dict[str, object]:
        return {
            "group": group,
            "n": len(rows),
            "s2b_precision": median(rows, "s2b_precision"),
            "s2b_recall": median(rows, "s2b_recall"),
            "s2b_f1": median(rows, "s2b_f1"),
            "ss_precision": median(rows, "ss_precision"),
            "ss_recall": median(rows, "ss_recall"),
            "ss_f1": median(rows, "ss_f1"),
            "recall_difference": median(rows, "recall_difference"),
            "f1_difference": median(rows, "f1_difference"),
        }

    by_species = defaultdict(list)
    by_depth = defaultdict(list)
    for row in paired:
        by_species[row["species"]].append(row)
        by_depth[row["depth"]].append(row)
    aggregate_rows = [grouped(species, rows) for species, rows in sorted(by_species.items())]
    write_tsv(OUT / "single_by_species.tsv", grouped_fields, aggregate_rows)
    depth_rows = [grouped(depth, by_depth[depth]) for depth in ("0.5", "1", "3", "5", "10")]
    write_tsv(OUT / "single_by_depth.tsv", grouped_fields, depth_rows)

    multi_fields = [
        "sample", "n_truth_clusters", "n_strainscan_reports", "strainscan_artifact_dirs",
        "s2b_precision", "s2b_recall", "s2b_f1", "ss_precision", "ss_recall", "ss_f1",
    ]
    multi_rows = []
    for manifest_row in manifest:
        if manifest_row["kind"] != "multi":
            continue
        sample = manifest_row["sample"]
        rerun = rerun_multi.get(sample)
        native = native_multi_by_sample.get(sample)
        if rerun is None or native is None:
            continue
        species_dir = WORK / "results/multi" / sample
        species_dirs = [path for path in species_dir.iterdir() if path.is_dir()] if species_dir.exists() else []
        reports = [path for path in species_dirs if (path / "final_report.txt").is_file()]
        truth_path = Path(manifest_row["truth"])
        if not truth_path.is_absolute():
            truth_path = PAPER / truth_path
        n_truth = len(read_tsv(truth_path))
        multi_rows.append(
            {
                "sample": sample,
                "n_truth_clusters": n_truth,
                "n_strainscan_reports": len(reports),
                "strainscan_artifact_dirs": len(species_dirs),
                "s2b_precision": native["precision"],
                "s2b_recall": native["recall"],
                "s2b_f1": native["f1"],
                "ss_precision": rerun["precision"],
                "ss_recall": rerun["recall"],
                "ss_f1": rerun["f1"],
            }
        )
    write_tsv(OUT / "multi_persample.tsv", multi_fields, multi_rows)

    multi_by_depth = defaultdict(list)
    for row in multi_rows:
        multi_by_depth[row["sample"].split("_")[1]].append(row)
    multi_fields2 = [
        "group", "n", "n_truth_clusters", "n_strainscan_reports", "strainscan_artifact_dirs",
        "s2b_precision", "s2b_recall", "s2b_f1", "ss_precision", "ss_recall", "ss_f1",
    ]

    def multi_group(group: str, rows: list[dict[str, object]]) -> dict[str, object]:
        def avg(field: str) -> float:
            values = [float(row[field]) for row in rows]  # type: ignore[arg-type]
            return sum(values) / len(values)

        return {
            "group": group,
            "n": len(rows),
            "n_truth_clusters": avg("n_truth_clusters"),
            "n_strainscan_reports": avg("n_strainscan_reports"),
            "strainscan_artifact_dirs": avg("strainscan_artifact_dirs"),
            "s2b_precision": avg("s2b_precision"),
            "s2b_recall": avg("s2b_recall"),
            "s2b_f1": avg("s2b_f1"),
            "ss_precision": avg("ss_precision"),
            "ss_recall": avg("ss_recall"),
            "ss_f1": avg("ss_f1"),
        }

    multi_aggregate = [
        multi_group(depth, multi_by_depth[depth]) for depth in ("low", "med", "high")
    ]
    write_tsv(OUT / "multi_by_depth.tsv", multi_fields2, multi_aggregate)

    summary = {
        "manifest_single_runs": attempted_single,
        "published_matched_single_runs": selected_single,
        "completed_single_runs": len(rerun_completed),
        "paired_native_runs": len(paired),
        "no_call_published_matched_single_runs": selected_single - len(rerun_completed),
        "unattempted_full_manifest_single_runs": attempted_single - selected_single,
        "multi_samples": len(multi_rows),
        "strainscan_multi_reports": sum(int(row["n_strainscan_reports"]) for row in multi_rows),  # type: ignore[arg-type]
        "strainscan_multi_expected_reports": 12 * 15,
    }
    write_tsv(
        OUT / "rerun_summary.tsv",
        ["metric", "value"],
        [{"metric": key, "value": value} for key, value in summary.items()],
    )
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
