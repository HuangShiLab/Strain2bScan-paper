#!/usr/bin/env python3
"""Run the current Strain2bScan build on the complete 2,025-sample single-species simulation pool.

The script rebuilds Strain2bScan species databases from the committed simulation pool, reruns
every matched single-species sample, scores cluster sets in Strain2bScan's own cluster space, and
compares precision/recall/F1 with the archived native benchmark. It intentionally does not modify
the frozen primary benchmark; outputs are written to a separate work directory and summary table.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
DEFAULT_SOFTWARE = PAPER.parent / "Strain2bScan"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def all_single_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if row["kind"] == "single"]


def metrics(truth: set[str], pred: set[str]) -> tuple[float, float, float, int, int, int]:
    tp, fp, fn = len(truth & pred), len(pred - truth), len(truth - pred)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, recall, f1, tp, fp, fn


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_db(binary: Path, species: str, source: Path, outdir: Path, threads: int) -> Path:
    db = outdir / species / "db.tsv"
    db.parent.mkdir(parents=True, exist_ok=True)
    if db.exists() and db.stat().st_size:
        return db
    env = {**os.environ, "STRAIN2BSCAN_THREADS": str(threads)}
    subprocess.run(
        [str(binary), "cluster", "--genomes", str(source), "--enzyme", "all",
         "--similarity", "0.95", "--out", str(db)],
        check=True, env=env,
    )
    if not db.exists() or not db.stat().st_size:
        raise RuntimeError(f"database build failed for {species}")
    return db


def profile_sample(binary: Path, db: Path, row: dict[str, str], outdir: Path,
                   threads: int, run_index: int, fixed_gate: bool = False) -> set[str]:
    out = outdir / "predictions" / row["sample"] / "pred.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    fifo = outdir / "fifo" / f"run_{run_index:04d}.fq"
    fifo.parent.mkdir(parents=True, exist_ok=True)
    if fifo.exists():
        fifo.unlink()
    os.mkfifo(fifo)
    env = {**os.environ, "STRAIN2BSCAN_THREADS": str(threads)}
    # Start the reader first so opening the FIFO for writing does not deadlock.
    proc = subprocess.Popen(
        [str(binary), "profile", "--db", str(db), "--reads", str(fifo), "--out", str(out)]
             + (["--fixed-gate"] if fixed_gate else []),
        env=env,
    )
    with fifo.open("wb") as fifo_handle:
        writer = subprocess.Popen(["gzip", "-dc", row["r1"], row["r2"]], stdout=fifo_handle,
                              stderr=subprocess.DEVNULL)
    try:
        returncode = proc.wait()
        writer.wait()
        if returncode:
            raise subprocess.CalledProcessError(returncode, proc.args)
    finally:
        if writer.poll() is None:
            writer.terminate()
            writer.wait()
        if fifo.exists():
            fifo.unlink()
    calls = set()
    for line in out.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        calls.add(line.split("\t", 1)[0])
    return calls


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path,
                    default=DEFAULT_SOFTWARE / "target/release/strain2bscan")
    ap.add_argument("--source-repo", type=Path, default=DEFAULT_SOFTWARE)
    ap.add_argument("--work", type=Path, default=PAPER / "work/sim2025_current")
    ap.add_argument("--summary", type=Path,
                    default=PAPER / "results/sim2025_current/software_regression.tsv")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0, help="smoke-test only the first N runs")
    ap.add_argument("--fixed-gate", action="store_true",
                    help="use the archived benchmark compatibility mode (fixed singleton/floor gate)")
    args = ap.parse_args()

    binary = args.binary.resolve()
    if not binary.is_file():
        raise SystemExit(f"binary not found: {binary}; run cargo build --release")
    source_repo = args.source_repo.resolve()
    commit = subprocess.check_output(["git", "-C", str(source_repo), "rev-parse", "HEAD"], text=True).strip()
    manifest = read_tsv(PAPER / "work/sim2025_current/run_manifest_2025.tsv")
    rows = all_single_rows(manifest)
    if args.limit:
        rows = rows[:args.limit]
    if len(rows) != 2025 and not args.limit:
        raise SystemExit(f"expected 2,025 single-species rows, got {len(rows)}")

    pool_rows = read_tsv(PAPER / "figure_raw_data/sim_pool_manifest.tsv")
    genome_to_cluster = {(r["species"], r["accession"]): r["cluster"] for r in pool_rows}
    by_species = defaultdict(list)
    for row in rows:
        by_species[row["species"]].append(row)

    old = {r["sample"]: r for r in read_tsv(
        PAPER / "figure_raw_data/sim_headtohead/strain2bscan_single_native_persample.tsv")}

    output_rows = []
    n_done = 0
    for species, species_rows in sorted(by_species.items()):
        print(f"[db] {species}", flush=True)
        db = build_db(binary, species, PAPER / f"figure_raw_data/sim_genome_pool/{species}",
                      args.work / "db", args.threads)
        for row in species_rows:
            truth = {
                genome_to_cluster[(species, fields[0])]
                for fields in (line.rstrip("\n").split("\t") for line in Path(row["truth"]).read_text().splitlines())
                if fields and fields[0] and not fields[0].startswith("#")
            }
            pred = profile_sample(binary, db, row, args.work, args.threads, n_done + 1,
                                  fixed_gate=args.fixed_gate)
            cp, cr, cf, tp, fp, fn = metrics(truth, pred)
            old_row = old[row["sample"]]
            op, orr, of = (float(old_row[k]) for k in ("precision", "recall", "f1"))
            exact = all(round(a, 4) == round(b, 4) for a, b in ((cp, op), (cr, orr), (cf, of)))
            output_rows.append({
                "sample": row["sample"], "species": species,
                "truth_n": len(truth), "pred_n": len(pred), "tp": tp, "fp": fp, "fn": fn,
                "current_precision": cp, "current_recall": cr, "current_f1": cf,
                "archived_precision": op, "archived_recall": orr, "archived_f1": of,
                "precision_difference": cp - op, "recall_difference": cr - orr,
                "f1_difference": cf - of, "metrics_exact": exact,
            })
            n_done += 1
            print(f"[profile {n_done:3d}/{len(rows)}] {row['sample']} "
                  f"P={cp:.3f} R={cr:.3f} F1={cf:.3f}", flush=True)

    # A compact source table in the native-benchmark schema makes downstream aggregation direct.
    source_fields = ["species","sample","strategy","k","rep","depth","n_truth_cl","n_pred_cl",
                     "TP","FP","FN","precision","recall","f1"]
    source_path = args.summary.parent / "current_single_persample.tsv"
    with source_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=source_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in output_rows:
            parts = row["sample"].split("__", 1)[1].split("_")
            strategy, k, rep, depth = parts[0], parts[1].removeprefix("k"), parts[2].removeprefix("rep"), parts[3].removeprefix("d")
            writer.writerow({
                "species": row["species"], "sample": row["sample"], "strategy": strategy,
                "k": k, "rep": rep, "depth": depth,
                "n_truth_cl": row["truth_n"], "n_pred_cl": row["pred_n"],
                "TP": row["tp"], "FP": row["fp"], "FN": row["fn"],
                "precision": row["current_precision"], "recall": row["current_recall"],
                "f1": row["current_f1"],
            })

    fields = list(output_rows[0])
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output_rows)

    mismatches = [r for r in output_rows if not r["metrics_exact"]]
    provenance = {
        "configuration": "current default depth-adaptive gates",
        "source_commit": commit,
        "binary": str(binary),
        "binary_sha256": sha256(binary),
        "runs": len(output_rows),
        "archived_metric_differences": len(mismatches),
        "max_abs_precision_difference": max((abs(r["precision_difference"]) for r in output_rows), default=0.0),
        "max_abs_recall_difference": max((abs(r["recall_difference"]) for r in output_rows), default=0.0),
        "max_abs_f1_difference": max((abs(r["f1_difference"]) for r in output_rows), default=0.0),
    }
    (args.summary.parent / "current_software_regression_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
