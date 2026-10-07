#!/usr/bin/env python3
"""Deterministically rerun the StrainScan arm of the 15-species head-to-head.

The original scratchpad drivers were lost. This script reconstructs the run manifest from the committed
simulation reads/truth and the pinned genome-pool manifest. It uses the StrainScan v1.0.14 container named
``ss`` (image digest recorded in ``results/benchmark_configuration.json``).

Examples
--------
Build one database and run one smoke sample::

  python3 scripts/rerun_strainscan_headtohead.py --build --profile --species Phocaeicola_dorei --limit 1

Build all databases and run the complete 204 single-species plus 12 multi-species manifest::

  python3 scripts/rerun_strainscan_headtohead.py --build --profile
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
READS = PAPER / "figure_raw_data/sim_single_species"
MULTI = PAPER / "figure_raw_data/sim_multi_species"
POOL = PAPER / "figure_raw_data/sim_genome_pool"
WORK = PAPER / "work/strainscan_headtohead_rerun"
CONFIG = PAPER / "results/benchmark_configuration.json"
CONTAINER = "ss"
PYBIN = "/opt/conda/envs/ss/bin/python"
TOOLDIR = "/host/StrainScan"
PAPERDIR = "/host/Strain2bScan-paper"

def sh(args: list[str], *, log: Path | None = None) -> None:
    print("+", " ".join(args), flush=True)
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("w") as fh:
            subprocess.run(args, stdout=fh, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(args, check=True)

def docker_args(*args: str) -> list[str]:
    return ["docker", "exec", "-w", TOOLDIR, "-e", "PATH=/opt/conda/envs/ss/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin", CONTAINER, *args]

def ensure_container() -> None:
    status = subprocess.run(["docker", "inspect", "-f", "{{.State.Running}}", CONTAINER],
                            capture_output=True, text=True)
    if status.returncode or status.stdout.strip() != "true":
        sh(["docker", "start", CONTAINER])

def species_list() -> list[str]:
    with (PAPER / "figure_raw_data/sim_pool_summary.tsv").open() as fh:
        return [r["species"] for r in csv.DictReader(fh, delimiter="\t")]

def manifest() -> list[dict[str, str]]:
    rows = []
    species = species_list()
    for sp in species:
        for strategy in (["same"] if sp == "Mycobacterium_tuberculosis" else ["diff"]):
            for k in (2, 3, 5):
                for rep in range(1, 6):
                    for depth in ("0.5", "1", "3", "5", "10"):
                        stem = f"{sp}__{strategy}_k{k}_rep{rep}_d{depth}"
                        r1 = READS / sp / "reads" / f"{stem}_R1.fastq.gz"
                        r2 = READS / sp / "reads" / f"{stem}_R2.fastq.gz"
                        truth = READS / sp / "truth" / f"{stem}.truth.tsv"
                        if not all(p.exists() for p in (r1, r2, truth)):
                            continue
                        rows.append({"kind": "single", "sample": stem, "species": sp,
                                     "r1": r1.relative_to(PAPER).as_posix(),
                                     "r2": r2.relative_to(PAPER).as_posix(),
                                     "truth": truth.relative_to(PAPER).as_posix()})
    for depth in ("low", "med", "high"):
        for n in range(1, 5):
            sample = f"depth_{depth}_sample{n:02d}"
            r1 = MULTI / f"depth_{depth}/reads/sample{n:02d}_R1.fastq.gz"
            r2 = MULTI / f"depth_{depth}/reads/sample{n:02d}_R2.fastq.gz"
            truth = MULTI / f"depth_{depth}/truth/sample{n:02d}.truth.tsv"
            if all(p.exists() for p in (r1, r2, truth)):
                rows.append({"kind": "multi", "sample": sample, "species": "multi",
                             "r1": r1.relative_to(PAPER).as_posix(),
                             "r2": r2.relative_to(PAPER).as_posix(),
                             "truth": truth.relative_to(PAPER).as_posix()})
    return rows

def host_to_container(path: str | Path) -> str:
    path = Path(path)
    return str(path if path.is_absolute() else PAPER / path).replace(str(PAPER), PAPERDIR, 1)

def build(species: list[str], threads: int) -> None:
    ensure_container()
    WORK.joinpath("db").mkdir(parents=True, exist_ok=True)
    WORK.joinpath("build_logs").mkdir(parents=True, exist_ok=True)
    for sp in species:
        db = WORK / "db" / sp
        if (db / "Cluster_Result/hclsMap_95_recls.txt").exists():
            print(f"[skip db] {sp}")
            continue
        args = docker_args(PYBIN, "StrainScan_build.py",
                           "-i", f"{PAPERDIR}/figure_raw_data/sim_genome_pool/{sp}",
                           "-o", f"{PAPERDIR}/work/strainscan_headtohead_rerun/db/{sp}",
                           "-t", str(threads))
        sh(args, log=WORK / "build_logs" / f"{sp}.log")
        if not (db / "Cluster_Result/hclsMap_95_recls.txt").exists():
            raise RuntimeError(f"StrainScan DB build did not create cluster mapping for {sp}")

def run(rows: list[dict[str, str]]) -> None:
    ensure_container()
    failed = []
    for row in rows:
        out = WORK / "results" / row["kind"] / row["sample"]
        final = out / "final_report.txt"
        if final.exists() and final.stat().st_size:
            print(f"[skip profile] {row['sample']}")
            continue
        args = docker_args(PYBIN, "StrainScan.py",
                           "-i", host_to_container(row["r1"]),
                           "-j", host_to_container(row["r2"]),
                           "-d", f"{PAPERDIR}/work/strainscan_headtohead_rerun/db/{row['species']}",
                           "-o", f"{PAPERDIR}/work/strainscan_headtohead_rerun/results/{row['kind']}/{row['sample']}")
        # StrainScan's low-depth Layer-1 mode is required for 0.5x and 1x samples.
        if row["kind"] == "single" and (row["sample"].endswith("_d0.5") or row["sample"].endswith("_d1")):
            args.extend(["-l", "1"])
        result = subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if result.returncode:
            failed.append(row["sample"])
            print(f"[failed] {row['sample']} (returncode {result.returncode})", flush=True)
            continue
        if not final.exists():
            failed.append(row["sample"])
            print(f"[failed] {row['sample']} (no final report)", flush=True)
            continue
        print(f"[done] {row['sample']}", flush=True)
    (WORK / "failed_runs.txt").write_text("\n".join(failed) + ("\n" if failed else ""))
    print(f"profile completed: {len(rows)-len(failed)}/{len(rows)}; failed: {len(failed)}")

def write_manifest(rows: list[dict[str, str]]) -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    with (WORK / "run_manifest.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["kind", "sample", "species", "r1", "r2", "truth"], delimiter="\t")
        writer.writeheader(); writer.writerows(rows)

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--profile", action="store_true")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--species", action="append", help="restrict DB build/run to this species")
    ap.add_argument("--limit", type=int, default=0, help="profile only the first N manifest rows")
    args = ap.parse_args()
    if not args.build and not args.profile:
        ap.error("choose --build and/or --profile")
    cfg = json.loads(CONFIG.read_text())
    if not _container_id().startswith(cfg["comparators"]["StrainScan"]["container_id"]):
        print("warning: StrainScan container ID differs from benchmark_configuration.json", file=sys.stderr)
    rows = manifest()
    write_manifest(rows)
    if args.species:
        wanted = set(args.species)
        if args.profile:
            # For a smoke run, include multi only when explicitly requested.
            rows = [r for r in rows if r["species"] in wanted]
        else:
            rows = []
    if args.limit:
        rows = rows[:args.limit]
    else:
        # Preserve the published matched subset: diff-cluster rep1 for all depths across species,
        # same-cluster rep1 for near-clonal M. tuberculosis, plus multi samples 01-04 per depth.
        matched = []
        for row in rows:
            if row["kind"] == "multi":
                sample = row["sample"].rsplit("_", 1)[-1]
                if sample in {f"{n:02d}" for n in range(1, 5)}:
                    matched.append(row)
            elif "_rep1_" in row["sample"]:
                species = row["species"]
                if species == "Mycobacterium_tuberculosis":
                    if "__same_k" in row["sample"]:
                        matched.append(row)
                elif "__diff_k" in row["sample"]:
                    matched.append(row)
        rows = matched
    print(f"manifest: {len(rows)} runs")
    if args.build:
        build(species_list() if not args.species else list(args.species), args.threads)
    if args.profile and rows:
        run(rows)
    return 0

def _container_id() -> str:
    result = subprocess.run(["docker", "inspect", "-f", "{{.Id}}", CONTAINER], capture_output=True, text=True)
    return result.stdout.strip().removeprefix("sha256:")

if __name__ == "__main__":
    raise SystemExit(main())
