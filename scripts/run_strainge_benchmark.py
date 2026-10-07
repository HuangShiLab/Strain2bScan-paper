#!/usr/bin/env python3
"""Run and score the StrainGE comparator rerun.

The benchmark uses the StrainGE-recommended reference-database workflow:
k-merize all panel genomes, remove near-subsets, cluster at Jaccard 0.90,
then create one species-level StrainGST database. Samples are k-merized once
and searched against each relevant species database. Accuracy is reported in
the StrainGE 0.90-reference space; no Strain2bScan clusters are used.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/strainge_rerun"
OUT = ROOT / "results/strainge_rerun"
DB = WORK / "databases"
RUNS = WORK / "runs"
DEFAULT_STRAINGE = Path("/Users/macstudio/Downloads/Strain2bLong/tools/envs/strainge/bin/straingst")
TIME_RE = re.compile(r"^\s*([0-9.]+) real\s+", re.M)
RSS_RE = re.compile(r"^\s*([0-9]+)\s+maximum resident set size\s+", re.M)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run_logged(cmd: list[str], log: Path, force: bool = False) -> dict[str, float | int]:
    log.parent.mkdir(parents=True, exist_ok=True)
    if log.exists() and TIME_RE.search(log.read_text(errors="replace")) and not force:
        text = log.read_text(errors="replace")
        wall = float(TIME_RE.search(text).group(1)) if TIME_RE.search(text) else 0.0
        rss = int(RSS_RE.search(text).group(1)) if RSS_RE.search(text) else 0
        return {"wall_s": wall, "maxrss_bytes": rss, "cached": 1}
    started = time.perf_counter()
    with log.open("w") as fh:
        proc = subprocess.run(["/usr/bin/time", "-l", *cmd], stdout=fh, stderr=subprocess.STDOUT)
    wall = time.perf_counter() - started
    if proc.returncode:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}; see {log}")
    text = log.read_text(errors="replace")
    match = RSS_RE.search(text)
    return {"wall_s": wall, "maxrss_bytes": int(match.group(1)) if match else 0, "cached": 0}


def parse_clusters(path: Path) -> dict[str, str]:
    genome_to_rep = {}
    with path.open() as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 1:
                continue
            rep = fields[0]
            for genome in fields:
                label = Path(genome).name
                if label.endswith(".hdf5"):
                    label = label[:-5]
                if label:
                    genome_to_rep[label] = rep
    return genome_to_rep


def write_json(path: Path, obj: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build_one_db(name: str, hdf5s: list[Path], dbdir: Path, scope: str,
                 n_expected: int | None = None) -> dict[str, object]:
    dbdir.mkdir(parents=True, exist_ok=True)
    sim_path = dbdir / f"{name}.kmersim.tsv"
    cluster_path = dbdir / f"{name}.clusters.tsv"
    keep_path = dbdir / f"{name}.keep.txt"
    keep_hdf5_path = dbdir / f"{name}.keep.hdf5.txt"
    db_path = dbdir / f"{name}.straingst.hdf5"
    if not (sim_path.exists() and cluster_path.exists() and keep_path.exists() and keep_hdf5_path.exists() and db_path.exists()):
        timing_sim = run_logged(
            [str(DEFAULT_STRAINGE), "kmersim", "--all-vs-all", "-t", "4", "-S", "jaccard", "-S", "subset",
             "-o", str(sim_path), *map(str, hdf5s)], dbdir / "logs" / f"{name}.kmersim.log")
        timing_cluster = run_logged(
            [str(DEFAULT_STRAINGE), "cluster", "-i", str(sim_path), "-d", "-C", "0.99", "-c", "0.90",
             "--clusters-out", str(cluster_path), "-o", str(keep_path), *map(str, hdf5s)],
            dbdir / "logs" / f"{name}.cluster.log")
        timing_db = run_logged(
            [str(DEFAULT_STRAINGE), "createdb", "-f", str(keep_path), "-o", str(db_path)],
            dbdir / "logs" / f"{name}.createdb.log")
    else:
        timing_sim = run_logged([str(DEFAULT_STRAINGE)], dbdir / "logs" / f"{name}.kmersim.log")
        timing_cluster = run_logged([str(DEFAULT_STRAINGE)], dbdir / "logs" / f"{name}.cluster.log")
        timing_db = run_logged([str(DEFAULT_STRAINGE)], dbdir / "logs" / f"{name}.createdb.log")
    kept = [x.strip() for x in keep_path.read_text().splitlines() if x.strip()]
    clusters = parse_clusters(cluster_path)
    return {
        "scope": scope,
        "species": name,
        "n_input_genomes": len(hdf5s),
        "n_skipped_empty_genomes": max(0, (n_expected or len(hdf5s)) - len(hdf5s)),
        "n_kept_references": len(kept),
        "n_clustered_genomes": len(clusters),
        "build_wall_s": sum(float(x["wall_s"]) for x in (timing_sim, timing_cluster, timing_db)),
        "build_maxrss_bytes": max(int(x["maxrss_bytes"]) for x in (timing_sim, timing_cluster, timing_db)),
        "database": str(db_path.relative_to(ROOT)),
        "database_sha256": sha256(db_path),
        "clusters": str(cluster_path.relative_to(ROOT)),
        "keep_list": str(keep_path.relative_to(ROOT)),
    }


def build_sim_databases() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    source = ROOT / "figure_raw_data/sim_genome_pool"
    groups: dict[str, list[Path]] = defaultdict(list)
    expected: dict[str, int] = {}
    for fasta in sorted(source.glob("*/*.fna")):
        hdf5 = DB / "kmersets" / f"{fasta.stem}.hdf5"
        species = fasta.parent.name
        expected[species] = expected.get(species, 0) + 1
        if not hdf5.exists():
            if fasta.stat().st_size <= 1024:
                continue
            raise FileNotFoundError(f"missing k-mer set: {hdf5}; run kmerize_one.sh first")
        groups[fasta.parent.name].append(hdf5)
    skipped = [(species, expected[species] - len(groups[species])) for species in groups]
    rows = [build_one_db(species, paths, DB / "sim", "simulation_15_species")
            for species, paths in sorted(groups.items())]
    for row, (_, count) in zip(rows, sorted(groups.items())):
        row["n_skipped_empty_genomes"] = expected[row["species"]] - len(groups[row["species"]])
    with (OUT / "skipped_empty_simulation_genomes.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["species", "n_skipped_empty_genomes"])
        writer.writerows([(species, count) for species, count in skipped if count])
    fields = list(rows[0])
    with (OUT / "simulation_database_build.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def build_mock_databases() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    source = Path("/Users/macstudio/Downloads/Strain2bScan-raw-data/MSA_all164")
    groups: dict[str, list[Path]] = defaultdict(list)
    for fasta in sorted(source.glob("*.fna")):
        species = fasta.name.split("__", 1)[0]
        hdf5 = DB / "mock_kmersets" / f"{fasta.stem}.hdf5"
        if not hdf5.exists():
            timing = run_logged(
                [str(DEFAULT_STRAINGE), "kmerize", "-o", str(hdf5), str(fasta)],
                DB / "logs" / f"mock_kmerize_{fasta.stem}.log")
        else:
            timing = {"wall_s": 0, "maxrss_bytes": 0}
        groups[species].append(hdf5)
    rows = [build_one_db(species, paths, DB / "mock", "ATCC_DNA_mock_164") for species, paths in sorted(groups.items())]
    fields = list(rows[0]) + ["kmerize_wall_s", "kmerize_maxrss_bytes"]
    for row in rows:
        row["kmerize_wall_s"] = sum(
            float(run_logged([str(DEFAULT_STRAINGE)], DB / "logs" / f"mock_kmerize_{p.stem}.log")["wall_s"])
            for p in groups[row["species"]]
        )
    with (OUT / "mock_database_build.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def selected_simulation_manifest() -> list[dict[str, str]]:
    rows = read_tsv(ROOT / "results/strainscan_rerun/run_manifest.tsv")
    return [row for row in rows if row["kind"] == "single" and "_rep1_" in row["sample"] or row["kind"] == "multi"]


def sample_kmerset(sample: str, r1: Path, r2: Path, force: bool = False) -> tuple[Path, dict[str, float | int]]:
    outdir = RUNS / "sample_kmersets" / sample
    outdir.mkdir(parents=True, exist_ok=True)
    # StrainGE infers FASTQ only from .fastq/.fastq.gz. Some ATCC raw files use
    # *_fq.gz, so expose stable .fastq.gz aliases without copying large reads.
    read1 = r1 if r1.name.endswith((".fastq", ".fastq.gz")) else outdir / "read1.fastq.gz"
    read2 = r2 if r2.name.endswith((".fastq", ".fastq.gz")) else outdir / "read2.fastq.gz"
    if read1 != r1:
        read1.unlink(missing_ok=True)
        read1.symlink_to(r1)
    if read2 != r2:
        read2.unlink(missing_ok=True)
        read2.symlink_to(r2)
    sample_h5 = outdir / "sample.hdf5"
    log = outdir / "kmerize.log"
    needs_kmerize = force or not sample_h5.exists()
    if needs_kmerize:
        sample_h5.unlink(missing_ok=True)
        timing = run_logged(
            [str(DEFAULT_STRAINGE), "kmerize", "-o", str(sample_h5), str(read1), str(read2)],
            log, force=True)
    else:
        timing = run_logged([str(DEFAULT_STRAINGE)], log)
    return sample_h5, timing


def parse_time_log(log: Path) -> dict[str, float | int]:
    text = log.read_text(errors="replace")
    wall = float(TIME_RE.search(text).group(1)) if TIME_RE.search(text) else 0.0
    match = RSS_RE.search(text)
    return {"wall_s": wall, "maxrss_bytes": int(match.group(1)) if match else 0}


def run_simulation_manifest() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = selected_simulation_manifest()
    db_build = {row["species"]: row for row in read_tsv(OUT / "simulation_database_build.tsv")}
    output_rows = []
    for row in manifest:
        sample = row["sample"]
        r1, r2 = ROOT / row["r1"], ROOT / row["r2"]
        outdir = RUNS / "simulation" / sample
        outdir.mkdir(parents=True, exist_ok=True)
        sample_h5, kmer_timing = sample_kmerset(sample, r1, r2)
        species_names = [row["species"]] if row["kind"] == "single" else sorted(db_build)
        search_wall = 0.0
        maxrss = int(kmer_timing["maxrss_bytes"])
        status = "completed"
        for species in species_names:
            db_path = ROOT / db_build[species]["database"]
            log = outdir / f"straingst_{species}.log"
            result = outdir / f"straingst_{species}.tsv"
            iterations = "5" if row["kind"] == "single" else "8"
            try:
                timing = run_logged(
                    [str(DEFAULT_STRAINGE), "run", str(db_path), str(sample_h5), "-o", str(result), "-i", iterations],
                    log, force=not result.exists())
                search_wall += float(timing["wall_s"])
                maxrss = max(maxrss, int(timing["maxrss_bytes"]))
            except Exception as exc:
                print(exc, file=sys.stderr)
                status = "failed"
        output_rows.append({
            "benchmark": "simulation", "kind": row["kind"], "sample": sample,
            "species": row["species"], "depth": sample.rsplit("_d", 1)[-1] if row["kind"] == "single" else sample.split("_")[1],
            "n_databases": len(species_names), "kmerize_wall_s": float(kmer_timing["wall_s"]),
            "search_wall_s_sum": search_wall, "total_wall_s": float(kmer_timing["wall_s"]) + search_wall,
            "maxrss_bytes": maxrss, "status": status,
            "output_dir": str(outdir.relative_to(ROOT)),
        })
        write_json(outdir / "efficiency.json", output_rows[-1])
        sample_h5.unlink(missing_ok=True)
    fields = list(output_rows[0])
    with (OUT / "simulation_run_efficiency.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output_rows)


MOCKS = [
    ("MSA1002_99", "MSA1002", "/Users/macstudio/Downloads/Strain2bScan-raw-data/MSA1002/shotgun/WMS_MSA1002_99_100ng_1_R1_fq.gz", "/Users/macstudio/Downloads/Strain2bScan-raw-data/MSA1002/shotgun/WMS_MSA1002_99_100ng_1_R2_fq.gz"),
    ("MSA1003_0", "MSA1003", "/Users/macstudio/Downloads/Strain2bScan-raw-data/MSA1003/shotgun/WMS_MSA1003_0_100ng_1_R1_fq.gz", "/Users/macstudio/Downloads/Strain2bScan-raw-data/MSA1003/shotgun/WMS_MSA1003_0_100ng_1_R2_fq.gz"),
    ("MSA1005_0", "MSA1005", "/Users/macstudio/Downloads/Strain2bScan-raw-data/MSA1005/shotgun/WMS_MSA1005_0_100ng_1_R1_fq.gz", "/Users/macstudio/Downloads/Strain2bScan-raw-data/MSA1005/shotgun/WMS_MSA1005_0_100ng_1_R2_fq.gz"),
    ("MSA1007_0", "MSA1007", "/Users/macstudio/Downloads/Strain2bScan-raw-data/MSA1007/shotgun/WMS_MSA1007_0_100ng_1_R1_fq.gz", "/Users/macstudio/Downloads/Strain2bScan-raw-data/MSA1007/shotgun/WMS_MSA1007_0_100ng_1_R2_fq.gz"),
]


def run_mock_manifest() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    db_build = {row["species"]: row for row in read_tsv(OUT / "mock_database_build.tsv")}
    output_rows = []
    for sample, mock, r1, r2 in MOCKS:
        outdir = RUNS / "mock" / sample
        outdir.mkdir(parents=True, exist_ok=True)
        sample_h5, kmer_timing = sample_kmerset(sample, Path(r1), Path(r2))
        search_wall = 0.0
        maxrss = int(kmer_timing["maxrss_bytes"])
        status = "completed"
        for species, metadata in sorted(db_build.items()):
            db_path = ROOT / metadata["database"]
            log = outdir / f"straingst_{species}.log"
            result = outdir / f"straingst_{species}.tsv"
            try:
                timing = run_logged(
                    [str(DEFAULT_STRAINGE), "run", str(db_path), str(sample_h5), "-o", str(result), "-i", "32"],
                    log, force=not result.exists())
                search_wall += float(timing["wall_s"])
                maxrss = max(maxrss, int(timing["maxrss_bytes"]))
            except Exception as exc:
                print(exc, file=sys.stderr)
                status = "failed"
        output_rows.append({
            "benchmark": "ATCC_DNA_mock", "sample": sample, "mock": mock,
            "n_databases": len(db_build), "kmerize_wall_s": float(kmer_timing["wall_s"]),
            "search_wall_s_sum": search_wall, "total_wall_s": float(kmer_timing["wall_s"]) + search_wall,
            "maxrss_bytes": maxrss, "status": status,
            "output_dir": str(outdir.relative_to(ROOT)),
        })
        write_json(outdir / "efficiency.json", output_rows[-1])
        sample_h5.unlink(missing_ok=True)
    fields = list(output_rows[0])
    with (OUT / "mock_run_efficiency.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output_rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("step", choices=["build-sim-db", "build-mock-db", "run-sim", "run-mock"])
    args = parser.parse_args()
    if args.step == "build-sim-db":
        build_sim_databases()
    elif args.step == "build-mock-db":
        build_mock_databases()
    elif args.step == "run-sim":
        run_simulation_manifest()
    elif args.step == "run-mock":
        run_mock_manifest()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
