#!/usr/bin/env python3
"""Fill peak-RSS and CPU fields in StrainGE efficiency tables from /usr/bin/time logs."""
from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/strainge_rerun"
RUNS = ROOT / "work/strainge_rerun/runs"
DB = ROOT / "work/strainge_rerun/databases"
WALL = re.compile(r"^\s*([0-9.]+) real\s+", re.M)
RSS = re.compile(r"^\s*([0-9]+)\s+maximum resident set size\s+", re.M)
CPU = re.compile(r"^\s*[0-9.]+\s+real\s+([0-9.]+)\s+user\s+([0-9.]+)\s+sys\s*", re.M)


def parse_log(path: Path) -> dict[str, float | int]:
    text = path.read_text(errors="replace")
    wall = float(WALL.search(text).group(1)) if WALL.search(text) else 0.0
    rss = int(RSS.search(text).group(1)) if RSS.search(text) else 0
    cpu = CPU.search(text)
    return {
        "wall_s": wall,
        "maxrss_bytes": rss,
        "cpu_s": float(cpu.group(1)) + float(cpu.group(2)) if cpu else 0.0,
    }


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def repair_database(scope: str, dbroot: Path) -> None:
    path = OUT / f"{scope}_database_build.tsv"
    fields, rows = read_tsv(path)
    for field in ("build_cpu_s",):
        if field not in fields:
            fields.append(field)
    for row in rows:
        logs = [dbroot / "logs" / f'{row["species"]}.{stage}.log' for stage in ("kmersim", "cluster", "createdb")]
        parsed = [parse_log(log) for log in logs if log.exists()]
        row["build_wall_s"] = sum(item["wall_s"] for item in parsed)
        row["build_maxrss_bytes"] = max((item["maxrss_bytes"] for item in parsed), default=0)
        row["build_cpu_s"] = sum(item["cpu_s"] for item in parsed)
    write_tsv(path, fields, rows)


def repair_runs(scope: str) -> None:
    path = OUT / f"{scope}_run_efficiency.tsv"
    fields, rows = read_tsv(path)
    for field in ("kmerize_cpu_s", "search_cpu_s_sum", "total_cpu_s"):
        if field not in fields:
            fields.append(field)
    for row in rows:
        output_dir = ROOT / row["output_dir"]
        kmer_log = RUNS / "sample_kmersets" / row["sample"] / "kmerize.log"
        if kmer_log.exists():
            kmer = parse_log(kmer_log)
            if kmer["wall_s"] == 0:
                # A later interrupted rerun can truncate a valid timing log. The
                # first completed pass remains authoritative; its values are in
                # the runner-written efficiency JSON.
                backup_path = ROOT / row["output_dir"] / "efficiency.json"
                if backup_path.exists():
                    import json
                    backup = json.loads(backup_path.read_text())
                    kmer["wall_s"] = float(backup.get("kmerize_wall_s", 0))
                    kmer["maxrss_bytes"] = int(backup.get("maxrss_bytes", 0))
                    kmer["cpu_s"] = ""
        else:
            kmer = {"wall_s": 0, "maxrss_bytes": 0, "cpu_s": 0}
        searches = [parse_log(log) for log in sorted(output_dir.glob("straingst_*.log"))]
        row["kmerize_wall_s"] = kmer["wall_s"]
        row["kmerize_cpu_s"] = kmer["cpu_s"]
        row["search_wall_s_sum"] = sum(item["wall_s"] for item in searches)
        row["search_cpu_s_sum"] = sum(item["cpu_s"] for item in searches)
        row["total_wall_s"] = kmer["wall_s"] + sum(item["wall_s"] for item in searches)
        kmer_cpu = kmer["cpu_s"]
        row["total_cpu_s"] = (float(kmer_cpu) if kmer_cpu != "" else 0.0) + sum(item["cpu_s"] for item in searches)
        row["maxrss_bytes"] = max([int(kmer["maxrss_bytes"]), *[item["maxrss_bytes"] for item in searches]])
        expected_outputs = int(row.get("n_databases", len(searches)))
        actual_outputs = sum(path.exists() and path.stat().st_size > 0 for path in output_dir.glob("straingst_*.tsv"))
        row["status"] = "completed" if actual_outputs >= expected_outputs else "failed"
    write_tsv(path, fields, rows)


def main() -> int:
    repair_database("simulation", DB / "sim")
    repair_database("mock", DB / "mock")
    repair_runs("simulation")
    repair_runs("mock")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
