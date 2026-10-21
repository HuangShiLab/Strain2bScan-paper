#!/usr/bin/env python3
"""Summarize the uniform full-depth 32-sample saliva WMS rebuild against native calls."""
from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WMS = ROOT / "results/wms32_rebuild/wms32_strain_calls_all.tsv"
NATIVE = ROOT / "results/wms32_rebuild/native20_matched/native_strain_calls.tsv"
OUTDIR = ROOT / "results/wms32_rebuild"


def sample_label(raw: str) -> str:
    """Convert WMS_S11_11 or S11-11 to the shared S11-11 label."""
    return raw.removeprefix("WMS_").replace("_", "-")


def read_calls(path: Path, sample_column: str) -> dict[str, set[tuple[str, str]]]:
    calls: dict[str, set[tuple[str, str]]] = defaultdict(set)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            calls[sample_label(row[sample_column])].add((row["species"], row["cluster"]))
    return calls


wms = read_calls(WMS, "#sample")
native = read_calls(NATIVE, "#sample")
if set(wms) != set(native) or len(wms) != 32:
    raise SystemExit(f"sample mismatch: WMS={len(wms)}, native={len(native)}")

rows = [["sample", "wms_calls", "native_calls", "shared_calls",
         "wms_only_calls", "native_only_calls", "wms_species"]]
for sample in sorted(wms):
    w, n = wms[sample], native[sample]
    shared = w & n
    rows.append([
        sample, len(w), len(n), len(shared), len(w - n), len(n - w),
        len({species for species, _ in w}),
    ])

with (OUTDIR / "per_sample_summary.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerows(rows)

wms_calls = sum(int(row[1]) for row in rows[1:])
native_calls = sum(int(row[2]) for row in rows[1:])
shared_calls = sum(int(row[3]) for row in rows[1:])
summary = {
    "n_samples": len(rows) - 1,
    "wms_calls": wms_calls,
    "median_wms_calls_per_sample": statistics.median(int(row[1]) for row in rows[1:]),
    "native_calls": native_calls,
    "median_native_calls_per_sample": statistics.median(int(row[2]) for row in rows[1:]),
    "wms_calls_shared_with_native": shared_calls,
    "wms_calls_shared_fraction": shared_calls / wms_calls,
    "native_calls_shared_with_wms": shared_calls,
    "native_calls_shared_fraction": shared_calls / native_calls,
}
(OUTDIR / "full_wms_profile_summary.json").write_text(
    json.dumps(summary, indent=2) + "\n"
)

provenance = {
    "source_commit": "6ae760c2f3bd2f39efc0e6aa06c53a84366cb16e",
    "panel_species": 19,
    "enzyme": "BcgI",
    "thresholds": {
        "min_species_markers": 20,
        "min_species_detect": 2,
        "min_support": 2,
        "min_coverage": 0.01,
        "min_abundance": 0,
    },
    "wms_manifest": "results/wms32_rebuild/MANIFEST.tsv",
    "native_source": "results/wms32_rebuild/native20_matched/native_strain_calls.tsv",
    "comparison_key": ["species", "cluster"],
    "interpretation": "modality and sequencing-depth comparison, not ground-truth validation",
}
(OUTDIR / "profiling_provenance.json").write_text(
    json.dumps(provenance, indent=2) + "\n"
)
