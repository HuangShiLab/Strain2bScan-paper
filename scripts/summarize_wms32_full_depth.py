#!/usr/bin/env python3
"""Summarize the uniform full-depth 32-sample saliva WMS rebuild against native calls."""
from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parents[1]
WMS = ROOT / "results/wms32_rebuild/wms32_strain_calls_all.tsv"
NATIVE = ROOT / "results/wms32_rebuild/native20_matched/native_strain_calls.tsv"
OUTDIR = ROOT / "results/wms32_rebuild"


def sample_label(raw: str) -> str:
    """Convert WMS_S11_11 or S11-11 to the shared S11-11 label."""
    return raw.removeprefix("WMS_").replace("_", "-")


def read_calls(path: Path, sample_column: str) -> dict[str, dict[tuple[str, str], float]]:
    calls: dict[str, dict[tuple[str, str], float]] = defaultdict(dict)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["species"], row["cluster"])
            calls[sample_label(row[sample_column])][key] = float(row["sample_fraction"])
    return calls


wms_data = read_calls(WMS, "#sample")
native_data = read_calls(NATIVE, "#sample")
wms = {sample: set(data) for sample, data in wms_data.items()}
native = {sample: set(data) for sample, data in native_data.items()}
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

native_rows = []
per_sample_medians = []
for sample in sorted(native_data):
    shared_values = []
    native_only_values = []
    for key, sample_fraction in native_data[sample].items():
        group = "shared" if key in wms_data[sample] else "native_only"
        native_rows.append({"sample": sample, "group": group,
                            "sample_fraction": sample_fraction})
        (shared_values if group == "shared" else native_only_values).append(sample_fraction)
    per_sample_medians.append({
        "sample": sample,
        "shared_median_sample_fraction": statistics.median(shared_values),
        "native_only_median_sample_fraction": statistics.median(native_only_values),
    })

with (OUTDIR / "native_call_shared_status.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["sample", "group", "sample_fraction"],
                            delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(native_rows)

with (OUTDIR / "per_sample_shared_native_median.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(per_sample_medians[0]),
                            delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(per_sample_medians)

def _median_group(group: str) -> float:
    values = [row["sample_fraction"] for row in native_rows if row["group"] == group]
    return statistics.median(values)

p_value = wilcoxon([row["shared_median_sample_fraction"] for row in per_sample_medians],
                   [row["native_only_median_sample_fraction"] for row in per_sample_medians]).pvalue
abundance_summary = [
    {"group": "shared", "n_calls": sum(row["group"] == "shared" for row in native_rows),
     "median_sample_fraction": _median_group("shared")},
    {"group": "native_only", "n_calls": sum(row["group"] == "native_only" for row in native_rows),
     "median_sample_fraction": _median_group("native_only")},
    {"group": "paired_per_sample_median_wilcoxon_p", "n_calls": len(per_sample_medians),
     "median_sample_fraction": p_value},
]
with (OUTDIR / "shared_native_abundance_summary.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["group", "n_calls", "median_sample_fraction"],
                            delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(abundance_summary)

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
