#!/usr/bin/env python3
"""Verify synchronized rhesus macaque HPC outputs against manuscript cohort counts."""
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/monkey_spatial_148/hpc_outputs"

def rows(rel):
    with open(BASE / rel, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

meta = rows("analysis/full_148_20261027/full_148_sample_metadata.tsv")
summary = rows("analysis/full_148_20261027/full148_sample_strain_summary.tsv")
calls = rows("analysis/full_148_20261027/full148_2brad_strain_calls_long.tsv")
clusters = rows("analysis/full_148_20261027/full148_2brad_cluster_GI_control_region_summary.tsv")
paired = rows("analysis/full_148_20261027/full148_mucosal_luminal_paired_tests.tsv")
ba = rows("analysis/functional_coupling_20261027/coupling_compositional_bileacid_li.tsv")
functions = rows("analysis/functional_coupling_20261027/coupling_compositional_functions_li.tsv")
panel = rows("work/monkey_panel_20261023/selected_species_manifest.tsv")
db = rows("analysis/strain_profiling_20261023/panel_db_metrics.tsv")

def group(rowset, value, key="sample_group"):
    return [r for r in rowset if r.get(key) == value]

checks = {
    "metadata_rows": len(meta),
    "GI_samples": len(group(meta, "GI_spatial")),
    "body_controls": len(group(meta, "body_control")),
    "GI_positive": sum(r.get("has_2brad_strain_call") == "True" for r in group(summary, "GI_spatial")),
    "luminal_positive": sum(r.get("has_2brad_strain_call") == "True" for r in group(summary, "GI_spatial") if r.get("compartment") == "luminal"),
    "mucosal_positive": sum(r.get("has_2brad_strain_call") == "True" for r in group(summary, "GI_spatial") if r.get("compartment") == "mucosal"),
    "body_positive": sum(r.get("has_2brad_strain_call") == "True" for r in group(summary, "body_control")),
    "GI_call_rows": len(group(calls, "GI_spatial")),
    "body_call_rows": len(group(calls, "body_control")),
    "paired_tests": len(paired),
    "paired_bonferroni": sum(float(r["p_bonferroni"]) < 0.05 for r in paired),
    "LI_enriched_clusters": sum(
        r.get("peak_region") == "LI"
        and float(r.get("region_specificity", 0)) >= 0.8
        and float(r.get("max_region_prevalence", 0)) >= 0.3
        for r in clusters
    ),
    "BA_tests": len(ba),
    "BA_q_lt0.05": sum(float(r["q"]) < 0.05 for r in ba),
    "BA_negative_q_lt0.05": sum(float(r["q"]) < 0.05 and float(r["rho"]) < 0 for r in ba),
    "function_tests": len(functions),
    "function_q_lt0.05": sum(float(r["q"]) < 0.05 for r in functions),
    "function_positive_q_lt0.05": sum(float(r["q"]) < 0.05 and float(r["rho"]) > 0 for r in functions),
    "function_negative_q_lt0.05": sum(float(r["q"]) < 0.05 and float(r["rho"]) < 0 for r in functions),
    "panel_species": len(panel),
    "selected_MAGs": sum(int(r["n_mags"]) for r in panel),
    "cluster_units": sum(int(r["clusters"]) for r in db),
    "panel_markers": sum(int(r["markers"]) for r in db),
}

expected = {
    "metadata_rows": 148,
    "GI_samples": 128,
    "body_controls": 20,
    "GI_positive": 113,
    "luminal_positive": 56,
    "mucosal_positive": 57,
    "body_positive": 2,
    "GI_call_rows": 5285,
    "body_call_rows": 3,
    "paired_tests": 153,
    "paired_bonferroni": 13,
    "LI_enriched_clusters": 89,
    "BA_tests": 876,
    "BA_q_lt0.05": 410,
    "BA_negative_q_lt0.05": 409,
    "function_tests": 3066,
    "function_q_lt0.05": 1233,
    "function_positive_q_lt0.05": 1019,
    "function_negative_q_lt0.05": 214,
    "panel_species": 18,
    "selected_MAGs": 315,
    "cluster_units": 303,
    "panel_markers": 30008,
}
failed = [(key, checks[key], expected[key]) for key in expected if checks[key] != expected[key]]
out = ROOT / "results/monkey_spatial_148/hpc_output_verification.tsv"
with out.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["check", "value"])
    for key, value in checks.items():
        writer.writerow([key, value])

if failed:
    raise SystemExit("Monkey HPC output verification failed: " + "; ".join(
        f"{key}={actual}, expected={expected}" for key, actual, expected in failed))
print(f"Monkey HPC output verification passed ({len(checks)} checks): {out}")
