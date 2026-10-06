#!/usr/bin/env python3
"""Check manuscript claims against frozen benchmark evidence and tracked tables.

This is a release gate for the primary ATCC mock, runtime, saliva, and cohort claims.
It checks source-derived numbers, forbidden stale phrasing, required qualification,
artifact checksums, figure-manifest checksums, and the primary plotting configuration.
"""
import csv
import hashlib
import json
import math
import re
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
errors = []

def fail(msg):
    errors.append(msg)

def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()

def load_tsv(path):
    with (ROOT / path).open() as fh:
        lines = [line for line in fh if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))

def read(*parts):
    return "\n".join((ROOT / p).read_text() for p in parts)

# 1. Frozen configuration and artifact integrity.
full_manuscript = (ROOT / "manuscript/full_manuscript.md").read_text()
if full_manuscript.count("\n## Results\n") != 1:
    fail("assembled manuscript must contain exactly one top-level Results heading")
if re.search(r"^##\s+Results\s*\n+###\s+Results\s*$", full_manuscript, flags=re.M):
    fail("assembled manuscript contains a duplicated nested Results heading")
for heading in ("## Introduction", "## Results", "## Discussion", "## Methods",
                "## Figure legends", "## Tables", "## References"):
    if full_manuscript.count("\n" + heading + "\n") != 1:
        fail(f"assembled manuscript must contain exactly one {heading!r} heading")
config_path = ROOT / "results/benchmark_configuration.json"
try:
    cfg = json.loads(config_path.read_text())
except Exception as exc:
    fail(f"cannot parse benchmark configuration: {exc}")
    cfg = {}
if cfg.get("configuration_version") != "mock-primary-v1":
    fail("unexpected benchmark configuration version")
for item in cfg.get("derived_evidence", []):
    p = item["path"]
    if not (ROOT / p).exists():
        fail(f"missing frozen evidence: {p}")
    elif sha(p) != item["sha256"]:
        fail(f"checksum mismatch: {p}")
release = cfg.get("strain2bscan", {})
release_archive = release.get("source_archive")
if not release_archive or not (ROOT / release_archive).exists():
    fail("benchmark source archive is missing")
elif sha(release_archive) != release.get("source_archive_sha256"):
    fail("benchmark source archive checksum mismatch")
if release.get("benchmark_tag") != "benchmark-v0.1.0-f26f234":
    fail("unexpected benchmark release tag")
manifest_path = ROOT / "results/mock_benchmark_f26f234/checksums.sha256"
if manifest_path.exists():
    for line in manifest_path.read_text().splitlines():
        if not line.strip():
            continue
        expected, rel = line.split("  ", 1)
        prefix = "results/mock_benchmark_f26f234/"
        if rel.startswith(prefix):
            rel = rel[len(prefix):]
        p = ROOT / "results/mock_benchmark_f26f234" / rel
        if not p.exists():
            fail(f"missing frozen mock artifact: {rel}")
        elif hashlib.sha256(p.read_bytes()).hexdigest() != expected:
            fail(f"checksum mismatch in frozen mock artifact: {rel}")

# 2. Figure manifest integrity.
fig_rows = load_tsv("figures/numbered/MANIFEST.tsv")
for row in fig_rows:
    for ext in ("png", "pdf"):
        p = Path("figures/numbered") / row[ext]
        if not p.exists():
            fail(f"missing numbered figure: {p}")
        elif sha(p) != row[ext + "_sha256"]:
            fail(f"figure checksum mismatch: {p}")

# 3. Recompute source-table ranges used in the revised text.
mock = load_tsv("data/fig6_fig12_metrics.tsv")
def rows(kind=None, mock_name=None, tool=None, variant=None):
    out = []
    for r in mock:
        if kind and r["kind"] != kind: continue
        if mock_name and r["mock"] != mock_name: continue
        if tool and r["tool"] != tool: continue
        if variant and r["variant"] != variant: continue
        out.append(r)
    return out

native = rows("2bRAD", "MSA1002", "Strain2bScan", "164")
low01 = [r for r in native if r["sample"].endswith("0_0.01ng_1")]
low1 = [r for r in native if r["sample"].endswith("0_0.1ng_1")]
host99 = [r for r in rows("WMS", "MSA1002", "Strain2bScan", "164") if r["sample"].endswith("99_100ng_1")]
clean_msa1003 = rows("WMS", "MSA1003", "Strain2bScan", "164")
clean_msa1005 = rows("WMS", "MSA1005", "Strain2bScan", "164")
clean_msa1007 = rows("WMS", "MSA1007", "Strain2bScan", "164")
if len(low01) != 1 or len(low1) != 1 or len(host99) != 1:
    fail("expected exactly one primary native MSA1002 row at 0.01 ng, 0.1 ng, and 99% host")
else:
    expected = {
        "P@1e-4 0.01 ng": float(low01[0]["precision_at_1e_4"]),
        "R@1e-4 0.01 ng": float(low01[0]["recall_at_1e_4"]),
        "F1@1e-4 0.01 ng": float(low01[0]["f1_at_1e_4"]),
        "AUPR 0.01 ng": float(low01[0]["aupr"]),
        "F1@1e-4 0.1 ng": float(low1[0]["f1_at_1e_4"]),
        "F1@1e-4 99% host": float(host99[0]["f1_at_1e_4"]),
        "BC similarity 99% host": 1 - float(host99[0]["bray_curtis"]),
    }
    wanted = {"P@1e-4 0.01 ng": 0.833, "R@1e-4 0.01 ng": 0.500,
              "F1@1e-4 0.01 ng": 0.625, "AUPR 0.01 ng": 0.500,
              "F1@1e-4 0.1 ng": 0.952, "F1@1e-4 99% host": 1.000,
              "BC similarity 99% host": 0.698}
    for k, v in wanted.items():
        if not math.isclose(expected[k], v, abs_tol=5e-4):
            fail(f"source-table mismatch for {k}: got {expected[k]:.4g}, text rule expects {v:.4g}")

for name, rr, lo, hi in [
    ("MSA1003", clean_msa1003, 0.90, 1.00),
    ("MSA1005", clean_msa1005, 0.60, 0.65),
    ("MSA1007", clean_msa1007, 0.70, 0.76),
]:
    f1 = [float(r["f1_at_1e_4"]) for r in rr]
    if len(f1) != 3:
        fail(f"expected 3 primary {name} WMS rows, got {len(f1)}")
    elif min(f1) < lo or max(f1) > hi:
        fail(f"primary {name} WMS F1@1e-4 outside claimed {lo:.2f}-{hi:.2f}: {f1}")

scaling = load_tsv("results/parallel_and_build_scaling.tsv")
speed = {r["metric"]: float(r["speedup"]) for r in scaling if r.get("metric") in {"db_build", "profile_sample4"}}
if not (math.isclose(speed["db_build"], 4.6, abs_tol=0.05) and math.isclose(speed["profile_sample4"], 5.8, abs_tol=0.05)):
    fail(f"scaling table changed: {speed}")
throughput = load_tsv("results/community_throughput.tsv")
tp = [float(r["speedup_x"]) for r in throughput]
if min(tp) != 121 or max(tp) != 146:
    fail(f"community throughput range changed: {min(tp)}-{max(tp)}")
for r in throughput:
    if "projected" not in r.get("strainscan_projected_s", "").lower() and "projected" not in " ".join(r.keys()).lower():
        # The label is in the comments, so absence here is not an error, but text must say projected.
        pass

saliva_ml = load_tsv("results/saliva_temporal_ml.tsv")
ml = {r["metric"]: r["strain"] for r in saliva_ml}
if ml.get("leave_one_timepoint_out_subject_accuracy") != "1.0000" or ml.get("n_subjects") != "8":
    fail(f"unexpected saliva temporal ML table: {ml}")
concordance = load_tsv("results/saliva_concordance.tsv")
if len(concordance) != 3 or sum(int(r["wms_confirmed_by_2brad"]) for r in concordance) != 65:
    fail("saliva concordance table no longer supports the three-sample 65/65 statement")
isolate = load_tsv("results/realworld_cohort_benchmark/isolate_validation_scorecard.tsv")
if len(isolate) != 6 or sum(r["correct"] == "True" for r in isolate) != 6:
    fail("isolate scorecard no longer supports the 6/6 self-call statement")

# 4. Primary plotting configuration.
plot = (ROOT / "scripts/plot_figs_h.py").read_text()
wms_start = plot.index("def build_wms():")
wms_end = plot.index("# ---------------- assemble Fig 6")
wms_block = plot[wms_start:wms_end]
brad_start = plot.index("def build_brad():")
brad_end = plot.index("# ---------------- assemble Fig S")
brad_block = plot[brad_start:brad_end]
if "Strain2bScan-port" in wms_block:
    fail("primary Figure 12 includes non-primary Strain2bScan-port rows")
if '"120"' in brad_block:
    fail("primary Figure 6 includes 20-species tree rows")
if '"120"' in wms_block:
    fail("primary Figure 12 includes 20-species tree rows")
if "164-tracegap" in wms_block:
    fail("primary Figure 12 includes trace-gap rows")

# 5. Required and forbidden manuscript phrasing.
text = read("manuscript/abstract.md", "manuscript/results.md", "manuscript/methods.md",
            "manuscript/discussion.md", "manuscript/figures.md",
            "manuscript/figure_legends_detailed.md")
required = [
    "250 bp", "ART", "error-modelled", "F1 at 1e-4 = 0.909", "AUPR = 1.0",
    "candidate strain-cluster calls", "panel/read compatibility", "prefix subsample",
    "121–146", "projected", "4.6×", "5.8×", "8 subjects",
]
for phrase in required:
    if phrase not in text:
        fail(f"required qualified claim absent from modular manuscript: {phrase!r}")
forbidden = [
    "error-free 150 bp", "150 bp paired-end", "Strain2bScan and StrainScan F1 = 1.0",
    "all three tools performed well", "confirmed by native", "130–146", "8.5×", "6.5×",
    "recovered all six", "all strains from\nthe ATCC mocks",
]
for phrase in forbidden:
    if phrase in text:
        fail(f"stale or unsupported phrasing still present: {phrase!r}")
for p in ("manuscript/abstract.md", "manuscript/results.md", "manuscript/discussion.md"):
    doc = (ROOT / p).read_text()
    for phrase in ("confirmed", "validated"):
        if phrase in doc:
            fail(f"unqualified validation language in {p}: {phrase!r}")

if errors:
    print("MANUSCRIPT CLAIM CHECK FAILED", file=sys.stderr)
    for e in errors:
        print(f" - {e}", file=sys.stderr)
    sys.exit(1)
print("Manuscript claim check passed: frozen artifacts, source-derived numbers, primary figures, and qualified claims are consistent.")
