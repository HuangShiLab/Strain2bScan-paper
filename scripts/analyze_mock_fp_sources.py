#!/usr/bin/env python3
"""Classify false-positive sources in Strain2bScan mock outputs.

This supports the supplementary analysis of cross-library contamination:
most WMS false positives on MSA-1005/1007 are not algorithmic shadows but
other mocks' ATCC strains appearing at trace levels (index hopping).

Inputs:
  work/mock_retest/Strain2bScan-raw-data/wms_analysis/wms164/*.pred
  work/mock_retest/Strain2bScan-raw-data/Ground_truth/*_ground_truth.txt
  work/mock_retest/Strain2bScan-raw-data/MSA_combined164_all_flat.members.tsv

Outputs:
  results/mock_fp_source_summary.tsv
  figures/supplementary/mock_fp_source_bar.png (optional)
"""
import csv, json, os, re
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAPER = Path(__file__).resolve().parent.parent
RAW = PAPER / "work" / "mock_retest" / "Strain2bScan-raw-data"
RES = PAPER / "results"
FIG = PAPER / "figures"
OUT_TSV = RES / "mock_fp_source_summary.tsv"
OUT_JSON = RES / "mock_fp_source_profiles.json"

MEMBERS = {}
for ln in open(RAW / "MSA_combined164_all_flat.members.tsv"):
    if ln.startswith("#") or not ln.strip():
        continue
    g, c = ln.rstrip("\n").split("\t")[:2]
    MEMBERS.setdefault(c, []).append(g)

def rep_of(members):
    a = [g for g in members if "ATCC_" in g]
    return a[0] if a else sorted(members)[0]

def species_of(g):
    s = re.sub(r"__(ATCC_|GCF_).*", "", g)
    if s != g:
        return s
    m = re.match(r"ATCC_(.+?)_ATCC_", g)
    if m:
        return m.group(1)
    return g

def atcc_token(g):
    if "_ATCC_" in g:
        return "ATCC_" + g.rsplit("_ATCC_", 1)[-1]
    m = re.search(r"ATCC_[A-Za-z0-9_]+$", g)
    return m.group(0) if m else None

def load_truth(mock):
    f = RAW / "Ground_truth" / f"{mock}_ground_truth.txt"
    rows = [l.rstrip("\n").split("\t") for l in open(f) if l.strip()]
    hdr = rows[0]
    out = {}
    if mock == "MSA1003":
        for r in rows[1:]:
            if len(r) >= 2 and r[0].strip():
                out[r[0].strip()] = float(r[1])
    else:
        i = hdr.index("seq_abd") if "seq_abd" in hdr else 2
        for r in rows[1:]:
            if len(r) > i and r[0].strip():
                out[r[0].strip()] = float(r[i])
    return out

def truth_species(mock):
    t = load_truth(mock)
    out = set()
    for atcc in t:
        num = re.sub(r"^ATCC[_ ]*", "", atcc).strip()
        pat = rf"(?:^|_){re.escape(num)}(_|$)"
        for g in (g for gs in MEMBERS.values() for g in gs):
            if "ATCC_" in g and re.search(pat, g):
                out.add(species_of(g))
                break
    return out

def all_truth_species_per_mock():
    mocks = ["MSA1002", "MSA1003", "MSA1005", "MSA1007"]
    return {m: truth_species(m) for m in mocks}

def read_pred(pred):
    out = {}
    in_tsv = False
    for ln in open(pred):
        if ln.startswith("#cluster"):
            in_tsv = True
            continue
        if not in_tsv or not ln.strip():
            continue
        p = ln.rstrip("\n").split("\t")
        if len(p) < 2:
            continue
        c, ab = p[0], p[1]
        try:
            ab = float(ab)
        except ValueError:
            continue
        if c not in MEMBERS:
            continue
        r = rep_of(MEMBERS[c])
        out[r] = out.get(r, 0.0) + ab
    return out

def classify_fp(g, true_species, truth_by_mock):
    sp = species_of(g)
    tok = atcc_token(g)
    # same-species decoy: species is in this sample's truth, but genome differs
    if sp in true_species:
        return "same_species_decoy"
    # other-mock ATCC: this species is a true member of another mock
    for m, ts in truth_by_mock.items():
        if sp in ts:
            return "other_mock_atcc"
    # GCF decoy of an other-mock species (no ATCC token)
    if tok is None:
        for m, ts in truth_by_mock.items():
            if sp in ts:
                return "other_mock_gcf"
    return "other"

def analyze_sample(sample, mock, truth_by_mock):
    pred = read_pred(RAW / "wms_analysis" / "wms164" / f"{sample}.pred")
    true_spec = truth_by_mock[mock]
    true_set = set()
    t = load_truth(mock)
    for atcc in t:
        num = re.sub(r"^ATCC[_ ]*", "", atcc).strip()
        pat = rf"(?:^|_){re.escape(num)}(_|$)"
        for g in (g for gs in MEMBERS.values() for g in gs):
            if "ATCC_" in g and re.search(pat, g):
                true_set.add(g)
                break
    true_tokens = {atcc_token(g) for g in true_set if atcc_token(g)}

    fps = []
    for g, ab in pred.items():
        if ab <= 0:
            continue
        tok = atcc_token(g)
        is_tp = (g in true_set) or (tok and tok in true_tokens)
        if is_tp:
            continue
        cls = classify_fp(g, true_spec, truth_by_mock)
        fps.append({"genome": g, "species": species_of(g), "abundance": ab, "class": cls})
    return sorted(fps, key=lambda x: -x["abundance"])

def main():
    if not (RAW / "MSA_combined164_all_flat.members.tsv").exists():
        print(f"ERROR: HPC mirror not found at {RAW}")
        return

    truth_by_mock = all_truth_species_per_mock()
    rows = []
    profiles = {}
    summary = {}
    SAMPLES = {
        "MSA1002": ["WMS_MSA1002_0_100ng_1", "WMS_MSA1002_0_100ng_2",
                      "WMS_MSA1002_90_100ng_1", "WMS_MSA1002_95_100ng_1", "WMS_MSA1002_99_100ng_1"],
        "MSA1003": [f"WMS_MSA1003_0_100ng_{r}" for r in (1, 2, 3)],
        "MSA1005": [f"WMS_MSA1005_0_100ng_{r}" for r in (1, 2, 3)],
        "MSA1007": [f"WMS_MSA1007_0_100ng_{r}" for r in (1, 2, 3)],
    }
    for mock, samples in SAMPLES.items():
        for sample in samples:
            fps = analyze_sample(sample, mock, truth_by_mock)
            if not fps:
                continue
            profiles[sample] = fps
            for fp in fps:
                rows.append({
                    "sample": sample,
                    "mock": mock,
                    "class": fp["class"],
                    "species": fp["species"],
                    "genome": fp["genome"],
                    "abundance": fp["abundance"],
                })
            summary[sample] = {"mock": mock}
            for cls in ["other_mock_atcc", "same_species_decoy", "other_mock_gcf", "other"]:
                n = sum(1 for fp in fps if fp["class"] == cls)
                ab = sum(fp["abundance"] for fp in fps if fp["class"] == cls)
                summary[sample][cls] = {"count": n, "total_abundance": ab}

    RES.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "mock", "class", "species", "genome", "abundance"], delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(OUT_JSON, "w") as fh:
        json.dump({"profiles": profiles, "summary": summary}, fh, indent=2)

    # Print summary table
    print("\n=== FP source summary per sample ===")
    print(f"{'sample':<25} {'total_FP':>8} {'other_mock_atcc':>18} {'same_species_decoy':>20} {'other':>8}")
    for sample, s in sorted(summary.items()):
        total = sum(v["count"] for v in s.values() if isinstance(v, dict))
        print(f"{sample:<25} {total:>8} {s['other_mock_atcc']['count']:>18} {s['same_species_decoy']['count']:>20} {s['other']['count']:>8}")

    print(f"\nwrote {OUT_TSV}")
    print(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
