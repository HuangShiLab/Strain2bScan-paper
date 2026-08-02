#!/usr/bin/env python3
"""Compare Strain2bScan-port simulation predictions against ground truth.

Inputs:
  work/mock_retest/Strain2bScan-port-results/sim/{single,multi}/{mode}/*.pred
  figure_raw_data/sim_single_species/<Species>/truth/<sample>.truth.tsv
  figure_raw_data/sim_multi_species/<depth>/truth/<sample>.truth.tsv

Outputs:
  results/sim_port_comparison.tsv   per-sample precision/recall
  stdout summary tables
"""

import csv
import glob
import os
import re
from collections import defaultdict
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
PRED_ROOT = PAPER / "work" / "mock_retest" / "Strain2bScan-port-results" / "sim"
TRUTH_SINGLE = PAPER / "figure_raw_data" / "sim_single_species"
TRUTH_MULTI = PAPER / "figure_raw_data" / "sim_multi_species"
OUT_DIR = PAPER / "results"


def read_truth_single(path):
    """Return set of true cluster IDs."""
    clusters = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and float(parts[2]) > 0:
                clusters.add(parts[1])
    return clusters


def read_pred_single(path):
    """Return set of predicted cluster IDs."""
    clusters = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts:
                clusters.add(parts[0])
    return clusters


def read_truth_multi(path):
    """Return set of true species names."""
    species = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 1 and float(parts[3]) > 0:
                species.add(parts[0])
    return species


def read_pred_multi(path):
    """Return set of predicted species names."""
    species = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts:
                species.add(parts[0])
    return species


def metrics(tp, fp, fn):
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return precision, recall, f1


def compare_single(mode):
    rows = []
    pred_dir = PRED_ROOT / "single" / mode
    for pred_path in sorted(pred_dir.glob("*.pred")):
        base = pred_path.stem
        species = base.split("__")[0]
        truth_path = TRUTH_SINGLE / species / "truth" / f"{base}.truth.tsv"
        if not truth_path.exists():
            print(f"[skip single {mode}] truth missing: {truth_path}")
            continue
        true_clusters = read_truth_single(truth_path)
        pred_clusters = read_pred_single(pred_path)
        tp = len(true_clusters & pred_clusters)
        fp = len(pred_clusters - true_clusters)
        fn = len(true_clusters - pred_clusters)
        p, r, f1 = metrics(tp, fp, fn)
        rows.append({
            "mode": mode,
            "dataset": "single",
            "sample": base,
            "n_true": len(true_clusters),
            "n_pred": len(pred_clusters),
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "precision": p,
            "recall": r,
            "F1": f1,
        })
    return rows


def compare_multi(mode):
    rows = []
    pred_dir = PRED_ROOT / "multi" / mode
    for pred_path in sorted(pred_dir.glob("*.pred")):
        base = pred_path.stem
        m = re.match(r"^(depth_[a-z]+)_(sample\d+)$", base)
        if not m:
            print(f"[skip multi {mode}] unparseable pred name: {base}")
            continue
        depth, sample = m.groups()
        truth_path = TRUTH_MULTI / depth / "truth" / f"{sample}.truth.tsv"
        if not truth_path.exists():
            print(f"[skip multi {mode}] truth missing: {truth_path}")
            continue
        true_species = read_truth_multi(truth_path)
        pred_species = read_pred_multi(pred_path)
        tp = len(true_species & pred_species)
        fp = len(pred_species - true_species)
        fn = len(true_species - pred_species)
        p, r, f1 = metrics(tp, fp, fn)
        rows.append({
            "mode": mode,
            "dataset": "multi",
            "sample": base,
            "n_true": len(true_species),
            "n_pred": len(pred_species),
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "precision": p,
            "recall": r,
            "F1": f1,
        })
    return rows


def aggregate(rows):
    by_mode_dataset = defaultdict(lambda: {"TP": 0, "FP": 0, "FN": 0})
    for r in rows:
        key = (r["mode"], r["dataset"])
        by_mode_dataset[key]["TP"] += r["TP"]
        by_mode_dataset[key]["FP"] += r["FP"]
        by_mode_dataset[key]["FN"] += r["FN"]
    out = []
    for (mode, dataset), agg in sorted(by_mode_dataset.items()):
        p, r, f1 = metrics(agg["TP"], agg["FP"], agg["FN"])
        out.append({
            "mode": mode,
            "dataset": dataset,
            "TP": agg["TP"],
            "FP": agg["FP"],
            "FN": agg["FN"],
            "precision": p,
            "recall": r,
            "F1": f1,
        })
    return out


def discover_modes(dataset):
    """Return sorted list of mode subdirectory names under sim/{dataset}."""
    d = PRED_ROOT / dataset
    if not d.exists():
        return []
    return sorted(p.name for p in d.iterdir() if p.is_dir())


def main():
    rows = []
    for mode in discover_modes("single"):
        rows.extend(compare_single(mode))
    for mode in discover_modes("multi"):
        rows.extend(compare_multi(mode))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "sim_port_comparison.tsv"
    fieldnames = ["mode", "dataset", "sample", "n_true", "n_pred",
                  "TP", "FP", "FN", "precision", "recall", "F1"]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out_path}")

    print("\n=== Per-sample results ===")
    print(f"{'mode':<10} {'dataset':<8} {'sample':<45} {'P':>6} {'R':>6} {'F1':>6}  TP/FP/FN")
    for r in rows:
        print(f"{r['mode']:<10} {r['dataset']:<8} {r['sample']:<45} "
              f"{r['precision']:>6.3f} {r['recall']:>6.3f} {r['F1']:>6.3f}  "
              f"{r['TP']}/{r['FP']}/{r['FN']}")

    print("\n=== Aggregated (micro-averaged) ===")
    agg = aggregate(rows)
    print(f"{'mode':<10} {'dataset':<8} {'P':>6} {'R':>6} {'F1':>6}  TP/FP/FN")
    for a in agg:
        print(f"{a['mode']:<10} {a['dataset']:<8} {a['precision']:>6.3f} {a['recall']:>6.3f} {a['F1']:>6.3f}  "
              f"{a['TP']}/{a['FP']}/{a['FN']}")


if __name__ == "__main__":
    main()
