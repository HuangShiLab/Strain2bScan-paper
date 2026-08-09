#!/usr/bin/env python3
"""Parse inStrain quick_profile outputs for the sim benchmark and compute per-sample
precision/recall/F1 relative to ground truth.

Inputs:
  work/mock_retest/Strain2bScan-raw-data/sim_benchmark/results/instrain/{single,multi}/
  figure_raw_data/sim_single_species/<Species>/truth/<sample>.truth.tsv
  figure_raw_data/sim_multi_species/<depth>/truth/sample01.truth.tsv

Outputs:
  results/instrain_sim.tsv
"""
import csv
import re
from pathlib import Path
from collections import defaultdict

PAPER = Path(__file__).resolve().parent.parent
INSTRAIN_ROOT = PAPER / "work" / "mock_retest" / "Strain2bScan-raw-data" / "sim_benchmark" / "results" / "instrain"
TRUTH_SINGLE = PAPER / "figure_raw_data" / "sim_single_species"
TRUTH_MULTI = PAPER / "figure_raw_data" / "sim_multi_species"
RES = PAPER / "results"


def read_truth_single(path):
    out = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 3:
                out[p[1]] = float(p[2])  # cluster -> relative_abundance
    return out


def read_truth_multi(path):
    out = defaultdict(float)
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 4:
                out[p[0]] += float(p[3])  # species -> summed relative_abundance
    return dict(out)


def read_coverm(path):
    """Read inStrain quick_profile coverm_raw.tsv; return dict genome->(coverage, breadth, rel_ab)."""
    if not path.exists():
        return {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            return {}
        out = {}
        for row in reader:
            genome = row.get("Genome") or row.get("genome") or row.get("#Genome")
            if not genome:
                continue
            cov = _float(row, ["Mean", "coverage", "Coverage", "mean"])
            breadth = _float(row, ["Covered Fraction", "breadth", "Breadth", "covered_fraction"])
            rel = _float(row, ["Relative Abundance", "relative_abundance", "rel_abundance", "rel_ab"])
            out[genome] = {"coverage": cov or 0.0, "breadth": breadth or 0.0, "rel_ab": rel or 0.0}
        return out


def _float(row, keys):
    for k in keys:
        if k in row and row[k]:
            try:
                return float(row[k])
            except ValueError:
                pass
    return None


def metrics(pred, truth, present_cov=0.1, present_breadth=0.05):
    pos = set(truth)
    # if inStrain reported relative abundance, use that for detection; otherwise coverage+breadth
    if any(v.get("rel_ab", 0) > 0 for v in pred.values()):
        det = {g for g, v in pred.items() if v.get("rel_ab", 0) > 0}
    else:
        det = {g for g, v in pred.items() if v.get("coverage", 0) >= present_cov and v.get("breadth", 0) >= present_breadth}
    tp = len(pos & det)
    fp = len(det - pos)
    fn = len(pos - det)
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return {"TP": tp, "FP": fp, "FN": fn, "precision": prec, "recall": rec, "F1": f1}


def main():
    rows = []
    # single species
    for sp_dir in sorted(INSTRAIN_ROOT.glob("single/*")):
        sample = sp_dir.name
        m = re.match(r"^([^_]+)_(.+)$", sample)
        species = m.group(1).replace("_", " ") if m else None
        truth_path = TRUTH_SINGLE / species / "truth" / f"{sample}.truth.tsv" if species else None
        if not truth_path or not truth_path.exists():
            print(f"[skip single {sample}] truth missing: {truth_path}")
            continue
        truth = read_truth_single(truth_path)
        coverm = sp_dir / "coverm_raw.tsv"
        pred = read_coverm(coverm)
        m = metrics(pred, truth)
        rows.append({"kind": "single", "sample": sample, "tool": "inStrain",
                     "n_truth": len(truth), "n_pred": len(pred), **m})

    # multi species
    for depth_dir in sorted(INSTRAIN_ROOT.glob("multi/*")):
        sample = depth_dir.name  # depth_low_sample01
        m = re.match(r"^(depth_[a-z]+)_sample01$", sample)
        depth = m.group(1) if m else None
        truth_path = TRUTH_MULTI / depth / "truth" / "sample01.truth.tsv" if depth else None
        if not truth_path or not truth_path.exists():
            print(f"[skip multi {sample}] truth missing: {truth_path}")
            continue
        truth = read_truth_multi(truth_path)
        coverm = depth_dir / "coverm_raw.tsv"
        pred = read_coverm(coverm)
        m = metrics(pred, truth)
        rows.append({"kind": "multi", "sample": sample, "tool": "inStrain",
                     "n_truth": len(truth), "n_pred": len(pred), **m})

    RES.mkdir(parents=True, exist_ok=True)
    cols = ["kind", "sample", "tool", "n_truth", "n_pred", "TP", "FP", "FN", "precision", "recall", "F1"]
    with open(RES / "instrain_sim.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(r[k], 4) if isinstance(r[k], float) else r[k]) for k in cols})
    print(f"wrote {RES / 'instrain_sim.tsv'} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
