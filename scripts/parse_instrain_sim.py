#!/usr/bin/env python3
"""Parse inStrain quick_profile outputs for the sim benchmark and compute per-sample
precision/recall/F1 relative to ground truth.

Inputs:
  work/mock_retest/Strain2bScan-raw-data/sim_benchmark/results/instrain/{single,multi}/
  figure_raw_data/sim_single_species/<Species>/truth/<sample>.truth.tsv
  figure_raw_data/sim_multi_species/<depth>/truth/sample01.truth.tsv
  results/genome_to_species.tsv

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


def load(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def read_truth_single(path):
    """Return cluster -> relative_abundance."""
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
    """Return species -> summed relative_abundance.

    Truth files use underscores in species names (e.g. Escherichia_coli);
    genome_to_species.tsv uses canonical spaces (Escherichia coli).
    Normalize to the space-separated form for comparison.
    """
    out = defaultdict(float)
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 4:
                species = p[0].replace("_", " ")
                out[species] += float(p[3])
    return dict(out)


def read_genome_coverage(path, min_cov=0.1, min_breadth=0.05):
    """Read inStrain genomeCoverage.csv; return set of detected genome IDs."""
    det = set()
    if not path.exists():
        return det
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            genome = row.get("genome")
            if not genome:
                continue
            try:
                cov = float(row.get("coverage", 0))
                breadth = float(row.get("breadth", 0))
            except ValueError:
                continue
            if cov >= min_cov and breadth >= min_breadth:
                det.add(genome)
    return det


def load_genome_to_species():
    """Map genome ID to species name from the multi-reference STB/FASTA."""
    path = RES / "genome_to_species.tsv"
    out = {}
    if path.exists():
        for row in load(path):
            out[row["genome"]] = row["species"]
    return out


def strain_level_metrics(n_detected, n_truth):
    """Single-species evaluation: detected genomes vs true strains.

    All detected genomes belong to the same species as the truth strains.
    We credit detection up to the true number of strains and count excess
    detections as false positives / over-splitting.
    """
    tp = min(n_detected, n_truth)
    fp = max(0, n_detected - n_truth)
    fn = max(0, n_truth - n_detected)
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return {"TP": tp, "FP": fp, "FN": fn, "precision": prec, "recall": rec, "F1": f1}


def species_level_metrics(detected_species, truth_species):
    """Multi-species evaluation: any detected genome -> species present."""
    tp = len(truth_species & detected_species)
    fp = len(detected_species - truth_species)
    fn = len(truth_species - detected_species)
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return {"TP": tp, "FP": fp, "FN": fn, "precision": prec, "recall": rec, "F1": f1}


def main():
    genome_to_species = load_genome_to_species()
    rows = []

    # single species
    for sp_dir in sorted(INSTRAIN_ROOT.glob("single/*")):
        sample = sp_dir.name
        # samples are named like "Escherichia_coli__diff_k2_rep1_d5"
        species = sample.split("__")[0] if "__" in sample else None
        truth_path = TRUTH_SINGLE / species / "truth" / f"{sample}.truth.tsv" if species else None
        if not truth_path or not truth_path.exists():
            print(f"[skip single {sample}] truth missing: {truth_path}")
            continue
        truth = read_truth_single(truth_path)
        detected = read_genome_coverage(sp_dir / "genomeCoverage.csv")
        m = strain_level_metrics(len(detected), len(truth))
        rows.append({"kind": "single", "sample": sample, "tool": "inStrain",
                     "n_truth": len(truth), "n_pred": len(detected), **m})

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
        detected_genomes = read_genome_coverage(depth_dir / "genomeCoverage.csv")
        detected_species = {genome_to_species.get(g, "unknown") for g in detected_genomes}
        detected_species.discard("unknown")
        m = species_level_metrics(detected_species, set(truth))
        rows.append({"kind": "multi", "sample": sample, "tool": "inStrain",
                     "n_truth": len(truth), "n_pred": len(detected_species), **m})

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
