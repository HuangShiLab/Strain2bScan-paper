#!/usr/bin/env python3
"""Sweep a global-abundance (or abundance) threshold across Strain2bScan
multi-profile .pred files and report species-level precision/recall vs truth.

Usage:
    python3 scripts/sweep_global_abundance.py \
        --preds scan_msa1005/default.pred scan_msa1005/no_floor.pred \
        --truth Ground_truth/MSA1005_ground_truth.txt \
        --genome-dir MSA_all164 \
        --out sweep_msa1005.tsv
"""
import argparse, csv, sys
from pathlib import Path


def truth_species(truth_file, genome_dir):
    """Return set of true species names from a ground-truth ATCC list."""
    atccs = []
    with open(truth_file) as f:
        for ln in f:
            first = ln.lstrip()
            if not first or first.startswith("base_counts") or first.startswith("genome"):
                continue
            atccs.append(ln.split("\t")[0].strip())
    species = set()
    gdir = Path(genome_dir)
    for atcc in atccs:
        matches = list(gdir.glob(f"*{atcc}*.fna"))
        if matches:
            species.add(matches[0].name.split("__")[0])
    return species


def predicted_species_at_threshold(pred_file, key, threshold):
    """Return species with a strain-resolved cluster whose key >= threshold."""
    species = set()
    with open(pred_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            key_idx = header.index(key)
        except ValueError:
            raise ValueError(f"{pred_file}: header {header!r} has no '{key}' column")
        species_idx = header.index("#species") if "#species" in header else None
        for ln in f:
            if not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            if species_idx is not None:
                sp = parts[species_idx]
                cl = parts[species_idx + 1]
            else:
                sp = parts[0]
                cl = parts[1]
            if cl.startswith("[") or cl.startswith("no cluster"):
                continue
            try:
                val = float(parts[key_idx])
            except ValueError:
                continue
            if val >= threshold:
                species.add(sp.strip())
    return species


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preds", nargs="+", required=True, help=".pred files")
    ap.add_argument("--truth", required=True, help="ground truth tsv")
    ap.add_argument("--genome-dir", required=True, help="dir with .fna files for ATCC->species mapping")
    ap.add_argument("--key", default="global_abundance", choices=["global_abundance", "abundance", "sample_fraction"])
    ap.add_argument("--thresholds", nargs="+", type=float,
                    default=[0.0, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    truth = truth_species(args.truth, args.genome_dir)
    print(f"Truth species ({len(truth)}): {sorted(truth)}", file=sys.stderr)

    rows = []
    for pred in args.preds:
        pred = Path(pred)
        for t in sorted(args.thresholds):
            pred_sp = predicted_species_at_threshold(pred, args.key, t)
            tp = len(pred_sp & truth)
            fp = len(pred_sp - truth)
            fn = len(truth - pred_sp)
            prec = tp / (tp + fp) if (tp + fp) else 0.0
            rec = tp / (tp + fn) if (tp + fn) else 0.0
            rows.append({
                "pred": pred.stem,
                "threshold": t,
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "precision": prec,
                "recall": rec,
                "n_pred": len(pred_sp),
            })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh,
                           fieldnames=["pred", "threshold", "TP", "FP", "FN", "precision", "recall", "n_pred"],
                           delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
