#!/usr/bin/env python3
"""Compare two fig6_fig12_metrics.tsv files, focusing on Strain2bScan rows."""
import csv, sys
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent.parent

FLOAT_COLS = ["precision", "recall", "f1", "aupr",
              "precision_at_1e_4", "recall_at_1e_4", "f1_at_1e_4",
              "bray_curtis", "l2"]
INT_COLS = ["TP", "FP", "FN"]

def load(path):
    rows = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        key = (r["kind"], r["mock"], r["sample"], r["tool"], r["variant"])
        rows.append((key, r))
    return dict(rows)

def main(old_path, new_path):
    old = load(old_path)
    new = load(new_path)
    keys = set(old) | set(new)
    changed = []
    unchanged = 0
    missing_old = []
    missing_new = []
    for k in sorted(keys):
        if k not in old:
            missing_old.append(k)
            continue
        if k not in new:
            missing_new.append(k)
            continue
        o, n = old[k], new[k]
        diffs = []
        for c in INT_COLS:
            if int(o[c]) != int(n[c]):
                diffs.append(f"{c}: {o[c]} -> {n[c]}")
        for c in FLOAT_COLS:
            ov = float(o[c]) if o[c] not in ("", "nan") else float("nan")
            nv = float(n[c]) if n[c] not in ("", "nan") else float("nan")
            if abs(ov - nv) > 1e-4:
                diffs.append(f"{c}: {ov:.6f} -> {nv:.6f}")
        if diffs:
            changed.append((k, diffs))
        else:
            unchanged += 1
    print(f"Rows: {len(keys)}  unchanged: {unchanged}  changed: {len(changed)}  missing_old: {len(missing_old)}  missing_new: {len(missing_new)}")
    if missing_old:
        print("\nMissing in old:")
        for k in missing_old:
            print(" ", k)
    if missing_new:
        print("\nMissing in new:")
        for k in missing_new:
            print(" ", k)
    if changed:
        print("\nChanged rows (delta > 1e-4):")
        for k, diffs in changed:
            print(f"\n{k}")
            for d in diffs:
                print(f"  {d}")
    return 1 if changed or missing_old or missing_new else 0

if __name__ == "__main__":
    old = Path(sys.argv[1]) if len(sys.argv) > 1 else PAPER / "data" / "fig6_fig12_metrics.tsv.bak_20261003_104241"
    new = Path(sys.argv[2]) if len(sys.argv) > 2 else PAPER / "data" / "fig6_fig12_metrics.tsv"
    sys.exit(main(old, new))
