#!/usr/bin/env python3
"""Compare newly generated .pred files to the original Strain2bScan-raw-data predictions."""
import sys, difflib
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent.parent.parent
OLD = PAPER / "work" / "mock_retest" / "Strain2bScan-raw-data"
NEW = Path(__file__).resolve().parent

MAPPING = [
    (NEW / "wms_analysis", OLD / "wms_analysis"),
    (NEW / "wms_analysis_tracegap", OLD / "wms_analysis_tracegap"),
    (NEW / "Strain2bScan-port-results" / "mock", OLD.parent / "Strain2bScan-port-results" / "mock"),
]

def main():
    total = same = diff = missing = 0
    for new_root, old_root in MAPPING:
        if not new_root.exists():
            continue
        for new_pred in new_root.rglob("*.pred"):
            rel = new_pred.relative_to(new_root)
            old_pred = old_root / rel
            total += 1
            if not old_pred.exists():
                missing += 1
                print(f"[missing old] {rel}")
                continue
            new_text = new_pred.read_text()
            old_text = old_pred.read_text()
            if new_text == old_text:
                same += 1
            else:
                diff += 1
                print(f"[diff] {rel}")
                for line in difflib.unified_diff(
                    old_text.splitlines(keepends=True),
                    new_text.splitlines(keepends=True),
                    fromfile=str(old_pred),
                    tofile=str(new_pred),
                    n=2,
                ):
                    sys.stdout.write(line)
    print(f"\nTotal: {total}  Same: {same}  Diff: {diff}  Missing old: {missing}")
    return 0 if diff == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
