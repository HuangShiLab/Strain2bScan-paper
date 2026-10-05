#!/usr/bin/env python3
"""Refresh SHA-256 values in figures/numbered/MANIFEST.tsv after regeneration."""
import csv, hashlib
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
path = ROOT / "figures/numbered/MANIFEST.tsv"
with path.open() as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
    fields = rows[0].keys()
for row in rows:
    for ext in ("png", "pdf"):
        row[f"{ext}_sha256"] = hashlib.sha256((ROOT / "figures/numbered" / row[ext]).read_bytes()).hexdigest()
with path.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
print(f"refreshed {len(rows)} checksums in {path}")
