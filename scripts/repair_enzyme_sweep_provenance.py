#!/usr/bin/env python3
"""Add row-specific prediction paths and hashes to the completed enzyme-sweep table."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "results/enzyme_sweep_225/persample.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


with TABLE.open(newline="") as handle:
    rows = list(csv.reader(handle, delimiter="\t"))

header = rows[0]
if "prediction_sha256" not in header:
    header.insert(header.index("prediction") + 1, "prediction_sha256")
prediction_i = header.index("prediction")
hash_i = header.index("prediction_sha256")
truth_i = header.index("truth")
database_i = header.index("database")

for row in rows[1:]:
    row[:] = [value.strip() for value in row]
    if len(row) == len(header) - 1:
        row.insert(hash_i, "")
    relative = Path("work/enzyme_sweep_225/predictions") / f"enzymes_{row[header.index('enzyme_count')]}" / f"{row[header.index('sample')]}.pred"
    path = ROOT / relative
    if not path.exists():
        raise SystemExit(f"missing archived prediction: {relative}")
    row[prediction_i] = relative.as_posix()
    row[hash_i] = sha256(path)
    for column in (truth_i, database_i):
        path = Path(row[column])
        row[column] = path.relative_to(ROOT).as_posix() if path.is_absolute() else row[column]

with TABLE.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerows(rows)

print(f"repaired {len(rows) - 1} rows in {TABLE}")
