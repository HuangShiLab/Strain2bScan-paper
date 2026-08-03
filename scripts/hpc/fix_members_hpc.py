#!/usr/bin/env python3
"""Normalize per-species cluster names in *_flat.members.tsv to global bare IDs."""
from collections import defaultdict
from pathlib import Path

ROOT = Path("/lustre1/g/aos_shihuang/LU/Strain2bScan-raw-data")
FILES = [
    ROOT / "MSA_combined164_all_flat.members.tsv",
    ROOT / "MSA1002_combined_all_flat.members.tsv",
]

def fix(path):
    rows = []
    prefixed = False
    with open(path) as fh:
        header = fh.readline()
        for ln in fh:
            if not ln.strip():
                continue
            g, c = ln.rstrip("\n").split("\t")[:2]
            rows.append((g, c))
            if "__C" in c:
                prefixed = True
    if not prefixed:
        print(f"{path}: already bare cluster IDs, no change")
        return
    # group by species, then by cluster number
    by_sp = defaultdict(lambda: defaultdict(list))
    for g, c in rows:
        sp, cid = c.rsplit("__", 1)
        by_sp[sp][int(cid[1:])].append(g)
    new_rows = []
    idx = 0
    for sp in sorted(by_sp):
        for cnum in sorted(by_sp[sp]):
            new_c = f"C{idx}"
            for g in by_sp[sp][cnum]:
                new_rows.append((g, new_c))
            idx += 1
    out = path.with_suffix(".members.tsv.new")
    with open(out, "w") as fh:
        fh.write(header)
        for g, c in new_rows:
            fh.write(f"{g}\t{c}\n")
    bak = path.with_name(path.name + ".bak")
    path.rename(bak)
    out.rename(path)
    print(f"{path}: normalized {len(new_rows)} genomes -> {idx} clusters")

if __name__ == "__main__":
    for f in FILES:
        fix(f)
