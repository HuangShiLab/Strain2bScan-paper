#!/usr/bin/env python3
"""Parse multi-profile stdout to extract per-species depth for detected-but-not-resolved species."""
import sys, re
from collections import defaultdict

def parse_stream(fh):
    sample = None
    rows = defaultdict(dict)  # sample -> species -> depth
    for ln in fh:
        ln = ln.rstrip("\n")
        m = re.search(r"RUN\s+(\S+)", ln)
        if m:
            sample = m.group(1)
            continue
        if sample is None:
            continue
        if "[detected, not strain-resolvable]" in ln or "[strain-resolved, no cluster above threshold]" in ln:
            parts = ln.split("\t")
            if not parts:
                continue
            sp = parts[0].strip().replace(".db", "")
            if not sp:
                continue
            dm = re.search(r"depth\s+([0-9.eE+-]+)x", ln)
            depth = float(dm.group(1)) if dm else 0.0
            rows[sample][sp] = depth
    return rows

def main():
    if len(sys.argv) != 3:
        print("usage: parse_detected_depth.py <stdout.log> <out.tsv>", file=sys.stderr)
        sys.exit(1)
    log_path, out_path = sys.argv[1:3]
    with open(log_path) as f:
        rows = parse_stream(f)
    if not rows:
        print("warning: no detected-not-resolved lines found", file=sys.stderr)
    samples = sorted(rows.keys())
    species = sorted({sp for d in rows.values() for sp in d})
    with open(out_path, "w") as f:
        f.write("sample\t" + "\t".join(species) + "\n")
        for s in samples:
            vals = [str(rows[s].get(sp, 0.0)) for sp in species]
            f.write(s + "\t" + "\t".join(vals) + "\n")
    print(f"wrote {out_path}: {len(samples)} samples x {len(species)} species")

if __name__ == "__main__":
    main()
