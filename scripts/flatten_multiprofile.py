#!/usr/bin/env python3
"""Convert Strain2bScan multi-profile .pred to a flat profile-style .pred.

Renames clusters to species__cluster, keeps coverage/support/depth/n_markers,
and uses global_abundance as the abundance column. Rows below --threshold or
rows that are detected-not-resolvable are dropped.
"""
import argparse, sys
from pathlib import Path


def flatten(in_path, out_path, threshold):
    with open(in_path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
    if "#species" not in hdr:
        print(f"warning: {in_path} is not a multi-profile pred; copying as-is", file=sys.stderr)
        Path(in_path).replace(out_path)
        return
    species_i = hdr.index("#species")
    cluster_i = species_i + 1
    abundance_i = hdr.index("global_abundance")
    coverage_i = hdr.index("coverage")
    support_i = hdr.index("support")
    depth_i = hdr.index("depth")
    n_markers_i = hdr.index("n_markers")

    kept = 0
    with open(in_path) as f, open(out_path, "w") as w:
        w.write("#cluster\tabundance\tcoverage\tsupport\tdepth\tn_markers\n")
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) <= cluster_i:
                continue
            species, cluster = p[species_i].strip(), p[cluster_i].strip()
            if cluster.startswith("[") or cluster.startswith("no cluster"):
                continue
            try:
                ga = float(p[abundance_i])
            except ValueError:
                continue
            if ga < threshold:
                continue
            kept += 1
            w.write(f"{species}__{cluster}\t{ga}\t{p[coverage_i]}\t{p[support_i]}\t{p[depth_i]}\t{p[n_markers_i]}\n")
    print(f"{in_path} -> {out_path}: kept {kept} clusters above global_abundance {threshold}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", "-i", dest="in_pred", required=True)
    ap.add_argument("--out", "-o", required=True)
    ap.add_argument("--threshold", "-t", type=float, default=0.001)
    args = ap.parse_args()
    flatten(args.in_pred, args.out, args.threshold)


if __name__ == "__main__":
    main()
