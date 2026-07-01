#!/usr/bin/env python3
"""Simulate strain-mixture mock metagenomes from a real genome panel (C. acnes-style design):
each sample mixes 2-5 random strains at log-normal depth (>=1x), error-free 150 bp reads.
Usage: python sim_strain_mock.py <genomes_dir> <out_dir> [n_samples] [seed]
Writes <out_dir>/reads/sample{i}.fq and <out_dir>/truth/sample{i}.truth.tsv (relative abundance).
"""
import random, glob, os, sys
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
def read_seq(p): return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()
def sim_reads(seq, depth, fh, tag, rng):
    L = len(seq); n = int(depth * L / 150)
    for i in range(n):
        p = rng.randint(0, max(0, L - 150)); r = seq[p:p + 150]
        if rng.random() < 0.5: r = "".join(COMP.get(c, "N") for c in reversed(r))
        fh.write(f"@{tag}_{i}\n{r}\n+\n{'I'*150}\n")
def main(gdir, outdir, nsamp=5, seed=42):
    rng = random.Random(seed)
    genomes = sorted(glob.glob(f"{gdir}/*.fna"))
    os.makedirs(f"{outdir}/reads", exist_ok=True); os.makedirs(f"{outdir}/truth", exist_ok=True)
    for s in range(1, nsamp + 1):
        k = rng.randint(2, 5); chosen = rng.sample(genomes, k)
        depths = [max(1.0, round(rng.lognormvariate(1.0, 0.6), 2)) for _ in chosen]
        with open(f"{outdir}/reads/sample{s}.fq", "w") as fq:
            for g, d in zip(chosen, depths):
                sim_reads(read_seq(g), d, fq, os.path.basename(g)[:-4], rng)
        tot = sum(depths)
        with open(f"{outdir}/truth/sample{s}.truth.tsv", "w") as t:
            t.write("#genome\tabundance\n")
            for g, d in zip(chosen, depths):
                t.write(f"{os.path.basename(g)[:-4]}\t{d/tot:.4f}\n")
        print(f"  sample{s}: {k} strains, depths {depths}", flush=True)
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 5,
         int(sys.argv[4]) if len(sys.argv) > 4 else 42)
