#!/usr/bin/env python3
"""Validate the simulated degraded genomes with CheckM2 (LINUX; CheckM2 has no osx-arm64 build).

Regenerates the degraded C. acnes truth genomes at each quality level (same seeds as
run_refqual.py), runs `checkm2 predict`, and tabulates DESIGNED vs CheckM2-MEASURED
completeness/contamination, so the supplementary figure's x-axis can be reported as CheckM2
scores. Run from the workdir.

Prereqs:
  conda create -n checkm2 -c bioconda -c conda-forge checkm2
  conda activate checkm2 && checkm2 database --download   # ~3 GB diamond DB (once)
"""
import os, glob, subprocess, csv, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from degrade import degrade, read_seq

GEN, TRUTH, WORK = "acnes/genomes", "acnes/truth", "refqual"
LEVELS = [(1.00, 0.00, 1), (0.95, 0.01, 20), (0.90, 0.02, 50),
          (0.80, 0.05, 100), (0.70, 0.08, 200), (0.50, 0.10, 400)]

truth_acc = set()
for i in range(1, 6):
    for l in open(f"{TRUTH}/sample{i}.truth.tsv"):
        if not l.startswith("#"):
            truth_acc.add(l.split("\t")[0])
contaminant = read_seq(sorted(glob.glob("multispecies/genomes/Escherichia_coli/*.fna"))[0])


def write_fasta(frags, path):
    with open(path, "w") as f:
        for i, fr in enumerate(frags):
            f.write(f">ctg{i}\n")
            for j in range(0, len(fr), 70):
                f.write(fr[j:j + 70] + "\n")


out = open(f"{WORK}/refqual_checkm2.tsv", "w")
out.write("completeness_designed\tcontamination_designed\t"
          "completeness_checkm2_mean\tcontamination_checkm2_mean\n")
for comp, contam, nc in LEVELS:
    tag = int(comp * 100)
    indir = f"{WORK}/checkm2_in_{tag}"
    os.makedirs(indir, exist_ok=True)
    for f in glob.glob(f"{indir}/*.fna"):
        os.remove(f)
    for acc in truth_acc:
        seq = read_seq(f"{GEN}/{acc}.fna")
        if comp >= 1.0:
            import shutil
            shutil.copy(f"{GEN}/{acc}.fna", f"{indir}/{acc}.fna")
        else:
            write_fasta(degrade(seq, comp, contam, nc, contaminant, sum(map(ord, acc)) + tag),
                        f"{indir}/{acc}.fna")
    od = f"{WORK}/checkm2_out_{tag}"
    subprocess.run(["checkm2", "predict", "--input", indir, "--output-directory", od,
                    "-x", "fna", "--force"], check=True)
    comps, conts = [], []
    for r in csv.DictReader(open(f"{od}/quality_report.tsv"), delimiter="\t"):
        comps.append(float(r["Completeness"]))
        conts.append(float(r["Contamination"]))
    out.write(f"{comp:.2f}\t{contam:.2f}\t{statistics.mean(comps):.1f}\t{statistics.mean(conts):.1f}\n")
    print(f"designed comp={comp:.2f} contam={contam:.2f} -> "
          f"CheckM2 comp={statistics.mean(comps):.1f}% contam={statistics.mean(conts):.1f}%", flush=True)
out.close()
print("done -> refqual/refqual_checkm2.tsv")
