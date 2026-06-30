#!/usr/bin/env python3
"""StrainScan arm of the reference-quality degradation benchmark — LINUX ONLY.

Mirrors run_refqual.py but builds a StrainScan k-mer DB at each quality level and profiles the
SAME fixed (high-quality) samples. StrainScan_build needs `dashing` (and sibeliaz), which have
no osx-arm64 builds, so this arm must run on Linux (e.g. the HPC), not the dev Mac.

Prereqs (Linux):
  conda env create -f env/strainscan-linux.yml      # or: conda install -c bioconda strainscan
  export STRAINSCAN=/path/to/StrainScan              # repo with StrainScan_build.py / StrainScan.py
  export STRAINSCAN_PY=/path/to/strainscan-env/bin/python
  run from a $WORKDIR containing acnes/{genomes,reads,truth}/ and
  multispecies/genomes/Escherichia_coli/  (regenerate with the other scripts first)

Design: vary REFERENCE DB quality, keep samples fixed (the clean test of "reference quality ->
profiling accuracy"). The 14 C. acnes truth strains are degraded (completeness down, contigs
up, contamination up); the 50 background genomes and the 5 samples are unchanged.
"""
import os, glob, subprocess, re, shutil, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from degrade import degrade, read_seq

SS = os.environ.get("STRAINSCAN", "StrainScan")
PY = os.environ.get("STRAINSCAN_PY", "python")
GEN, TRUTH, READS, WORK = "acnes/genomes", "acnes/truth", "acnes/reads", "refqual_ss"
os.makedirs(WORK, exist_ok=True)
LEVELS = [(1.00, 0.00, 1), (0.95, 0.01, 20), (0.90, 0.02, 50),
          (0.80, 0.05, 100), (0.70, 0.08, 200), (0.50, 0.10, 400)]


def strip_ver(acc):  # GCF_009737125.1 -> GCF_009737125 (StrainScan reports version-less)
    return re.sub(r"\.\d+$", "", acc)


# truth strains + per-sample truth (accession -> abundance)
truth_acc, sample_truth = set(), {}
for i in range(1, 6):
    rows = []
    for l in open(f"{TRUTH}/sample{i}.truth.tsv"):
        if l.startswith("#"):
            continue
        g, ab = l.strip().split("\t")
        truth_acc.add(g)
        rows.append((g, float(ab)))
    sample_truth[i] = rows
all_g = [os.path.basename(p)[:-4] for p in glob.glob(f"{GEN}/*.fna")]
background = [g for g in all_g if g not in truth_acc]
contaminant = read_seq(sorted(glob.glob("multispecies/genomes/Escherichia_coli/*.fna"))[0])


def write_fasta(frags, path):
    with open(path, "w") as f:
        for i, fr in enumerate(frags):
            f.write(f">ctg{i}\n")
            for j in range(0, len(fr), 70):
                f.write(fr[j:j + 70] + "\n")


def evaluate(detected, truth_rows, present=0.01):
    """Strain (accession) level P/R + Bray-Curtis over version-less accessions."""
    truth = {strip_ver(g): a for g, a in truth_rows}
    det = {strip_ver(g): a for g, a in detected.items()}
    tset = {k for k, a in truth.items() if a >= present}
    dset = {k for k, a in det.items() if a >= present}
    tp, fp, fn = len(tset & dset), len(dset - tset), len(tset - dset)
    labels = set(truth) | set(det)
    num = sum(abs(det.get(k, 0) - truth.get(k, 0)) for k in labels)
    den = sum(det.get(k, 0) + truth.get(k, 0) for k in labels)
    bc = num / den if den else 0.0
    return tp, fp, fn, bc


out = open(f"{WORK}/refqual_degradation_strainscan.tsv", "w")
out.write("completeness\tcontamination\tn_contigs\tprecision\trecall\tbray_curtis\n")
for comp, contam, nc in LEVELS:
    tag = int(comp * 100)
    pdir = f"{WORK}/panel_{tag}"
    os.makedirs(pdir, exist_ok=True)
    for f in glob.glob(f"{pdir}/*.fna"):
        os.remove(f)
    for acc in truth_acc:
        if comp >= 1.0:
            shutil.copy(f"{GEN}/{acc}.fna", f"{pdir}/{acc}.fna")
        else:
            seq = read_seq(f"{GEN}/{acc}.fna")
            write_fasta(degrade(seq, comp, contam, nc, contaminant, sum(map(ord, acc)) + tag),
                        f"{pdir}/{acc}.fna")
    for g in background:
        shutil.copy(f"{GEN}/{g}.fna", f"{pdir}/{g}.fna")

    db = f"{WORK}/DB_{tag}"
    shutil.rmtree(db, ignore_errors=True)
    subprocess.run([PY, f"{SS}/StrainScan_build.py", "-i", pdir, "-o", db], check=True)

    TP = FP = FN = 0
    bcs = []
    for i in range(1, 6):
        od = f"{WORK}/out_{tag}_s{i}"
        subprocess.run([PY, f"{SS}/StrainScan.py", "-i", f"{READS}/sample{i}.fq",
                        "-d", db, "-o", od], check=True)
        detected = {}
        rep = f"{od}/final_report.txt"
        if os.path.exists(rep):
            lines = [l for l in open(rep) if not l.startswith("Strain_ID")]
            depths = {}
            for l in lines:
                c = l.rstrip("\n").split("\t")
                if len(c) > 5 and c[1]:
                    depths[c[1]] = depths.get(c[1], 0) + float(c[5] or 0)
            s = sum(depths.values()) or 1.0
            detected = {k: v / s for k, v in depths.items()}
        tp, fp, fn, bc = evaluate(detected, sample_truth[i])
        TP += tp; FP += fp; FN += fn; bcs.append(bc)
    P = TP / (TP + FP) if TP + FP else 0
    R = TP / (TP + FN) if TP + FN else 0
    BC = sum(bcs) / len(bcs)
    out.write(f"{comp:.2f}\t{contam:.2f}\t{nc}\t{P:.3f}\t{R:.3f}\t{BC:.3f}\n")
    print(f"completeness={comp:.2f} contam={contam:.2f}: P={P:.2f} R={R:.2f} BC={BC:.3f}", flush=True)
out.close()
print("done -> refqual_ss/refqual_degradation_strainscan.tsv "
      "(plot alongside the Strain2bScan arm with plot_refqual.py)")
