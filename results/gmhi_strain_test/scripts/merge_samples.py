#!/usr/bin/env python3
"""Merge all technical replicate FASTQ files per individual for PRJEB12449."""
import csv, os
from collections import defaultdict
from pathlib import Path

META = "/Volumes/MoneyCat/Data/GMHI/GMHI_master_metadata.tsv"
FASTQ_DIR = Path("/Volumes/MoneyCat/Data/GMHI/PRJEB12449")
OUT_DIR = Path("/Volumes/MoneyCat/Data/GMHI/strain_test/merged_reads")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def find_fastq(run_acc):
    """Find FASTQ file(s) for a run accession, handling single- and paired-end."""
    files = []
    # Single-end
    se = FASTQ_DIR / f"{run_acc}.fastq.gz"
    if se.exists() and not se.name.startswith("._"):
        files.append(se)
    # Paired-end
    pe1 = FASTQ_DIR / f"{run_acc}_1.fastq.gz"
    pe2 = FASTQ_DIR / f"{run_acc}_2.fastq.gz"
    if pe1.exists() and not pe1.name.startswith("._"):
        files.append(pe1)
    if pe2.exists() and not pe2.name.startswith("._"):
        files.append(pe2)
    return files

# Group run accessions by sample accession
samples = defaultdict(list)
with open(META) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        if r.get("study_accession") != "PRJEB12449":
            continue
        sample_acc = r.get("sample_accession", "")
        run_acc = r.get("run_accession", "")
        if not sample_acc or not run_acc:
            continue
        fqs = find_fastq(run_acc)
        if fqs:
            samples[sample_acc].append((run_acc, fqs, r.get("phenotype_binary", ""), r.get("phenotype_detail", "")))

print(f"Found {len(samples)} unique samples with FASTQ files")

# Pick 5 individuals: 2 Healthy, 1 CRC, 1 Overweight, 1 Obesity
picked = {}
need = {"Healthy": 2}
need_detail = {"CRC": 1, "Overweight": 1, "Obesity": 1}
for sample_acc, runs in samples.items():
    if sample_acc in picked:
        continue
    pheno_bin = runs[0][2]
    pheno_det = runs[0][3]
    if pheno_bin == "Healthy" and need["Healthy"] > 0:
        picked[sample_acc] = runs
        need["Healthy"] -= 1
    elif pheno_bin == "Nonhealthy" and pheno_det in need_detail and need_detail[pheno_det] > 0:
        picked[sample_acc] = runs
        need_detail[pheno_det] -= 1
    if need["Healthy"] == 0 and all(v == 0 for v in need_detail.values()):
        break

print(f"Selected {len(picked)} samples:")
for sample_acc, runs in picked.items():
    n_files = sum(len(fqs) for _, fqs, _, _ in runs)
    print(f"  {sample_acc}: {runs[0][3]} ({len(runs)} runs, {n_files} files)")

# Merge all FASTQ files per sample using cat
manifest = []
for sample_acc, runs in picked.items():
    out_fq = OUT_DIR / f"{sample_acc}.fastq.gz"
    fqs = [fq for _, fq_list, _, _ in runs for fq in fq_list]
    if out_fq.exists() and out_fq.stat().st_size > 0:
        print(f"SKIP {sample_acc} (exists)")
    else:
        print(f"Merging {sample_acc} ({len(fqs)} files) -> {out_fq}")
        with open(out_fq, "wb") as out:
            for fq in fqs:
                with open(fq, "rb") as src:
                    while True:
                        chunk = src.read(1024 * 1024)
                        if not chunk:
                            break
                        out.write(chunk)
    total_bytes = sum(fq.stat().st_size for _, fq_list, _, _ in runs for fq in fq_list)
    manifest.append((sample_acc, runs[0][2], runs[0][3], str(out_fq), len(runs), len(fqs), total_bytes))

# Write manifest
with open(OUT_DIR / "manifest.tsv", "w") as f:
    f.write("sample_accession\tphenotype_binary\tphenotype_detail\tmerged_fastq\tn_runs\tn_files\ttotal_bytes\n")
    for row in manifest:
        f.write("\t".join(map(str, row)) + "\n")

print(f"\nMerged reads in {OUT_DIR}")
