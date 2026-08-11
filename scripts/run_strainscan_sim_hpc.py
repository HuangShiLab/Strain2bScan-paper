#!/usr/bin/env python3
"""Build StrainScan DBs from the sim genome pool and profile the samples used in FigS5.

Designed to run on the HPC cluster as one or more SLURM array jobs:

  # build DBs for all species in the pool
  sbatch --array=0-14 strainscan_build_sim.slurm

  # profile all (sample, species) combinations
  sbatch --array=0-45 strainscan_profile_sim.slurm

The driver script itself is mode-agnostic: it reads SLURM_ARRAY_TASK_ID and
BUILD_MODE (0/1) from the environment.
"""
import os
import sys
import subprocess
from pathlib import Path

BASE = Path("/lustre1/g/aos_shihuang/LU/Strain2bScan-raw-data/sim_benchmark")
SSREPO = Path("/lustre1/g/aos_shihuang/tools/StrainScan")
PY = Path("/group/aos_shihuang/conda/envs/strainscan/bin/python")

DB_ROOT = BASE / "strainscan_dbs"
RES_ROOT = BASE / "strainscan_results"
POOL = BASE / "sim_genome_pool"

SINGLE_SAMPLES = [
    ("Escherichia_coli", "Escherichia_coli__diff_k2_rep1_d5"),
    ("Cutibacterium_acnes", "Cutibacterium_acnes__diff_k2_rep1_d5"),
    ("Staphylococcus_epidermidis", "Staphylococcus_epidermidis__diff_k2_rep1_d5"),
    ("Prevotella_copri", "Prevotella_copri__diff_k2_rep1_d5"),
]
MULTI_DEPTHS = ["depth_low", "depth_med", "depth_high"]


def all_species():
    return sorted(p.name for p in POOL.iterdir() if p.is_dir())


def build_one(species: str, threads: int = 8):
    DB_ROOT.mkdir(parents=True, exist_ok=True)
    out_dir = DB_ROOT / species
    work_dir = DB_ROOT / f"build_{species}"
    work_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(PY), str(SSREPO / "StrainScan_build.py"),
        "-i", str(POOL / species),
        "-o", str(out_dir),
        "-t", str(threads),
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SSREPO)
    print(f"[build {species}] {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, cwd=str(work_dir), check=True, env=env)
    print(f"[build {species}] done -> {out_dir}", flush=True)


def profile_one(sample_key: str, species: str, reads: tuple[str, str | None]):
    """Profile one sample against one species DB.

    sample_key is used to build the output directory tree, e.g.
    'single/Escherichia_coli__diff_k2_rep1_d5' or 'multi/depth_low_sample01'.
    """
    db_dir = DB_ROOT / species
    if not (db_dir / "Kmer_Sets_L1").exists():
        raise FileNotFoundError(f"DB missing for {species}: {db_dir}")
    out_dir = RES_ROOT / sample_key / species
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(PY), str(SSREPO / "StrainScan.py"),
        "-i", reads[0],
        "-d", str(db_dir),
        "-o", str(out_dir),
    ]
    if reads[1]:
        cmd += ["-j", reads[1]]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SSREPO)
    print(f"[profile {sample_key}/{species}] {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, cwd=str(out_dir), check=True, env=env)
    print(f"[profile {sample_key}/{species}] done", flush=True)


def single_reads(species: str, sample: str):
    rdir = BASE / "sim_single_species" / species / "reads"
    r1 = str(rdir / f"{sample}_R1.fastq.gz")
    r2 = str(rdir / f"{sample}_R2.fastq.gz")
    return (r1, r2)


def multi_reads(depth: str, sample: str):
    rdir = BASE / "sim_multi_species" / depth / "reads"
    r1 = str(rdir / f"{sample}_R1.fastq.gz")
    r2 = str(rdir / f"{sample}_R2.fastq.gz")
    return (r1, r2)


def build_mode():
    species = all_species()
    task = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
    if task >= len(species):
        print(f"task {task} >= n_species {len(species)}; nothing to do", flush=True)
        return
    build_one(species[task])


def profile_mode():
    species = all_species()
    jobs = []
    # multi-species: every species DB against every multi sample
    for depth in MULTI_DEPTHS:
        for sp in species:
            jobs.append((f"multi/{depth}_sample01", sp, multi_reads(depth, "sample01")))
    # single-species: only the four samples used in FigS5
    for sp, sample in SINGLE_SAMPLES:
        jobs.append((f"single/{sample}", sp, single_reads(sp, sample)))
    task = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
    if task >= len(jobs):
        print(f"task {task} >= n_jobs {len(jobs)}; nothing to do", flush=True)
        return
    sample_key, sp, reads = jobs[task]
    profile_one(sample_key, sp, reads)


if __name__ == "__main__":
    mode = os.environ.get("MODE", "profile")
    if mode == "build":
        build_mode()
    elif mode == "profile":
        profile_mode()
    else:
        raise ValueError(f"MODE must be build or profile, got {mode}")
