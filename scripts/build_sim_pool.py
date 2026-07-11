#!/usr/bin/env python3
"""Pin the 15-species HQ genome pool and build per-species Strain2bScan cluster DBs (14-enzyme, 0.95),
extracting cluster membership for same/different-cluster strain selection in the simulation."""
import os, glob, subprocess, shutil

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
BIN = "/Users/macstudio/Downloads/Strain2bScan/target/release/strain2bscan"
QC = "/Users/macstudio/Downloads/Strain2bScan-paper/data/genome_qc_16s_panel.tsv"
GEN = f"{WORK}/genomes_16s"
POOL = f"{WORK}/sim/pool"; DBS = f"{WORK}/sim/dbs"
SET14 = "CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"
os.makedirs(POOL, exist_ok=True); os.makedirs(DBS, exist_ok=True)

# HQ accessions per species
hq = {}
for line in open(QC).read().splitlines()[1:]:
    f = line.split("\t")
    if f[7] == "1":
        hq.setdefault(f[0], set()).add(f[1])   # species -> {GCA_num}

manifest = ["species\taccession\tcluster"]
summary = ["species\tn_hq_genomes\tn_clusters"]
for sp in sorted(hq):
    src = f"{GEN}/{sp}"
    genomes = [g for g in sorted(glob.glob(f"{src}/*.fna")) if os.path.basename(g)[:-4] in hq[sp]]
    if len(genomes) < 5:
        print(f"skip {sp} ({len(genomes)} HQ)"); continue
    pdir = f"{POOL}/{sp}"; os.makedirs(pdir, exist_ok=True)
    for g in glob.glob(f"{pdir}/*.fna"): os.remove(g)
    for g in genomes: os.symlink(g, f"{pdir}/{os.path.basename(g)}")
    db = f"{DBS}/{sp}.tsv"
    subprocess.run([BIN, "cluster", "--genomes", pdir, "--enzyme", SET14, "--out", db, "--similarity", "0.95"],
                   capture_output=True, text=True, env={**os.environ, "STRAIN2BSCAN_THREADS": "16"})
    g2c = {}
    mem = db.replace(".tsv", ".members.tsv")
    for l in open(mem):
        if l.startswith("#"): continue
        parts = l.rstrip("\n").split("\t")
        if len(parts) >= 2: g2c[parts[0]] = parts[1]
    ncl = len(set(g2c.values()))
    for g in genomes:
        acc = os.path.basename(g)[:-4]
        manifest.append(f"{sp}\t{acc}\t{g2c.get(acc,'?')}")
    summary.append(f"{sp}\t{len(genomes)}\t{ncl}")
    print(f"{sp}: {len(genomes)} HQ genomes -> {ncl} clusters", flush=True)

open(f"{WORK}/sim/pool_manifest.tsv", "w").write("\n".join(manifest) + "\n")
open(f"{WORK}/sim/pool_summary.tsv", "w").write("\n".join(summary) + "\n")
print("\nwrote sim/pool_manifest.tsv + sim/pool_summary.tsv")
