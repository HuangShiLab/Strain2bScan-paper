#!/usr/bin/env python3
"""Exploratory profiling of clinical oral 2bRAD samples (PRJNA1131785 S_# cohort) vs the oral panel.
No case/control labels are available in the public metadata, so this is a demonstration that
Strain2bScan resolves strain-level profiles on the clinical cohort, not a differential test."""
import subprocess, os, glob, gzip, time

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
BIN  = "/Users/macstudio/Downloads/Strain2bScan/target/release/strain2bscan"
DBS  = f"{WORK}/dbs/oral"

def profile(reads):
    tmp = f"{WORK}/results/_tmp_clin.fq"
    with gzip.open(reads, "rb") as fi, open(tmp, "wb") as fo:
        while True:
            b = fi.read(1 << 22)
            if not b: break
            fo.write(b)
    t0 = time.time()
    out = subprocess.run([BIN, "multi-profile", "--dbs", DBS, "--reads", tmp, "--enzyme", "BcgI",
                          "--min-species-markers", "50", "--min-species-detect", "5"],
                         capture_output=True, text=True, env={**os.environ, "STRAIN2BSCAN_THREADS": "16"})
    dt = time.time() - t0; os.remove(tmp)
    return out.stdout, dt

rows = [["sample", "markers", "species_resolved", "strain_calls", "runtime_s", "top_species"]]
print(f"{'sample':14s} {'markers':>9s} {'species':>7s} {'strains':>7s} {'t_s':>5s}  top species")
for f in sorted(glob.glob(f"{WORK}/reads/clinical_2brad/*_1.fastq.gz")):
    alias = os.path.basename(f).split("__")[0]
    so, dt = profile(f)
    nmark = "?"; species = {}; nstrain = 0
    for line in so.splitlines():
        if line.startswith("sample:"): nmark = line.split()[1]
        if not line.startswith("  "): continue
        p = line.strip().split("\t")
        if "[detected" in line or "[strain-resolved, no cluster" in line: continue
        if len(p) >= 5:
            nstrain += 1; species[p[0]] = species.get(p[0], 0.0) + float(p[4])
    top = sorted(species, key=species.get, reverse=True)[:3]
    print(f"{alias:14s} {nmark:>9s} {len(species):>7d} {nstrain:>7d} {dt:>5.1f}  {', '.join(s.replace('_',' ') for s in top)}")
    rows.append([alias, nmark, len(species), nstrain, f"{dt:.1f}", "; ".join(top)])
with open(f"{WORK}/results/clinical_exploratory.tsv", "w") as w:
    for r in rows: w.write("\t".join(str(x) for x in r) + "\n")
print(f"\nwrote results/clinical_exploratory.tsv")
