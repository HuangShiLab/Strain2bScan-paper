#!/usr/bin/env python3
"""Single-species multi-strain simulated metagenomes (StrainScan-style), one consistent dataset.
Per species (15): 2/3/5-strain communities with StrainScan uneven coverage profiles, under two
difficulty strategies (strains from DIFFERENT clusters / from the SAME cluster), N reps, each
simulated at a depth ladder (0.5/1/3/5/10x) with ART (art_illumina -p -l 250 -m 600 -s 150).
The SAME community (strain set + relative abundances) is simulated across all depths.
Outputs paired FASTQ + per-sample truth (strain, cluster, relative abundance, coverage) into
figure_raw_data/sim_single_species/<Species>/. Resumable (skips existing R1.gz)."""
import os, glob, random, subprocess, gzip, shutil, sys
from concurrent.futures import ThreadPoolExecutor, as_completed

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
ART = f"{WORK}/mamba/envs/art/bin/art_illumina"
POOL = f"{WORK}/sim/pool"
OUT = "/Users/macstudio/Downloads/Strain2bScan-paper/figure_raw_data/sim_single_species"
TMP = f"{WORK}/sim/_arttmp"
PROFILES = {2: [100, 10], 3: [100, 50, 10], 5: [100, 70, 50, 20, 10]}
DEPTHS = [0.5, 1, 3, 5, 10]
REPS = int(os.environ.get("REPS", "5"))
os.makedirs(TMP, exist_ok=True)

# cluster membership per species
clusters = {}   # species -> {cluster: [acc,...]}
for line in open(f"{WORK}/sim/pool_manifest.tsv").read().splitlines()[1:]:
    sp, acc, cl = line.split("\t")
    clusters.setdefault(sp, {}).setdefault(cl, []).append(acc)

def pick_diff(cl, k, rng):
    cs = [c for c in cl if cl[c]]
    if len(cs) < k: return None
    chosen_cl = rng.sample(cs, k)
    return [(rng.choice(cl[c]), c) for c in chosen_cl]

def pick_same(cl, k, rng):
    big = [c for c in cl if len(cl[c]) >= k]
    if not big: return None
    c = rng.choice(big)
    return [(a, c) for a in rng.sample(cl[c], k)]

def gz(src, dst):
    # system gzip -3 is much faster than python gzip level 6
    with open(dst, "wb") as fo:
        subprocess.run(["gzip", "-3", "-c", src], stdout=fo)

def simulate(sp, strat, k, rep):
    cl = clusters[sp]
    rng = random.Random(hash((sp, strat, k, rep)) & 0xFFFFFFFF)
    strains = (pick_diff if strat == "diff" else pick_same)(cl, k, rng)
    if strains is None: return 0
    prof = PROFILES[k]; tot = sum(prof)
    relab = [p / tot for p in prof]
    rdir = f"{OUT}/{sp}/reads"; tdir = f"{OUT}/{sp}/truth"
    os.makedirs(rdir, exist_ok=True); os.makedirs(tdir, exist_ok=True)
    made = 0
    for D in DEPTHS:
        name = f"{sp}__{strat}_k{k}_rep{rep}_d{D}"
        r1 = f"{rdir}/{name}_R1.fastq.gz"; r2 = f"{rdir}/{name}_R2.fastq.gz"
        if os.path.exists(r1) and os.path.getsize(r1) > 0:
            made += 1; continue
        cat1 = f"{TMP}/{name}_1.fq"; cat2 = f"{TMP}/{name}_2.fq"
        open(cat1, "w").close(); open(cat2, "w").close()
        truth = []
        for j, ((acc, clid), ra) in enumerate(zip(strains, relab)):
            cov = ra * D * k                      # mean per-strain coverage = D
            g = f"{POOL}/{sp}/{acc}.fna"
            pref = f"{TMP}/{name}_s{j}"
            subprocess.run([ART, "-p", "-l", "250", "-f", f"{cov:.4f}", "-m", "600", "-s", "150",
                            "-i", g, "-o", pref, "-na", "-q", "-rs", str(1000 + j)],
                           capture_output=True)
            for a, cat in ((f"{pref}1.fq", cat1), (f"{pref}2.fq", cat2)):
                if os.path.exists(a):
                    with open(a) as fi, open(cat, "a") as fo: shutil.copyfileobj(fi, fo, 1 << 22)
                    os.remove(a)
            truth.append((acc, clid, f"{ra:.4f}", f"{cov:.4f}"))
        gz(cat1, r1); gz(cat2, r2); os.remove(cat1); os.remove(cat2)
        with open(f"{tdir}/{name}.truth.tsv", "w") as t:
            t.write("#strain\tcluster\trelative_abundance\tcoverage\tmean_depth\tn_strains\tstrategy\n")
            for acc, clid, ra, cov in truth:
                t.write(f"{acc}\t{clid}\t{ra}\t{cov}\t{D}\t{k}\t{strat}\n")
        made += 1
    return made

tasks = [(sp, strat, k, rep) for sp in sorted(clusters)
         for k in (2, 3, 5) for strat in ("diff", "same") for rep in range(1, REPS + 1)]
total = 0; done = 0
with ThreadPoolExecutor(max_workers=8) as ex:
    futs = {ex.submit(simulate, *t): t for t in tasks}
    for f in as_completed(futs):
        sp, strat, k, rep = futs[f]
        total += f.result(); done += 1
        if done % 25 == 0:
            print(f"  ...{done}/{len(tasks)} configs, {total} samples", flush=True)
print(f"SINGLE_SPECIES_DONE  {total} samples from {len(tasks)} configs -> {OUT}", flush=True)
