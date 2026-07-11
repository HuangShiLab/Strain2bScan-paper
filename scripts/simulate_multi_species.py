#!/usr/bin/env python3
"""Multi-species multi-strain simulated community from the SAME 15-species pool, across 3 depth
gradients (low/med/high). Each sample contains all 15 species; per species 1-2 strains at
log-normal relative abundance; ART PE250. Outputs paired FASTQ + truth (species, strain, cluster,
relative abundance, coverage) into figure_raw_data/sim_multi_species/depth_<grad>/. Resumable."""
import os, glob, random, subprocess, shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
ART = f"{WORK}/mamba/envs/art/bin/art_illumina"
POOL = f"{WORK}/sim/pool"
OUT = "/Users/macstudio/Downloads/Strain2bScan-paper/figure_raw_data/sim_multi_species"
TMP = f"{WORK}/sim/_arttmp_ms"
GRADIENTS = {"low": 1.0, "med": 5.0, "high": 10.0}   # mean per-strain coverage
N = int(os.environ.get("NSAMP", "20"))
os.makedirs(TMP, exist_ok=True)

clusters = {}   # species -> {cluster: [acc]}
for line in open(f"{WORK}/sim/pool_manifest.tsv").read().splitlines()[1:]:
    sp, acc, cl = line.split("\t")
    clusters.setdefault(sp, {}).setdefault(cl, []).append(acc)
acc2cl = {}
for sp, cl in clusters.items():
    for c, accs in cl.items():
        for a in accs: acc2cl[(sp, a)] = c
species = sorted(clusters)

def gz(src, dst):
    with open(dst, "wb") as fo:
        subprocess.run(["gzip", "-3", "-c", src], stdout=fo)

def sample(grad, idx):
    depth = GRADIENTS[grad]
    rng = random.Random(hash((grad, idx)) & 0xFFFFFFFF)
    rdir = f"{OUT}/depth_{grad}/reads"; tdir = f"{OUT}/depth_{grad}/truth"
    os.makedirs(rdir, exist_ok=True); os.makedirs(tdir, exist_ok=True)
    name = f"sample{idx:02d}"
    r1 = f"{rdir}/{name}_R1.fastq.gz"
    if os.path.exists(r1) and os.path.getsize(r1) > 0: return 1
    # pick strains: all species, 1-2 strains each; weight ~ lognormal
    picks = []   # (sp, acc, weight)
    for sp in species:
        allg = [a for c in clusters[sp].values() for a in c]
        k = 2 if (len(allg) >= 2 and rng.random() < 0.3) else 1
        for a in rng.sample(allg, k):
            picks.append((sp, a, rng.lognormvariate(0, 1.0)))
    tot = sum(w for _, _, w in picks); nstr = len(picks)
    cat1 = f"{TMP}/{grad}_{name}_1.fq"; cat2 = f"{TMP}/{grad}_{name}_2.fq"
    open(cat1, "w").close(); open(cat2, "w").close()
    truth = []
    for j, (sp, acc, w) in enumerate(picks):
        ra = w / tot; cov = ra * depth * nstr
        if cov < 0.02: cov = 0.02
        pref = f"{TMP}/{grad}_{name}_s{j}"
        subprocess.run([ART, "-p", "-l", "250", "-f", f"{cov:.4f}", "-m", "600", "-s", "150",
                        "-i", f"{POOL}/{sp}/{acc}.fna", "-o", pref, "-na", "-q", "-rs", str(2000 + j)],
                       capture_output=True)
        for a, cat in ((f"{pref}1.fq", cat1), (f"{pref}2.fq", cat2)):
            if os.path.exists(a):
                with open(a) as fi, open(cat, "a") as fo: shutil.copyfileobj(fi, fo, 1 << 22)
                os.remove(a)
        truth.append((sp, acc, acc2cl[(sp, acc)], f"{ra:.5f}", f"{cov:.4f}"))
    gz(cat1, r1); gz(cat2, f"{rdir}/{name}_R2.fastq.gz"); os.remove(cat1); os.remove(cat2)
    with open(f"{tdir}/{name}.truth.tsv", "w") as t:
        t.write("#species\tstrain\tcluster\trelative_abundance\tcoverage\tmean_depth\n")
        for sp, acc, cl, ra, cov in truth:
            t.write(f"{sp}\t{acc}\t{cl}\t{ra}\t{cov}\t{depth}\n")
    return 1

tasks = [(g, i) for g in GRADIENTS for i in range(1, N + 1)]
done = 0
with ThreadPoolExecutor(max_workers=8) as ex:
    futs = [ex.submit(sample, g, i) for g, i in tasks]
    for f in as_completed(futs):
        done += f.result()
        if done % 10 == 0: print(f"  ...{done}/{len(tasks)} samples", flush=True)
print(f"MULTI_SPECIES_DONE  {done} samples ({N}/gradient x 3) -> {OUT}", flush=True)
