#!/usr/bin/env python3
"""Fig 4 (revised): impact of reference-genome completeness on strain identification, across ALL 15
species. For each species, hold the (simulated) sample reads fixed and degrade the truth strains'
reference genomes across a completeness ladder (100->50%, with co-varying contamination + fragmentation),
rebuild the DB, re-profile, and measure precision/recall/abundance error. Reuses the consistent
15-species pool + single-species simulated samples."""
import os, glob, subprocess, re, shutil, gzip, sys
sys.path.insert(0, "/Users/macstudio/Downloads/Strain2bScan-paper/scripts")
from degrade import degrade, read_seq

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
BIN = "/Users/macstudio/Downloads/Strain2bScan/target/release/strain2bscan"
POOL = f"{WORK}/sim/pool"
SIM = "/Users/macstudio/Downloads/Strain2bScan-paper/figure_raw_data/sim_single_species"
OUTDIR = f"{WORK}/sim/refqual15"; os.makedirs(OUTDIR, exist_ok=True)
SET14 = "CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"
LEVELS = [(1.00, 0.00, 1), (0.95, 0.01, 20), (0.90, 0.02, 50), (0.80, 0.05, 100), (0.70, 0.08, 200), (0.50, 0.10, 400)]
DEPTH = "5"   # fixed sample depth for the completeness sweep
CONTAINMENT = os.environ.get("CONTAINMENT", "") not in ("", "0")   # cluster on max-containment
CLUSTER_FLAGS = ["--containment"] if CONTAINMENT else []
OUTNAME = "refqual_15species_containment.tsv" if CONTAINMENT else "refqual_15species.tsv"
env = {**os.environ, "STRAIN2BSCAN_THREADS": "16"}

species = sorted(os.path.basename(d.rstrip("/")) for d in glob.glob(f"{POOL}/*/"))

def pick_config(sp):
    for strat, k in [("diff", 5), ("same", 5), ("diff", 3), ("same", 3), ("diff", 2), ("same", 2)]:
        if glob.glob(f"{SIM}/{sp}/reads/{sp}__{strat}_k{k}_rep1_d{DEPTH}_R1.fastq.gz"):
            return strat, k
    return None

def write_fasta(frags, path):
    with open(path, "w") as f:
        for i, fr in enumerate(frags):
            f.write(f">ctg{i}\n")
            for j in range(0, len(fr), 70): f.write(fr[j:j+70] + "\n")

def combined_reads(sp, name):
    """decompress R1+R2 into one temp fastq for profiling."""
    tmp = f"{OUTDIR}/_reads.fq"
    with open(tmp, "wb") as fo:
        for m in ("R1", "R2"):
            p = f"{SIM}/{sp}/reads/{name}_{m}.fastq.gz"
            with gzip.open(p, "rb") as fi: shutil.copyfileobj(fi, fo, 1 << 22)
    return tmp

out = open(f"{OUTDIR}/{OUTNAME}", "w")
out.write("species\tconfig\tn_truth_strains\tcompleteness\tcontamination\tn_contigs\tn_clusters\tprecision\trecall\tbray_curtis\n")
for sp in species:
    cfg = pick_config(sp)
    if not cfg: print(f"skip {sp} (no d{DEPTH} sample)", flush=True); continue
    strat, k = cfg
    # the 5 reps as the fixed sample set + their truth
    samples = []   # (name, [(acc, cluster, relab)])
    truth_acc = set()
    for rep in range(1, 6):
        name = f"{sp}__{strat}_k{k}_rep{rep}_d{DEPTH}"
        tf = f"{SIM}/{sp}/truth/{name}.truth.tsv"
        if not os.path.exists(tf): continue
        rows = []
        for l in open(tf):
            if l.startswith("#"): continue
            p = l.rstrip("\n").split("\t"); rows.append((p[0], float(p[2])))  # acc, relab
            truth_acc.add(p[0])
        samples.append((name, rows))
    panel = [os.path.basename(p)[:-4] for p in glob.glob(f"{POOL}/{sp}/*.fna")]
    background = [g for g in panel if g not in truth_acc]
    # cross-species contaminant (E. coli, or S. aureus if target is E. coli)
    csp = "Escherichia_coli" if sp != "Escherichia_coli" else "Staphylococcus_aureus"
    contaminant = read_seq(sorted(glob.glob(f"{POOL}/{csp}/*.fna"))[0])
    print(f"{sp}: {strat}_k{k}, {len(truth_acc)} truth strains, {len(background)} background, {len(samples)} samples", flush=True)
    for comp, contam, nc in LEVELS:
        tag = int(comp * 100); pdir = f"{OUTDIR}/{sp}_panel_{tag}"
        shutil.rmtree(pdir, ignore_errors=True); os.makedirs(pdir)
        for acc in truth_acc:
            if comp >= 1.0:
                shutil.copy(f"{POOL}/{sp}/{acc}.fna", f"{pdir}/{acc}.fna")
            else:
                seq = read_seq(f"{POOL}/{sp}/{acc}.fna")
                write_fasta(degrade(seq, comp, contam, nc, contaminant, sum(ord(c) for c in acc) + tag), f"{pdir}/{acc}.fna")
        for g in background:
            os.symlink(os.path.abspath(f"{POOL}/{sp}/{g}.fna"), f"{pdir}/{g}.fna")
        db = f"{OUTDIR}/{sp}_db_{tag}.tsv"
        r = subprocess.run([BIN, "cluster", "--genomes", pdir, "--enzyme", SET14, "--out", db, "--similarity", "0.95"] + CLUSTER_FLAGS,
                           capture_output=True, text=True, env=env)
        ncl = int(m.group(1)) if (m := re.search(r"into (\d+) cluster", r.stdout)) else 0
        g2c = {}
        for l in open(db.replace(".tsv", ".members.tsv")):
            if not l.startswith("#"):
                a, c = l.rstrip("\n").split("\t"); g2c[a] = c
        TP = FP = FN = 0; bcs = []
        for name, rows in samples:
            reads = combined_reads(sp, name)
            pred = f"{OUTDIR}/_pred.tsv"
            subprocess.run([BIN, "profile", "--db", db, "--reads", reads, "--out", pred], capture_output=True, text=True, env=env)
            ab = {}
            for acc, relab in rows:
                c = g2c.get(acc, acc); ab[c] = ab.get(c, 0) + relab
            tf = f"{OUTDIR}/_truth.tsv"; open(tf, "w").write("".join(f"{c}\t{a}\n" for c, a in ab.items()))
            ev = subprocess.run([BIN, "evaluate", "--pred", pred, "--truth", tf, "--present", "0.01"],
                                capture_output=True, text=True, env=env).stdout
            TP += int(re.search(r"TP=(\d+)", ev).group(1)); FP += int(re.search(r"FP=(\d+)", ev).group(1)); FN += int(re.search(r"FN=(\d+)", ev).group(1))
            bcs.append(float(re.search(r"Bray-Curtis=([0-9.]+)", ev).group(1)))
            os.remove(reads)
        P = TP / (TP + FP) if TP + FP else 0; R = TP / (TP + FN) if TP + FN else 0; BC = sum(bcs) / len(bcs) if bcs else 0
        out.write(f"{sp}\t{strat}_k{k}\t{len(truth_acc)}\t{comp:.2f}\t{contam:.2f}\t{nc}\t{ncl}\t{P:.3f}\t{R:.3f}\t{BC:.3f}\n"); out.flush()
        print(f"  comp={comp:.2f} contam={contam:.2f} contigs={nc}: clusters={ncl} P={P:.2f} R={R:.2f} BC={BC:.3f}", flush=True)
    # cleanup this species' panels + DBs to save disk
    for p in glob.glob(f"{OUTDIR}/{sp}_panel_*") + glob.glob(f"{OUTDIR}/{sp}_db_*"):
        shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) else os.remove(p)
print("REFQUAL15_DONE", flush=True); out.close()
