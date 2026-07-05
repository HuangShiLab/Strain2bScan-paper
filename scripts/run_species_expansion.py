#!/usr/bin/env python3
"""Expand the single-species head-to-head to 3 StrainScan-native species (Akkermansia
muciniphila, Prevotella copri, M. tuberculosis) using StrainScan's own pre-built DBs
(profiling only -- no dashing needed, works on macOS; run scripts/dl_strainscan_dbs.py first).
For each species:
  1. take a genome subset FROM StrainScan's own DB accession list (guarantees overlap,
     avoiding the earlier C.acnes DB-version-mismatch problem)
  2. resolve versions + download from NCBI
  3. build our own Strain2bScan cluster DB on that subset
  4. simulate 5 strain-mixture mock samples (log-normal depth >=1x, 2-5 truth strains)
  5. profile with REAL StrainScan (their pre-built DB, wall-clock capped -- see SS_TIMEOUT_S)
     and Strain2bScan (our DB)
  6. tabulate precision/recall/time/memory for both

NOTE on M. tuberculosis: in our run, StrainScan did not complete profiling a single sample
within >3.3 hours / >25.5 GB RAM against its own 792-strain DB (killed to protect the host;
see docs/species_expansion.md). SS_TIMEOUT_S below caps any repeat of that at a sane wall-clock
budget and records "DNF" instead of hanging.
"""
import os, sys, json, io, re, glob, random, shutil, subprocess, statistics, urllib.request, zipfile

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
BIN = os.environ.get("STRAIN2BSCAN_BIN", "strain2bscan")
SET14 = "CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"
SSPY = os.environ.get("STRAINSCAN_PY", "python")
SSREPO = os.environ.get("STRAINSCAN_REPO", "StrainScan")
SS_TIMEOUT_S = int(os.environ.get("SS_TIMEOUT_S", "600"))  # 10 min/sample safety cap
WORK = "species_expansion"
os.makedirs(WORK, exist_ok=True)

SPECIES = [
    dict(name="Akkermansia_muciniphila", short="akk", db="ss_db3/akk_db/DB_Akk",
         recls="ss_db3/akk_db/DB_Akk/Cluster_Result/hclsMap_95_recls.txt", n_subset=40, n_truth=5, seed=101),
    dict(name="Prevotella_copri", short="prev", db="ss_db3/prev_db/DB_Pre",
         recls="ss_db3/prev_db/DB_Pre/Cluster_Result/hclsMap_95_recls.txt", n_subset=40, n_truth=5, seed=102),
    dict(name="Mycobacterium_tuberculosis", short="mtb", db="ss_db3/mtb_db/DB_Mtb",
         recls="ss_db3/mtb_db/DB_Mtb/Cluster_Result/hclsMap_95_recls.txt", n_subset=40, n_truth=5, seed=103),
]

def get(url, timeout=90):
    req = urllib.request.Request(url, headers={"User-Agent": "s2bs-expansion"})
    return urllib.request.urlopen(req, timeout=timeout).read()

def resolve_version(acc):
    d = json.loads(get(f"{API}/genome/accession/{acc}/dataset_report"))
    r = d.get("reports", [])
    return r[0]["current_accession"] if r else None

def download_genome(versioned_acc, out):
    if os.path.exists(out):
        return True
    z = zipfile.ZipFile(io.BytesIO(get(f"{API}/genome/accession/{versioned_acc}/download?include_annotation_type=GENOME_FASTA")))
    fna = [n for n in z.namelist() if n.endswith("_genomic.fna")]
    if not fna:
        return False
    with z.open(fna[0]) as fh, open(out, "wb") as o:
        o.write(fh.read())
    return True

def read_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()

COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
def sim_reads(seq, depth, fh, tag, rng):
    L = len(seq); n = int(depth * L / 150)
    for i in range(n):
        p = rng.randint(0, max(0, L - 150)); r = seq[p:p+150]
        if rng.random() < 0.5:
            r = "".join(COMP.get(c, "N") for c in reversed(r))
        fh.write(f"@{tag}_{i}\n{r}\n+\n{'I'*150}\n")

def timed_run(cmd, env=None, timeout=None):
    """Returns (stdout, seconds, peak_rss_mb, ok). ok=False on timeout (killed)."""
    try:
        r = subprocess.run(["/usr/bin/time", "-l"] + cmd, capture_output=True, text=True, env=env, timeout=timeout)
    except subprocess.TimeoutExpired:
        return "", float(timeout), None, False
    t = re.search(r"([0-9.]+) real", r.stderr)
    m = re.search(r"(\d+)\s+maximum resident", r.stderr)
    return r.stdout, (float(t.group(1)) if t else 0), (int(m.group(1))/1048576 if m else 0), True

results = []
for sp in SPECIES:
    print(f"\n########## {sp['name']} ##########", flush=True)
    d = f"{WORK}/{sp['short']}"
    gdir = f"{d}/genomes"; os.makedirs(gdir, exist_ok=True)
    rng = random.Random(sp["seed"])

    # 1) pick subset of accessions from StrainScan's own DB list. Each line = one cluster;
    #    the last tab-separated field is a COMMA-separated list of that cluster's genomes
    #    (one genome per singleton cluster, several for merged clusters).
    all_accs = sorted(set(
        acc.strip()
        for line in open(sp["recls"])
        for acc in line.rstrip("\n").split("\t")[-1].split(",")
        if acc.strip().startswith("GC")
    ))
    subset = rng.sample(all_accs, min(sp["n_subset"], len(all_accs)))
    print(f"  panel: {len(all_accs)} genomes in StrainScan DB, using subset of {len(subset)}", flush=True)

    # 2) resolve + download
    for acc in subset:
        out = f"{gdir}/{acc}.fna"
        if os.path.exists(out):
            continue
        try:
            v = resolve_version(acc)
            if v:
                download_genome(v, out)
        except Exception as e:
            print(f"    {acc} FAIL: {e}", flush=True)
    have = [a for a in subset if os.path.exists(f"{gdir}/{a}.fna")]
    print(f"  downloaded {len(have)}/{len(subset)} genomes", flush=True)
    if len(have) < sp["n_truth"] + 5:
        print("  SKIP: too few genomes downloaded", flush=True); continue

    # 3) build Strain2bScan DB
    db_path = f"{d}/s2bscan.db.tsv"
    env16 = {**os.environ, "STRAIN2BSCAN_THREADS": "16"}
    bout, build_t, build_rss, _ = timed_run(
        [BIN, "cluster", "--genomes", gdir, "--enzyme", SET14, "--out", db_path, "--similarity", "0.95"], env16)
    ncl = int(re.search(r"into (\d+) cluster", bout).group(1)) if "into" in bout else 0
    print(f"  Strain2bScan build: {build_t:.2f}s, {build_rss:.0f}MB, {ncl} clusters", flush=True)
    g2c = {}
    for l in open(db_path.replace(".tsv", ".members.tsv")):
        if not l.startswith("#"):
            g, c = l.strip().split("\t"); g2c[g] = c

    # 4) simulate 5 samples
    truth_pool = rng.sample(have, sp["n_truth"])
    sdir = f"{d}/samples"; tdir = f"{d}/truth"
    os.makedirs(sdir, exist_ok=True); os.makedirs(tdir, exist_ok=True)
    for i in range(1, 6):
        k = rng.randint(2, min(5, len(truth_pool)))
        chosen = rng.sample(truth_pool, k)
        depths = [max(1.0, round(rng.lognormvariate(1.0, 0.6), 2)) for _ in chosen]
        tot = sum(depths)
        with open(f"{sdir}/sample{i}.fq", "w") as fq:
            for acc, dep in zip(chosen, depths):
                sim_reads(read_seq(f"{gdir}/{acc}.fna"), dep, fq, acc, rng)
        with open(f"{tdir}/sample{i}.truth.tsv", "w") as t:
            t.write("#genome\tabundance\n")
            for acc, dep in zip(chosen, depths):
                t.write(f"{acc}\t{dep/tot:.4f}\n")

    # 5a) profile with Strain2bScan
    TP=FP=FN=0; s2b_times=[]; s2b_rss=[]
    for i in range(1, 6):
        pred = f"{d}/s2b_s{i}.pred"
        po, pt, prss, _ = timed_run([BIN, "profile", "--db", db_path, "--reads", f"{sdir}/sample{i}.fq", "--out", pred], env16)
        s2b_times.append(pt); s2b_rss.append(prss)
        ab = {}
        for l in open(f"{tdir}/sample{i}.truth.tsv"):
            if l.startswith("#"): continue
            g, a = l.strip().split("\t"); c = g2c.get(g, g); ab[c] = ab.get(c, 0) + float(a)
        truth_f = f"{d}/s2b_s{i}.truth"
        open(truth_f, "w").write("".join(f"{c}\t{a}\n" for c, a in ab.items()))
        ev = subprocess.run([BIN, "evaluate", "--pred", pred, "--truth", truth_f, "--present", "0.01"],
                            capture_output=True, text=True).stdout
        TP += int(re.search(r"TP=(\d+)", ev).group(1)); FP += int(re.search(r"FP=(\d+)", ev).group(1)); FN += int(re.search(r"FN=(\d+)", ev).group(1))
    s2b_P = TP/(TP+FP) if TP+FP else 0; s2b_R = TP/(TP+FN) if TP+FN else 0

    # 5b) profile with REAL StrainScan (their own pre-built DB), wall-clock capped per sample
    ss_times=[]; ss_rss=[]; ss_TP=ss_FP=ss_FN=0; ss_dnf=False
    for i in range(1, 6):
        od = f"{d}/ss_out_s{i}"
        shutil.rmtree(od, ignore_errors=True)
        _, t, m, ok = timed_run([SSPY, f"{SSREPO}/StrainScan.py", "-i", f"{sdir}/sample{i}.fq",
                                  "-d", sp["db"], "-o", od], timeout=SS_TIMEOUT_S)
        if not ok:
            print(f"  sample{i}: StrainScan DNF (killed after {SS_TIMEOUT_S}s)", flush=True)
            ss_dnf = True
            ss_times.append(SS_TIMEOUT_S)
            continue
        ss_times.append(t); ss_rss.append(m if m is not None else 0)
        rep = f"{od}/final_report.txt"
        detected = set()
        if os.path.exists(rep):
            for l in open(rep).read().splitlines()[1:]:
                f = l.split("\t")
                if len(f) > 1 and f[1]:
                    detected.add(f[1])
        truth_accs = set(l.split("\t")[0] for l in open(f"{tdir}/sample{i}.truth.tsv") if not l.startswith("#"))
        ss_TP += len(detected & truth_accs); ss_FP += len(detected - truth_accs); ss_FN += len(truth_accs - detected)

    if ss_dnf:
        ss_profile_s, ss_rss_mb, ss_P_s, ss_R_s = f"DNF(>{SS_TIMEOUT_S}s)", "DNF", "DNF", "DNF"
    else:
        ss_P = ss_TP/(ss_TP+ss_FP) if ss_TP+ss_FP else 0; ss_R = ss_TP/(ss_TP+ss_FN) if ss_TP+ss_FN else 0
        ss_profile_s, ss_rss_mb = round(statistics.mean(ss_times),2), round(statistics.mean(ss_rss),0)
        ss_P_s, ss_R_s = round(ss_P,3), round(ss_R,3)

    row = dict(species=sp["name"], n_genomes_subset=len(have), n_clusters=ncl,
               build_s=round(build_t,2), build_rss_mb=round(build_rss,0),
               s2b_profile_s=round(statistics.mean(s2b_times),2), s2b_rss_mb=round(statistics.mean(s2b_rss),0),
               s2b_precision=round(s2b_P,3), s2b_recall=round(s2b_R,3),
               ss_profile_s=ss_profile_s, ss_rss_mb=ss_rss_mb, ss_precision=ss_P_s, ss_recall=ss_R_s)
    results.append(row)
    print(f"  Strain2bScan: {row['s2b_profile_s']}s/{row['s2b_rss_mb']}MB P={row['s2b_precision']} R={row['s2b_recall']}", flush=True)
    print(f"  StrainScan:   {row['ss_profile_s']}s/{row['ss_rss_mb']}MB P={row['ss_precision']} R={row['ss_recall']}", flush=True)

cols = list(results[0].keys()) if results else []
with open(f"{WORK}/species_expansion.tsv", "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in results:
        o.write("\t".join(str(r[c]) for c in cols) + "\n")
print(f"\ndone -> {WORK}/species_expansion.tsv")
