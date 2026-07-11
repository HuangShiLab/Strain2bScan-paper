#!/usr/bin/env python3
"""Recompute between-strain distance concordance with whole-genome, for 2bRAD tags AND 16S, across
the 15 slide-deck species. Whole-genome = bottom-3000 canonical 21-mer MinHash (Mash). 2bRAD =
Strain2bScan `build` BcgI tag sets. 16S = longest barrnap-extracted 16S gene, exact canonical
21-mer Jaccard. All distances use the Mash transform D(J) = -ln(2J/(1+J)); per species we correlate
the pairwise distance vector against whole-genome (Spearman). Run with PYTHONHASHSEED=0."""
import subprocess, os, glob, math, heapq, itertools, sys
from concurrent.futures import ThreadPoolExecutor

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
BIN = os.path.expanduser("~/Downloads/Strain2bScan/target/release/strain2bscan")
BARRNAP = f"{WORK}/barrnap/bin/barrnap"
GEN = f"{WORK}/genomes_16s"
OUT_TSV = os.path.expanduser("~/Downloads/Strain2bScan-paper/results/mash_2brad_vs_16s_recomputed.tsv")
K, SKETCH = 21, 3000
COMP = str.maketrans("ACGT", "TGCA")

def read_fna(path):
    seqs, cur = [], []
    for line in open(path):
        if line.startswith(">"):
            if cur: seqs.append("".join(cur)); cur = []
        else: cur.append(line.strip())
    if cur: seqs.append("".join(cur))
    return seqs

def canon_kmers(seq, exact=False):
    seq = seq.upper(); hs = set()
    for i in range(len(seq) - K + 1):
        km = seq[i:i+K]
        if not exact and ("N" in km): continue
        rc = km.translate(COMP)[::-1]
        c = km if km < rc else rc
        hs.add(hash(c) & 0xFFFFFFFFFFFFFFFF)
    return hs

def wgs_sketch(path):
    hs = set()
    for seq in read_fna(path): hs |= canon_kmers(seq)
    return set(heapq.nsmallest(SKETCH, hs))

def sixteen_s_kmers(path):
    """barrnap -> longest 16S copy -> canonical 21-mer set (exact)."""
    out = subprocess.run([BARRNAP, "--quiet", "--kingdom", "bac", path],
                         capture_output=True, text=True,
                         env={**os.environ, "PATH": f"{os.path.dirname(BARRNAP)}:" + os.environ["PATH"]})
    # parse GFF for 16S features, extract sequence from the genome
    contigs = {}
    name = None
    for line in open(path):
        if line.startswith(">"): name = line[1:].split()[0]; contigs[name] = []
        else: contigs[name].append(line.strip())
    contigs = {k: "".join(v) for k, v in contigs.items()}
    best = ""
    for line in out.stdout.splitlines():
        if line.startswith("#") or "\t" not in line: continue
        f = line.split("\t")
        if len(f) < 9 or "16S" not in f[8]: continue
        ctg, s, e, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
        seq = contigs.get(ctg, "")[s:e]
        if strand == "-": seq = seq.upper().translate(COMP)[::-1]
        if len(seq) > len(best): best = seq
    return canon_kmers(best, exact=True) if len(best) > 1000 else None

def jacc_minhash(a, b):
    union = heapq.nsmallest(SKETCH, a | b)
    return sum(1 for h in union if h in a and h in b) / len(union) if union else 0.0

def jacc(a, b):
    return len(a & b) / len(a | b) if (a | b) else 0.0

def dist(J):
    if J <= 0: return None
    x = 2 * J / (1 + J)
    return -math.log(x) if 0 < x < 1 else 0.0

def build_tags(gdir, out):
    subprocess.run([BIN, "build", "--genomes", gdir, "--enzyme", "BcgI", "--out", out],
                   capture_output=True, text=True, env={**os.environ, "STRAIN2BSCAN_THREADS": "16"})
    tags = {}
    for line in open(out):
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[1]: tags[p[0]] = set(p[1].split(","))
    return tags

def spearman(x, y):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0]*len(v); i = 0
        while i < len(v):
            j = i
            while j+1 < len(v) and v[order[j+1]] == v[order[i]]: j += 1
            avg = (i+j)/2 + 1
            for t in range(i, j+1): r[order[t]] = avg
            i = j+1
        return r
    xr, yr = rank(x), rank(y); n = len(x)
    mx, my = sum(xr)/n, sum(yr)/n
    sx = sum((a-mx)**2 for a in xr); sy = sum((b-my)**2 for b in yr)
    sxy = sum((a-mx)*(b-my) for a, b in zip(xr, yr))
    return sxy/math.sqrt(sx*sy) if sx > 0 and sy > 0 else 0.0

species_dirs = sorted(glob.glob(f"{GEN}/*/"))
print("species\tn_genomes\tn_pairs\tspearman_wgs_2brad\tspearman_wgs_16s\tn_16s_found", flush=True)
rows = []
for sd in species_dirs:
    sp = os.path.basename(sd.rstrip("/"))
    genomes = sorted(glob.glob(f"{sd}/*.fna"))
    if len(genomes) < 8: print(f"# skip {sp} ({len(genomes)} genomes)", flush=True); continue
    accs = [os.path.basename(g)[:-4] for g in genomes]
    tags = build_tags(sd, f"{WORK}/16s_tags_{sp}.tsv")
    with ThreadPoolExecutor(max_workers=8) as ex:
        wgs = dict(zip(accs, ex.map(wgs_sketch, genomes)))
        s16 = dict(zip(accs, ex.map(sixteen_s_kmers, genomes)))
    n16 = sum(1 for a in accs if s16.get(a))
    common = [a for a in accs if a in tags and wgs.get(a) and s16.get(a)]
    dw2, d2, dw16, d16 = [], [], [], []
    for a, b in itertools.combinations(common, 2):
        Jw = jacc_minhash(wgs[a], wgs[b]); Dw = dist(Jw)
        Ta, Tb = tags[a], tags[b]; D2 = dist(jacc(Ta, Tb))
        D16 = dist(jacc(s16[a], s16[b]))
        if Dw is None: continue
        if D2 is not None: dw2.append(Dw); d2.append(D2)
        if D16 is not None: dw16.append(Dw); d16.append(D16)
    sp2 = spearman(dw2, d2) if len(dw2) > 2 else float("nan")
    sp16 = spearman(dw16, d16) if len(dw16) > 2 else float("nan")
    row = [sp, len(common), len(dw2), round(sp2, 3), round(sp16, 3), n16]
    print("\t".join(str(x) for x in row), flush=True)
    rows.append(row)

os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
with open(OUT_TSV, "w") as f:
    f.write("# 2bRAD-tag AND 16S vs whole-genome between-strain distance concordance, recomputed with\n")
    f.write("# Strain2bScan. WGS: bottom-3000 canonical 21-mer MinHash (Mash). 2bRAD: Strain2bScan build\n")
    f.write("# BcgI tag sets. 16S: longest barrnap 16S gene, exact canonical 21-mer Jaccard. distance\n")
    f.write("# = -ln(2J/(1+J)); Spearman of pairwise distance vector vs whole-genome. PYTHONHASHSEED=0.\n")
    f.write("species\tn_genomes\tn_pairs\tspearman_wgs_2brad\tspearman_wgs_16s\tn_16s_found\n")
    for r in rows: f.write("\t".join(str(x) for x in r) + "\n")
print(f"\n# wrote {OUT_TSV}", flush=True)
