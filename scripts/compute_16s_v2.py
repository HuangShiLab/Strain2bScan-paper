#!/usr/bin/env python3
"""16S-vs-2bRAD concordance, v2: restrict to COMPLETE / near-complete genomes (CheckM completeness
>=97%, contamination <=5%, assembly level Complete Genome/Chromosome) so incomplete assemblies do not
bias Mash distances, and add genome-level subsampling bootstrap CIs on the Spearman correlations.
Run with PYTHONHASHSEED=0."""
import subprocess, os, glob, math, heapq, itertools, random
from concurrent.futures import ThreadPoolExecutor

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
BIN = os.path.expanduser("~/Downloads/Strain2bScan/target/release/strain2bscan")
BARRNAP = f"{WORK}/barrnap/bin/barrnap"
GEN = f"{WORK}/genomes_16s"
OUT_TSV = os.path.expanduser("~/Downloads/Strain2bScan-paper/results/mash_2brad_vs_16s_recomputed.tsv")
K, SKETCH = 21, 3000
COMP = str.maketrans("ACGT", "TGCA")
random.seed(0)

# high-quality genome set per species
HQ = {}
for line in open(f"{WORK}/genome_qc.tsv").read().splitlines()[1:]:
    f = line.split("\t")
    if f[7] == "1": HQ.setdefault(f[0], set()).add(f[1])

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
        if not exact and "N" in km: continue
        rc = km.translate(COMP)[::-1]; c = km if km < rc else rc
        hs.add(hash(c) & 0xFFFFFFFFFFFFFFFF)
    return hs

def wgs_sketch(path):
    hs = set()
    for seq in read_fna(path): hs |= canon_kmers(seq)
    return set(heapq.nsmallest(SKETCH, hs))

def sixteen_s(path):
    out = subprocess.run([BARRNAP, "--quiet", "--kingdom", "bac", path], capture_output=True, text=True,
                         env={**os.environ, "PATH": f"{os.path.dirname(BARRNAP)}:" + os.environ["PATH"]})
    contigs = {}; name = None
    for line in open(path):
        if line.startswith(">"): name = line[1:].split()[0]; contigs[name] = []
        else: contigs[name].append(line.strip())
    contigs = {k: "".join(v) for k, v in contigs.items()}
    best = ""
    for line in out.stdout.splitlines():
        if line.startswith("#") or "\t" not in line: continue
        c = line.split("\t")
        if len(c) < 9 or "16S" not in c[8]: continue
        seq = contigs.get(c[0], "")[int(c[3])-1:int(c[4])]
        if c[6] == "-": seq = seq.upper().translate(COMP)[::-1]
        if len(seq) > len(best): best = seq
    return canon_kmers(best, exact=True) if len(best) > 1000 else None

def jacc_minhash(a, b):
    u = heapq.nsmallest(SKETCH, a | b)
    return sum(1 for h in u if h in a and h in b) / len(u) if u else 0.0

def jacc(a, b): return len(a & b) / len(a | b) if (a | b) else 0.0

def dist(J):
    if J <= 0: return None
    x = 2*J/(1+J); return -math.log(x) if 0 < x < 1 else 0.0

def build_tags(gdir, accs, out):
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
            avg = (i+j)/2+1
            for t in range(i, j+1): r[order[t]] = avg
            i = j+1
        return r
    xr, yr = rank(x), rank(y); n = len(x); mx, my = sum(xr)/n, sum(yr)/n
    sx = sum((a-mx)**2 for a in xr); sy = sum((b-my)**2 for b in yr)
    sxy = sum((a-mx)*(b-my) for a, b in zip(xr, yr))
    return sxy/math.sqrt(sx*sy) if sx > 0 and sy > 0 else float("nan")

def boot_ci(idx_pairs, dw, dcmp, genomes_n, B=500, frac=0.8):
    """Subsample genomes (frac, without replacement) B times; recompute Spearman -> 2.5/97.5 pctile."""
    m = max(6, int(frac*genomes_n)); vals = []
    pair_of = {}
    for k, (i, j) in enumerate(idx_pairs): pair_of[(i, j)] = k
    for _ in range(B):
        sub = set(random.sample(range(genomes_n), m))
        xs, ys = [], []
        for k, (i, j) in enumerate(idx_pairs):
            if i in sub and j in sub: xs.append(dw[k]); ys.append(dcmp[k])
        if len(xs) > 3:
            s = spearman(xs, ys)
            if not math.isnan(s): vals.append(s)
    vals.sort()
    if len(vals) < 10: return (float("nan"), float("nan"))
    return (vals[int(0.025*len(vals))], vals[int(0.975*len(vals))])

print("species\tn_hq_genomes\tn_pairs\tsp_2brad\t2brad_lo\t2brad_hi\tsp_16s\t16s_lo\t16s_hi\tn_16s", flush=True)
rows = []
for sd in sorted(glob.glob(f"{GEN}/*/")):
    sp = os.path.basename(sd.rstrip("/"))
    hq = HQ.get(sp, set())
    genomes = [g for g in sorted(glob.glob(f"{sd}/*.fna")) if os.path.basename(g)[:-4] in hq]
    if len(genomes) < 8: print(f"# skip {sp} ({len(genomes)} HQ)", flush=True); continue
    accs = [os.path.basename(g)[:-4] for g in genomes]
    # stage HQ genomes into a temp dir for tag building
    tdir = f"{WORK}/_hq_{sp}"; os.makedirs(tdir, exist_ok=True)
    for g in glob.glob(f"{tdir}/*.fna"): os.remove(g)
    for g in genomes: os.symlink(g, f"{tdir}/{os.path.basename(g)}")
    tags = build_tags(tdir, accs, f"{WORK}/hqtags_{sp}.tsv")
    with ThreadPoolExecutor(max_workers=8) as ex:
        wgs = dict(zip(accs, ex.map(wgs_sketch, genomes)))
        s16 = dict(zip(accs, ex.map(sixteen_s, genomes)))
    common = [a for a in accs if a in tags and wgs.get(a) and s16.get(a)]
    idx = {a: k for k, a in enumerate(common)}
    pairs, dw, d2, d16 = [], [], [], []
    for a, b in itertools.combinations(common, 2):
        Dw = dist(jacc_minhash(wgs[a], wgs[b]))
        if Dw is None: continue
        D2 = dist(jacc(tags[a], tags[b])); D16 = dist(jacc(s16[a], s16[b]))
        if D2 is None or D16 is None: continue
        pairs.append((idx[a], idx[b])); dw.append(Dw); d2.append(D2); d16.append(D16)
    s2 = spearman(dw, d2); s16v = spearman(dw, d16)
    lo2, hi2 = boot_ci(pairs, dw, d2, len(common))
    lo16, hi16 = boot_ci(pairs, dw, d16, len(common))
    row = [sp, len(common), len(dw), round(s2,3), round(lo2,3), round(hi2,3),
           round(s16v,3), round(lo16,3), round(hi16,3), len(common)]
    print("\t".join(str(x) for x in row), flush=True)
    rows.append(row)

with open(OUT_TSV, "w") as f:
    f.write("# 2bRAD-tag AND 16S vs whole-genome between-strain distance concordance (Strain2bScan).\n")
    f.write("# HIGH-QUALITY genomes only: CheckM completeness >=97%, contamination <=5%, assembly level\n")
    f.write("# Complete Genome/Chromosome (genome_qc.tsv). WGS: bottom-3000 21-mer MinHash; 2bRAD:\n")
    f.write("# Strain2bScan build BcgI tags; 16S: longest barrnap gene, 21-mer Jaccard. distance=-ln(2J/(1+J)).\n")
    f.write("# Spearman of pairwise distance vs whole-genome; CI = 2.5/97.5 pctile over 500 genome subsamples (80%).\n")
    f.write("species\tn_hq_genomes\tn_pairs\tsp_2brad\t2brad_lo\t2brad_hi\tsp_16s\t16s_lo\t16s_hi\tn_16s\n")
    for r in rows: f.write("\t".join(str(x) for x in r) + "\n")
print(f"\n# wrote {OUT_TSV}", flush=True)
