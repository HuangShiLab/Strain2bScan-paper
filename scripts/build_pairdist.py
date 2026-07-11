#!/usr/bin/env python3
"""Dump per-pair between-strain distances (whole-genome, 2bRAD, 16S) for each of the 15 species,
high-quality genomes only, for the scatter supplementary. Reuses cached 2bRAD tags (hqtags_*.tsv)
and caches barrnap 16S sequences. Run with PYTHONHASHSEED=0."""
import subprocess, os, glob, math, heapq, itertools
from concurrent.futures import ThreadPoolExecutor

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
BARRNAP = f"{WORK}/barrnap/bin/barrnap"
GEN = f"{WORK}/genomes_16s"
OUTDIR = f"{WORK}/results/pairdist"; os.makedirs(OUTDIR, exist_ok=True)
CACHE16 = f"{WORK}/cache16s"; os.makedirs(CACHE16, exist_ok=True)
K, SKETCH = 21, 3000
COMP = str.maketrans("ACGT", "TGCA")

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

def canon(seq, exact=False):
    seq = seq.upper(); hs = set()
    for i in range(len(seq)-K+1):
        km = seq[i:i+K]
        if not exact and "N" in km: continue
        rc = km.translate(COMP)[::-1]; c = km if km < rc else rc
        hs.add(hash(c) & 0xFFFFFFFFFFFFFFFF)
    return hs

def wgs_sketch(path):
    hs = set()
    for s in read_fna(path): hs |= canon(s)
    return set(heapq.nsmallest(SKETCH, hs))

def sixteen_s_seq(path, sp):
    acc = os.path.basename(path)[:-4]
    cf = f"{CACHE16}/{sp}__{acc}.fa"
    if os.path.exists(cf):
        s = "".join(l.strip() for l in open(cf) if not l.startswith(">"))
        return s or None
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
    open(cf, "w").write(f">{acc}\n{best}\n")
    return best or None

def jmh(a, b):
    u = heapq.nsmallest(SKETCH, a | b); return sum(1 for h in u if h in a and h in b)/len(u) if u else 0.0
def jac(a, b): return len(a & b)/len(a | b) if (a | b) else 0.0
def dist(J):
    if J <= 0: return None
    x = 2*J/(1+J); return -math.log(x) if 0 < x < 1 else 0.0

def load_tags(sp):
    p = f"{WORK}/hqtags_{sp}.tsv"; tags = {}
    if not os.path.exists(p): return tags
    for line in open(p):
        if line.startswith("#"): continue
        q = line.rstrip("\n").split("\t")
        if len(q) >= 2 and q[1]: tags[q[0]] = set(q[1].split(","))
    return tags

for sd in sorted(glob.glob(f"{GEN}/*/")):
    sp = os.path.basename(sd.rstrip("/")); hq = HQ.get(sp, set())
    genomes = [g for g in sorted(glob.glob(f"{sd}/*.fna")) if os.path.basename(g)[:-4] in hq]
    if len(genomes) < 8: continue
    accs = [os.path.basename(g)[:-4] for g in genomes]
    tags = load_tags(sp)
    with ThreadPoolExecutor(max_workers=8) as ex:
        wgs = dict(zip(accs, ex.map(wgs_sketch, genomes)))
        seq16 = dict(zip(accs, ex.map(lambda g: sixteen_s_seq(g, sp), genomes)))
    k16 = {a: (canon(s, exact=True) if s and len(s) > 1000 else None) for a, s in seq16.items()}
    common = [a for a in accs if a in tags and wgs.get(a) and k16.get(a)]
    with open(f"{OUTDIR}/{sp}.tsv", "w") as w:
        w.write("acc_a\tacc_b\twgs_dist\tbrad_dist\ts16_dist\n")
        for a, b in itertools.combinations(common, 2):
            Dw = dist(jmh(wgs[a], wgs[b]))
            if Dw is None: continue
            D2 = dist(jac(tags[a], tags[b])); D16 = dist(jac(k16[a], k16[b]))
            if D2 is None or D16 is None: continue
            w.write(f"{a}\t{b}\t{Dw:.5f}\t{D2:.5f}\t{D16:.5f}\n")
    print(f"{sp}: {len(common)} genomes -> {OUTDIR}/{sp}.tsv", flush=True)
print("PAIRDIST_DONE", flush=True)
