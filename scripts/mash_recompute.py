#!/usr/bin/env python3
"""Recompute the 2bRAD-vs-whole-genome strain-distance concordance with Strain2bScan.

For each species: whole-genome Mash distance (pure-Python bottom-k MinHash, k=21) vs 2bRAD-tag
distance (Strain2bScan `build` per-genome tag sets; both the 14-enzyme default and BcgI-alone).
Both distances use the Mash transform D(J) = -ln(2J/(1+J)); we correlate the pairwise vectors
(Pearson + Spearman). 16S is deferred (needs an rRNA finder; not installed). Run with
PYTHONHASHSEED=0 for reproducible sketches.
"""
import subprocess, os, glob, math, heapq, itertools

BIN = os.path.expanduser("~/Downloads/Strain2bScan/target/release/strain2bscan")
SC = "/private/tmp/claude-501/-Users-shihuang-Downloads-StrainScan-rust-claude/9f8a698a-0227-4dd2-8236-2536be092ed9/scratchpad"
FULLSET = "CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"
SPECIES = [
    ("Cutibacterium_acnes",        f"{SC}/rerun/acnes/genomes",                 0.985),
    ("Staphylococcus_aureus",      f"{SC}/rerun/staph/Staphylococcus_aureus",   0.970),
    ("Staphylococcus_epidermidis", f"{SC}/rerun/staph/Staphylococcus_epidermidis", 0.947),
]
K, SKETCH = 21, 3000
COMP = str.maketrans("ACGT", "TGCA")


def read_fna(path):
    seqs, cur = [], []
    for line in open(path):
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def wgs_sketch(path):
    hs = set()
    for seq in read_fna(path):
        seq = seq.upper()
        for i in range(len(seq) - K + 1):
            kmer = seq[i:i + K]
            rc = kmer.translate(COMP)[::-1]
            c = kmer if kmer < rc else rc
            hs.add(hash(c) & 0xFFFFFFFFFFFFFFFF)
    return set(heapq.nsmallest(SKETCH, hs))


def jacc_minhash(a, b):
    union = heapq.nsmallest(SKETCH, a | b)
    inter = sum(1 for h in union if h in a and h in b)
    return inter / len(union) if union else 0.0


def dist(J):
    if J <= 0:
        return 10.0
    x = 2 * J / (1 + J)
    return -math.log(x) if 0 < x < 1 else 0.0


def build_tags(gdir, enzyme, out):
    subprocess.run([BIN, "build", "--genomes", gdir, "--enzyme", enzyme, "--out", out],
                   capture_output=True, text=True, env={**os.environ, "STRAIN2BSCAN_THREADS": "16"})
    tags = {}
    for line in open(out):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[1]:
            tags[p[0]] = set(p[1].split(","))
    return tags


def pearson(x, y):
    n = len(x); mx = sum(x) / n; my = sum(y) / n
    sx = sum((a - mx) ** 2 for a in x); sy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sx * sy) if sx > 0 and sy > 0 else 0.0


def spearman(x, y):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v); i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for t in range(i, j + 1):
                r[order[t]] = avg
            i = j + 1
        return r
    return pearson(rank(x), rank(y))


OUT_TSV = os.path.expanduser("~/Downloads/Strain2bScan-paper/results/mash_2brad_vs_wgs_recomputed.tsv")
HDR = ["species", "enzyme_set", "n_genomes", "n_pairs_used", "n_pairs_undefined", "pearson", "spearman", "PPTX_2brad"]
print("\t".join(HDR))
allrows = []
for sp, gdir, pptx in SPECIES:
    genomes = sorted(glob.glob(f"{gdir}/*.fna"))
    accs = [os.path.basename(g)[:-4] for g in genomes]
    tags_by_set = {"14enz": build_tags(gdir, FULLSET, f"{SC}/2bfull_{sp}.tsv"),
                   "BcgI": build_tags(gdir, "BcgI", f"{SC}/2bbcgi_{sp}.tsv")}
    sk = {acc: wgs_sketch(g) for acc, g in zip(accs, genomes)}
    for label, tags in tags_by_set.items():
        common = [a for a in accs if a in tags and a in sk]
        dw, d2, undef = [], [], 0
        for a, b in itertools.combinations(common, 2):
            Jw = jacc_minhash(sk[a], sk[b])
            Ta, Tb = tags[a], tags[b]
            J2 = len(Ta & Tb) / len(Ta | Tb) if (Ta | Tb) else 0.0
            if Jw <= 0 or J2 <= 0:          # distance undefined in one space (beyond resolution)
                undef += 1
                continue
            dw.append(dist(Jw)); d2.append(dist(J2))
        row = [sp, label, len(common), len(dw), undef,
               round(pearson(dw, d2), 3), round(spearman(dw, d2), 3), pptx]
        print("\t".join(str(x) for x in row), flush=True)
        allrows.append(row)
os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
with open(OUT_TSV, "w") as f:
    f.write("# 2bRAD-tag vs whole-genome between-strain distance concordance, recomputed with Strain2bScan.\n")
    f.write(f"# WGS: bottom-{SKETCH} canonical {K}-mer MinHash Mash distance. 2bRAD: per-genome Strain2bScan\n")
    f.write("# `build` tag sets (14-enzyme default and BcgI-alone). distance = -ln(2J/(1+J)); pairs with\n")
    f.write("# zero overlap in either space excluded (n_pairs_undefined). 16S deferred (needs an rRNA finder,\n")
    f.write("# e.g. barrnap; not installed here). PPTX_2brad = prior Strain2bfunc value (slide 23) for reference.\n")
    f.write("\t".join(HDR) + "\n")
    for r in allrows:
        f.write("\t".join(str(x) for x in r) + "\n")
print(f"\n# wrote {OUT_TSV}")
