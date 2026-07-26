#!/usr/bin/env python3
"""Fair-play control for inStrain: dereplicate the 120-genome reference at 98% ANI.

Why: inStrain's documentation requires a DEREPLICATED genome set (they recommend dRep at ~98%
ANI). Our DB deliberately holds 6 conspecifics per species at 95-99.9% ANI, so bowtie2
multi-maps reads across near-identical references and inflates breadth on decoys -- that is
what produced inStrain's 48 false positives at 0% host. Scoring inStrain on a reference that
violates its stated assumption is not a fair comparison, so we also run it the way its authors
intend and report both.

Representative choice mirrors dRep's quality-based pick: ATCC reference genomes are complete
closed assemblies and outscore the draft spike-ins, so ATCC wins its cluster -- which is also
the generous choice for inStrain (it gets the exact target as its mapping reference).

Scoring stays in cluster space, identical to Strain2bScan/StrainScan: a 98% cluster is a true
positive when it contains an ATCC target.
"""
import subprocess, glob, os, sys, json
from collections import defaultdict

RAW = "/Users/macstudio/Downloads/Strain2bScan-raw-data"
SP = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
OUT = f"{SP}/instrain_derep"
os.makedirs(OUT, exist_ok=True)
ANI_CUT = 98.0

gens = sorted(glob.glob(f"{RAW}/MSA1002_all120/*.f*a"))
print(f"input genomes: {len(gens)}", file=sys.stderr)

# ---- pairwise ANI (skani triangle) -------------------------------------------------
mat = f"{OUT}/ani120.tsv"
if not os.path.exists(mat):
    print("running skani triangle...", file=sys.stderr)
    with open(mat, "w") as fh:
        subprocess.run(["skani", "triangle", "-t", "8", "--full-matrix", "-m", "200"] + gens,
                       stdout=fh, stderr=subprocess.DEVNULL, timeout=3600)

lines = [l.rstrip("\n") for l in open(mat) if l.strip()]
n = int(lines[0])
names, rows = [], []
for l in lines[1:]:
    p = l.split("\t")
    names.append(os.path.basename(p[0]))
    rows.append([float(x) if x else 0.0 for x in p[1:]])

# ---- greedy single-linkage clustering at 98% ANI ------------------------------------
parent = list(range(len(names)))
def find(a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]; a = parent[a]
    return a
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra

for i in range(len(names)):
    for j in range(i + 1, len(names)):
        if rows[i][j] >= ANI_CUT:
            union(i, j)

clust = defaultdict(list)
for i, nm in enumerate(names):
    clust[find(i)].append(nm)
print(f"98% ANI clusters: {len(clust)} (from {len(names)} genomes)", file=sys.stderr)

# ---- pick representative: ATCC if present (dRep would pick it -- complete assembly) --
def stem(f):
    return f.rsplit(".", 1)[0]

reps, members = [], {}
for k, ms in clust.items():
    atcc = [m for m in ms if "ATCC_" in m]
    rep = atcc[0] if atcc else sorted(ms)[0]
    reps.append(rep)
    members[stem(rep)] = [stem(m) for m in ms]

n_atcc_rep = sum(1 for r in reps if "ATCC_" in r)
print(f"representatives: {len(reps)} ({n_atcc_rep} are ATCC targets)", file=sys.stderr)

# ---- write ref98.fasta + .stb -------------------------------------------------------
fa, stb = f"{OUT}/ref98.fasta", f"{OUT}/ref98.stb"
byname = {os.path.basename(g): g for g in gens}
with open(fa, "w") as out, open(stb, "w") as sf:
    for rep in reps:
        gid = stem(rep)
        for ln in open(byname[rep]):
            if ln.startswith(">"):
                scaf = f"{gid}__{ln[1:].split()[0]}"
                out.write(f">{scaf}\n")
                sf.write(f"{scaf}\t{gid}\n")
            else:
                out.write(ln)
json.dump(members, open(f"{OUT}/ref98_members.json", "w"))
print(f"wrote {fa}", file=sys.stderr)
