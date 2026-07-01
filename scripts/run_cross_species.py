#!/usr/bin/env python3
"""Cross-species generalization: same 14-enzyme, coverage-gated pipeline on 3 species.
C. acnes uses the real MockMetagenomes samples; S. aureus/S. epidermidis use simulated
strain-mixture mocks (sim_strain_mock.py) from real NCBI panels. -> cross_species.tsv"""
import subprocess, re, os, sys, glob, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim_strain_mock
BIN = os.environ.get("STRAIN2BSCAN_BIN","strain2bscan")
SET = "CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"
env = {**os.environ, "STRAIN2BSCAN_THREADS": "16"}
def timed(cmd):
    r = subprocess.run(["/usr/bin/time", "-l", BIN] + cmd, capture_output=True, text=True, env=env)
    t = re.search(r'([0-9.]+) real', r.stderr)
    return r.stdout, (float(t.group(1)) if t else 0)
def run(name, gdir, reads, truth, seed):
    if reads is None:  # simulate mocks
        mock = f"{gdir}_mock"; sim_strain_mock.main(gdir, mock, 5, seed)
        reads, truth = f"{mock}/reads", f"{mock}/truth"
    ng = len(glob.glob(gdir + "/*.fna"))
    db = f"xspec/{name}.db.tsv"; os.makedirs("xspec", exist_ok=True)
    so, bt = timed(["cluster", "--genomes", gdir, "--enzyme", SET, "--out", db, "--similarity", "0.95"])
    ss = int(re.search(r'strain_specific=(\d+)', so).group(1)); ncl = int(re.search(r'into (\d+) cluster', so).group(1))
    g2c = {}
    for l in open(db.replace(".tsv", ".members.tsv")):
        if not l.startswith("#"): g, c = l.strip().split("\t"); g2c[g] = c
    TP = FP = FN = 0; bcs = []; pts = []
    for i in range(1, 6):
        pred = f"xspec/{name}_s{i}.pred"
        _, pt = timed(["profile", "--db", db, "--reads", f"{reads}/sample{i}.fq", "--out", pred]); pts.append(pt)
        ab = {}
        for l in open(f"{truth}/sample{i}.truth.tsv"):
            if l.startswith("#"): continue
            g, a = l.strip().split("\t"); c = g2c.get(g, g); ab[c] = ab.get(c, 0) + float(a)
        tf = f"xspec/{name}_s{i}.truth"; open(tf, "w").write("".join(f"{c}\t{a}\n" for c, a in ab.items()))
        ev = subprocess.run([BIN, "evaluate", "--pred", pred, "--truth", tf, "--present", "0.01"], capture_output=True, text=True).stdout
        TP += int(re.search(r'TP=(\d+)', ev).group(1)); FP += int(re.search(r'FP=(\d+)', ev).group(1)); FN += int(re.search(r'FN=(\d+)', ev).group(1))
        bcs.append(float(re.search(r'Bray-Curtis=([0-9.]+)', ev).group(1)))
    P = TP/(TP+FP) if TP+FP else 0; R = TP/(TP+FN) if TP+FN else 0
    return dict(species=name, n_genomes=ng, n_clusters=ncl, cluster_frac=round(ncl/ng, 2),
                strain_specific=ss, build_s=round(bt, 2), profile_s=round(statistics.mean(pts), 2),
                precision=round(P, 3), recall=round(R, 3), bray_curtis=round(statistics.mean(bcs), 3))
rows = [
    run("Cutibacterium_acnes", "acnes/genomes", "acnes/reads", "acnes/truth", 0),
    run("Staphylococcus_aureus", "staph/Staphylococcus_aureus", None, None, 1),
    run("Staphylococcus_epidermidis", "staph/Staphylococcus_epidermidis", None, None, 2),
]
cols = ["species","n_genomes","n_clusters","cluster_frac","strain_specific","build_s","profile_s","precision","recall","bray_curtis"]
os.makedirs("xspec", exist_ok=True)
with open("xspec/cross_species.tsv", "w") as o:
    o.write("\t".join(cols)+"\n")
    for r in rows: o.write("\t".join(str(r[c]) for c in cols)+"\n")
for r in rows: print("  "+" ".join(f"{c}={r[c]}" for c in ["species","n_genomes","n_clusters","strain_specific","precision","recall","bray_curtis"]), flush=True)
print("done -> xspec/cross_species.tsv")
