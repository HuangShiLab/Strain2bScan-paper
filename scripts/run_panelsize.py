#!/usr/bin/env python3
"""Test the prediction: P. copri's precision problem (recombination-driven false uniqueness)
should IMPROVE as the reference panel grows -- a larger panel is more likely to already contain,
in some OTHER genome, a given recombinant segment, reclassifying it shared not falsely 'unique'.

Design: strictly-nested panels 40 subset of 80 subset of 112 (StrainScan's whole P. copri
panel), the SAME 5 simulated samples across all three sizes. Only Strain2bScan is exercised
(this is a Strain2bScan-internal prediction; no StrainScan / no cross-tool comparison), so
genome source doesn't need to match StrainScan's DB.

Genomes are fetched from EBI/ENA (NCBI was IP-blocked during this run). GCF (RefSeq) accessions
map to GCA (GenBank) by identical numeric part; ENA's unversioned /fasta/<GCA> auto-resolves to
the current version. Content saved under the ORIGINAL GCF stem so truth/DB labels stay
consistent with data/accessions/prevotella_copri_40x.txt.

The 5 samples reproduce the original species_expansion run's RNG sequence (seed=102): the 40
baseline = random.Random(102).sample(sorted(all-112), 40), then truth_pool = sample(that, 5),
then per-sample draws -- so the 40-genome precision here should track the published 0.231
(ENA/GCA vs NCBI/GCF assemblies are near-identical for these prokaryotic WGS genomes).
"""
import io, gzip, os, re, random, subprocess, time, urllib.request

BIN = os.environ.get("STRAIN2BSCAN_BIN", "strain2bscan")
SET14 = "CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"
WORK = "panelsize"
ENA = "https://www.ebi.ac.uk/ena/browser/api/fasta"
os.makedirs(WORK, exist_ok=True)


def ena_download(gcf, out, retries=4):
    """Download the assembly FASTA from ENA for GCF's paired GCA, save under the GCF stem."""
    if os.path.exists(out):
        return True
    gca = "GCA_" + gcf[len("GCF_"):].split(".")[0]
    url = f"{ENA}/{gca}?download=true&gzip=true"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "s2bs-panelsize"})
            raw = urllib.request.urlopen(req, timeout=120).read()
            data = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
            if not data.startswith(b">"):
                raise ValueError("not FASTA")
            with open(out, "wb") as o:
                o.write(data)
            return True
        except Exception as e:
            if attempt == retries - 1:
                print(f"    FAIL {gcf} ({gca}): {e}", flush=True)
            time.sleep(2)
    return False


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


def timed_run(cmd, env=None):
    r = subprocess.run(["/usr/bin/time", "-l"] + cmd, capture_output=True, text=True, env=env)
    t = re.search(r"([0-9.]+) real", r.stderr); m = re.search(r"(\d+)\s+maximum resident", r.stderr)
    return r.stdout, (float(t.group(1)) if t else 0), (int(m.group(1)) / 1048576 if m else 0)


# --- 1) full 112-accession list (from StrainScan's own P. copri DB), sorted (as originally) ---
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_ACC = os.path.join(_REPO, "data", "accessions", "prevotella_copri_112.txt")
ALL112 = sorted(l.strip() for l in open(_ACC) if l.strip() and not l.startswith("#"))
print(f"full panel: {len(ALL112)} accessions", flush=True)

gdir = f"{WORK}/genomes"; os.makedirs(gdir, exist_ok=True)
print("downloading all 112 genomes from ENA ...", flush=True)
have = []
for i, gcf in enumerate(ALL112):
    if ena_download(gcf, f"{gdir}/{gcf}.fna"):
        have.append(gcf)
    if (i + 1) % 20 == 0:
        print(f"  {i+1}/{len(ALL112)} attempted, {len(have)} ok", flush=True)
print(f"  {len(have)}/{len(ALL112)} genomes downloaded", flush=True)

# --- 2) reproduce the original species_expansion RNG sequence (seed=102) ---
rng = random.Random(102)
subset40 = rng.sample(ALL112, 40)                 # the original 40-genome baseline (exact order)
missing = [a for a in subset40 if a not in set(have)]
if missing:
    print(f"  WARNING: {len(missing)} baseline genomes missing from ENA: {missing}", flush=True)
base40 = [a for a in subset40 if a in set(have)]   # keep order, drop any un-downloadable
truth_pool = rng.sample(base40, 5)
print(f"truth pool (5): {truth_pool}", flush=True)

sdir, tdir = f"{WORK}/samples", f"{WORK}/truth"
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
print("generated 5 samples", flush=True)

# --- 3) strictly-nested panels: 40 subset of 80 subset of full ---
extra = [a for a in ALL112 if a in set(have) and a not in set(base40)]   # deterministic (sorted) order
panels = [("40", base40),
          ("80", base40 + extra[:40]),
          (str(len(base40) + len(extra)), base40 + extra)]

env16 = {**os.environ, "STRAIN2BSCAN_THREADS": "16"}
results = []
for label, members in panels:
    pdir = f"{WORK}/panel_{label}"; os.makedirs(pdir, exist_ok=True)
    for f in os.listdir(pdir):
        os.remove(os.path.join(pdir, f))
    for a in members:
        os.symlink(os.path.abspath(f"{gdir}/{a}.fna"), f"{pdir}/{a}.fna")
    print(f"\n== panel {label}: {len(members)} genomes ==", flush=True)

    db_path = f"{WORK}/db_{label}.tsv"
    bout, bt, brss = timed_run([BIN, "cluster", "--genomes", pdir, "--enzyme", SET14,
                                "--out", db_path, "--similarity", "0.95"], env16)
    ncl = int(re.search(r"into (\d+) cluster", bout).group(1)) if "into" in bout else 0
    print(f"  build: {bt:.2f}s {brss:.0f}MB, {ncl} clusters", flush=True)
    g2c = {}
    for l in open(db_path.replace(".tsv", ".members.tsv")):
        if not l.startswith("#"):
            g, c = l.strip().split("\t"); g2c[g] = c

    TP = FP = FN = 0
    for i in range(1, 6):
        pred = f"{WORK}/s_{label}_{i}.pred"
        timed_run([BIN, "profile", "--db", db_path, "--reads", f"{sdir}/sample{i}.fq", "--out", pred], env16)
        ab = {}
        for l in open(f"{tdir}/sample{i}.truth.tsv"):
            if l.startswith("#"):
                continue
            g, a = l.strip().split("\t"); c = g2c.get(g, g); ab[c] = ab.get(c, 0) + float(a)
        truth_f = f"{WORK}/s_{label}_{i}.truth"
        open(truth_f, "w").write("".join(f"{c}\t{a}\n" for c, a in ab.items()))
        ev = subprocess.run([BIN, "evaluate", "--pred", pred, "--truth", truth_f, "--present", "0.01"],
                            capture_output=True, text=True).stdout
        TP += int(re.search(r"TP=(\d+)", ev).group(1)); FP += int(re.search(r"FP=(\d+)", ev).group(1)); FN += int(re.search(r"FN=(\d+)", ev).group(1))
    P = TP / (TP + FP) if TP + FP else 0; R = TP / (TP + FN) if TP + FN else 0
    print(f"  precision={P:.3f} recall={R:.3f}  (TP={TP} FP={FP} FN={FN})", flush=True)
    results.append(dict(panel_size=len(members), n_clusters=ncl,
                        precision=round(P, 3), recall=round(R, 3), TP=TP, FP=FP, FN=FN))

with open(f"{WORK}/panelsize_results.tsv", "w") as o:
    cols = list(results[0].keys())
    o.write("\t".join(cols) + "\n")
    for r in results:
        o.write("\t".join(str(r[c]) for c in cols) + "\n")
print(f"\ndone -> {WORK}/panelsize_results.tsv", flush=True)
for r in results:
    print(" ", r, flush=True)
