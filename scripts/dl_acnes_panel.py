#!/usr/bin/env python3
"""Set up the C. acnes data for the reference-quality experiment, in the current workdir:
  acnes/genomes/   64-genome panel (14 truth + 50 background) from NCBI by accession
  acnes/reads/     5 mock samples (sample{1..5}.fq) from MockMetagenomes4Benchmark
  acnes/truth/     per-sample ground truth (copied from the repo's data/truth/)
  multispecies/genomes/Escherichia_coli/  one E. coli genome (degradation contaminant)

Run from the workdir, e.g.  cd work && python3 ../scripts/dl_acnes_panel.py
Idempotent: skips files already present.
"""
import os, io, gzip, glob, shutil, zipfile, urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"
MOCK = ("https://raw.githubusercontent.com/HuangShiLab/MockMetagenomes4Benchmark/"
        "main/Cutibacterium_acnes_0.01/samples")
ECOLI = "GCF_000005845.2"  # E. coli K-12 MG1655 (stable complete reference) — contaminant


def fetch(url, timeout=300):
    req = urllib.request.Request(url, headers={"User-Agent": "strain2bscan-bench"})
    return urllib.request.urlopen(req, timeout=timeout).read()


def download_genome(acc, out):
    if os.path.exists(out):
        return
    z = zipfile.ZipFile(io.BytesIO(fetch(f"{API}/{acc}/download?include_annotation_type=GENOME_FASTA")))
    fna = next(n for n in z.namelist() if n.endswith("_genomic.fna"))
    with open(out, "wb") as f:
        f.write(z.open(fna).read())


os.makedirs("acnes/genomes", exist_ok=True)
os.makedirs("acnes/reads", exist_ok=True)
os.makedirs("acnes/truth", exist_ok=True)
os.makedirs("multispecies/genomes/Escherichia_coli", exist_ok=True)

# 1) 64-genome panel
accs = [l.strip() for l in open(f"{REPO}/data/accessions/cae_panel_64.txt") if l.strip()]
for i, acc in enumerate(accs, 1):
    download_genome(acc, f"acnes/genomes/{acc}.fna")
    if i % 10 == 0:
        print(f"  genomes {i}/{len(accs)}", flush=True)
print(f"panel: {len(glob.glob('acnes/genomes/*.fna'))} genomes")

# 2) E. coli contaminant
download_genome(ECOLI, f"multispecies/genomes/Escherichia_coli/{ECOLI}.fna")

# 3) 5 mock samples (merge left+right paired reads)
for i in range(1, 6):
    out = f"acnes/reads/sample{i}.fq"
    if os.path.exists(out):
        continue
    with open(out, "wb") as o:
        for side in ("left", "right"):
            o.write(gzip.decompress(fetch(f"{MOCK}/sample{i}.{side}.fq.gz")))
    print(f"  reads sample{i}", flush=True)

# 4) ground truth (committed in the repo)
for f in glob.glob(f"{REPO}/data/truth/sample*.truth.tsv"):
    shutil.copy(f, f"acnes/truth/{os.path.basename(f)}")

print("acnes data ready in", os.getcwd())
