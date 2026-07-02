#!/usr/bin/env python3
"""Download the 55-species x 4-strain multi-species panel (218 genomes) by PINNED accession,
into the current workdir: multispecies/genomes/<species>/<accession>.fna

Reproducible: reads the fixed (species, accession) pairs from the repo's
data/accessions/multispecies_55x4.tsv (originally resolved via the NCBI Datasets taxon API,
frozen here so re-runs always fetch the identical genomes). Run from the workdir.
"""
import io, os, zipfile, urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
OUT = "multispecies/genomes"


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "strain2bscan-bench"})
    return urllib.request.urlopen(req, timeout=90).read()


def download(acc, out):
    if os.path.exists(out):
        return True
    url = f"{BASE}/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA"
    z = zipfile.ZipFile(io.BytesIO(get(url)))
    fna = [n for n in z.namelist() if n.endswith("_genomic.fna")]
    if not fna:
        return False
    with z.open(fna[0]) as fh, open(out, "wb") as o:
        o.write(fh.read())
    return True


pairs = []
for l in open(f"{REPO}/data/accessions/multispecies_55x4.tsv"):
    if l.startswith("#") or not l.strip():
        continue
    sp, acc = l.rstrip("\n").split("\t")
    pairs.append((sp.replace(" ", "_"), acc))

tot = 0
by_species = {}
for sp, acc in pairs:
    by_species.setdefault(sp, []).append(acc)

for sp, accs in by_species.items():
    d = os.path.join(OUT, sp)
    os.makedirs(d, exist_ok=True)
    ok = 0
    for a in accs:
        try:
            if download(a, os.path.join(d, f"{a}.fna")):
                ok += 1
        except Exception as e:
            print(f"  {a} FAIL: {e}", flush=True)
    tot += ok
    print(f"{sp}: {ok}/{len(accs)} genomes", flush=True)

print(f"TOTAL genomes: {tot} across {len(by_species)} species")
