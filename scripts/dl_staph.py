#!/usr/bin/env python3
"""Download the S. aureus and S. epidermidis reference panels (for the cross-species
experiment) from NCBI, into the current workdir: staph/<species>/<accession>.fna

Reproducible: uses the pinned accession lists in the repo's data/accessions/
(saureus_panel.txt, sepidermidis_panel.txt). Run from the workdir.
"""
import os, io, zipfile, urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"
PANELS = {
    "Staphylococcus_aureus": "saureus_panel.txt",
    "Staphylococcus_epidermidis": "sepidermidis_panel.txt",
}


def download(acc, out):
    if os.path.exists(out):
        return
    url = f"{API}/{acc}/download?include_annotation_type=GENOME_FASTA"
    z = zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": "s2bs"}), timeout=90).read()))
    fna = next(n for n in z.namelist() if n.endswith("_genomic.fna"))
    with open(out, "wb") as f:
        f.write(z.open(fna).read())


for sp, listfile in PANELS.items():
    d = f"staph/{sp}"
    os.makedirs(d, exist_ok=True)
    accs = [l.strip() for l in open(f"{REPO}/data/accessions/{listfile}") if l.strip()]
    ok = 0
    for a in accs:
        try:
            download(a, f"{d}/{a}.fna")
            ok += 1
        except Exception as e:
            print(f"  FAIL {a}: {e}", flush=True)
    print(f"{sp}: {ok}/{len(accs)} genomes", flush=True)
print("done: staph/{Staphylococcus_aureus,Staphylococcus_epidermidis}/ ready")
