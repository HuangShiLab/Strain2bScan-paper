#!/usr/bin/env python3
"""Download a small gut bacterial reference panel for Strain2bScan testing."""
import subprocess, os, json, shutil
from pathlib import Path

OUT = Path("/Volumes/MoneyCat/Data/GMHI/strain_test/genomes")
OUT.mkdir(parents=True, exist_ok=True)
FLAT = OUT / "flat"
FLAT.mkdir(exist_ok=True)

SPECIES = [
    ("Faecalibacterium prausnitzii", 10),
    ("Bacteroides fragilis", 10),
    ("Prevotella copri", 10),
    ("Escherichia coli", 10),
    ("Bifidobacterium longum", 10),
    ("Akkermansia muciniphila", 10),
    ("Roseburia intestinalis", 10),
]

def run(cmd, check=True):
    print(" ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if check and r.returncode != 0:
        print("STDERR:", r.stderr)
        r.check_returncode()
    return r

for sp, limit in SPECIES:
    safe = sp.replace(" ", "_").lower()
    
    # Check if already done
    existing = sorted(FLAT.glob(f"{safe}_*.fna"))
    if len(existing) >= limit:
        print(f"SKIP {sp} ({len(existing)} genomes already in flat dir)")
        continue
    
    # Query genome accessions
    print(f"\n=== Querying {sp} ===")
    r = run([
        "datasets", "summary", "genome", "taxon", sp,
        "--assembly-level", "complete",
        "--assembly-source", "refseq",
        "--as-json-lines"
    ], check=False)
    if r.returncode != 0:
        print(f"  query failed: {r.stderr}")
        continue
    
    accessions = []
    for line in r.stdout.strip().splitlines():
        if not line.strip():
            continue
        try:
            rec = json.loads(line)
            acc = rec.get("accession", "")
            if acc:
                accessions.append(acc)
        except json.JSONDecodeError:
            continue
    
    if len(accessions) > limit:
        accessions = accessions[:limit]
    print(f"  selected {len(accessions)} accessions out of {len(accessions) if len(accessions) >= limit else 'fewer than limit'}")
    
    if not accessions:
        continue
    
    # Download selected accessions
    zip_path = OUT / f"{safe}.zip"
    print(f"=== Downloading {sp} ({len(accessions)} genomes) ===")
    cmd = ["datasets", "download", "genome", "accession"] + accessions + ["--include", "genome", "--filename", str(zip_path)]
    run(cmd)
    
    # Extract
    extract_dir = OUT / safe
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(exist_ok=True)
    run(["unzip", "-q", "-o", str(zip_path), "-d", str(extract_dir)])
    
    # Copy to flat dir
    data_dir = extract_dir / "ncbi_dataset" / "data"
    copied = 0
    for gdir in data_dir.iterdir():
        if not gdir.is_dir():
            continue
        fna = list(gdir.glob("*.fna"))
        if not fna:
            continue
        src = fna[0]
        dst = FLAT / f"{safe}_{src.name}"
        if not dst.exists():
            shutil.copy2(src, dst)
            copied += 1
    print(f"  copied {copied} genomes for {sp}")

# Summary
print(f"\nTotal genomes in {FLAT}: {len(list(FLAT.glob('*.fna')))}")
