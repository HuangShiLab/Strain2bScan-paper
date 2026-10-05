#!/usr/bin/env python3
"""Build per-species Strain2bScan k-mer databases for the gut panel."""
import os, subprocess, shutil
from pathlib import Path
from collections import defaultdict

FLAT = Path("/Volumes/MoneyCat/Data/GMHI/strain_test/genomes/flat")
BY_SPECIES = Path("/Volumes/MoneyCat/Data/GMHI/strain_test/genomes/by_species")
DBS = Path("/Volumes/MoneyCat/Data/GMHI/strain_test/dbs")
BIN = Path("/Users/macstudio/Downloads/Strain2bScan/target/release/strain2bscan")

BY_SPECIES.mkdir(parents=True, exist_ok=True)
DBS.mkdir(parents=True, exist_ok=True)

# Group genomes by species prefix
species_files = defaultdict(list)
for fna in FLAT.glob("*.fna"):
    prefix = fna.name.split("_", 1)[0]
    species_files[prefix].append(fna)

print(f"Building DBs for {len(species_files)} species")

for species, files in sorted(species_files.items()):
    sp_dir = BY_SPECIES / species
    sp_dir.mkdir(exist_ok=True)
    # Symlink genomes
    for fna in files:
        link = sp_dir / fna.name
        if not link.exists():
            os.symlink(fna, link)
    
    db_path = DBS / f"{species}.tsv"
    if db_path.exists() and db_path.stat().st_size > 0:
        print(f"SKIP {species} (DB exists)")
        continue
    
    print(f"\n=== Building DB for {species} ({len(files)} genomes) ===")
    cmd = [
        str(BIN), "build",
        "--genomes", str(sp_dir),
        "--enzyme", "recommended",
        "--marker-source", "kmer",
        "--kmer-size", "31",
        "--sketch-scale", "30",
        "--out", str(db_path)
    ]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)

print(f"\nAll DBs built in {DBS}")
