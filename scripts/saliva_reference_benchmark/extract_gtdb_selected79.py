#!/usr/bin/env python3
import gzip,re,sys
from pathlib import Path
from collections import defaultdict
base="/lustre1/g/aos_shihuang"
map_file=Path(base)/"selected79_gtdb_species_map.tsv"
outroot=Path(base)/"projects/strain2bscan_reference_benchmark_20261120/gtdb_pseudo"
map_genome={}; map_species={}
with map_file.open() as f:
    for line in f:
        genome,species=line.rstrip("\n").split("\t")
        map_genome[genome]=species
        map_species[species]=True

def safe_name(sp):
    stem=sp.removeprefix("s__")
    stem=re.sub(r"[^A-Za-z0-9_.]+","_",stem).strip("_")
    return "GTDB__"+stem

used=defaultdict(int)
safe_by_species={}
for sp in map_species:
    stem=safe_name(sp)
    if stem in used:
        used[stem]+=1
        stem=f"{stem}_{used[stem]}"
    else:
        used[stem]=1
    safe_by_species[sp]=stem
    (outroot/stem).mkdir(parents=True,exist_ok=True)

seen_headers=0; written_records=0; written_genomes=set(); bytes_seen=0
fasta=Path(base)/"databases/GTDB/GTDBr226/2bradm_gtdb_db/db/BcgI.species.fa.gz"
with gzip.open(fasta,"rt") as src, open(Path(base)/"projects/strain2bscan_reference_benchmark_20261120/results/gtdb_extraction.log","w") as log:
    header=None; seq=[]
    for line in src:
        bytes_seen += len(line)
        if line.startswith(">"):
            if header is not None:
                accession=header[1:].split("|",1)[0]
                if accession in map_genome:
                    sp=map_genome[accession]
                    outdir=outroot/safe_by_species[sp]
                    out=outdir/f"{accession}.fna"
                    with out.open("a") as outfh:
                        outfh.write(header+"\n"+"".join(seq)+"\n")
                    written_records += 1
                    written_genomes.add(accession)
            header=line.rstrip("\n"); seq=[]
            seen_headers += 1
            if seen_headers % 1000000 == 0:
                print(f"headers={seen_headers} written_records={written_records} genomes={len(written_genomes)}",flush=True)
        else:
            seq.append(line)
    if header is not None:
        accession=header[1:].split("|",1)[0]
        if accession in map_genome:
            sp=map_genome[accession]
            outdir=outroot/safe_by_species[sp]
            out=outdir/f"{accession}.fna"
            with out.open("a") as outfh:
                outfh.write(header+"\n"+"".join(seq)+"\n")
            written_records += 1
            written_genomes.add(accession)
    print(f"done headers={seen_headers} written_records={written_records} genomes={len(written_genomes)} bytes={bytes_seen}",flush=True)
    log.write(f"headers={seen_headers}\nrecords={written_records}\ngenomes={len(written_genomes)}\nbytes={bytes_seen}\n")

manifest=Path(base)/"projects/strain2bscan_reference_benchmark_20261120/results/gtdb_selected79_panel_manifest.tsv"
with manifest.open("w") as w:
    w.write("safe_species\tgtdb_species\tgenome_dir\tn_genome_files\n")
    for sp,stem in sorted(safe_by_species.items()):
        n=len(list((outroot/stem).glob("*.fna")))
        w.write(f"{stem}\t{sp}\t{outroot/stem}\t{n}\n")
print(f"manifest={manifest}",flush=True)
