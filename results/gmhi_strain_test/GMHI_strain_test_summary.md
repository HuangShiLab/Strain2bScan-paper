# GMHI PRJEB12449 strain-level profiling pilot

Date: 2026-09-25

## Objective

Test whether Strain2bScan can perform strain-level profiling on shotgun
metagenome (WMS) data from the GMHI cohort, using a small gut bacterial
reference panel.

## Reference panel

| Species | N complete RefSeq genomes downloaded |
|---|---|
| Faecalibacterium prausnitzii | 25 |
| Bacteroides fragilis | 42 |
| Prevotella copri | 24 |
| Escherichia coli | 10 |
| Bifidobacterium longum | 10 |
| Akkermansia muciniphila | 10 |
| Roseburia intestinalis | 2 |
| **Total** | **123** |

Marker space: k-mer sketch (k=31, scale=30) using `strain2bscan build`.
Per-species DBs were placed in `dbs/` and profiled with `multi-profile`.

## Samples

PRJEB12449 (Zeller et al. CRC/Overweight/Obesity/Healthy gut microbiome) is the
only fully downloaded BioProject (1,323 FASTQ files, ~1.0 TB).
Two Healthy individuals were selected and all technical replicates merged:

| Sample | Phenotype | Replicates | Files | Merged size |
|---|---|---|---|---|
| SAMEA3879525 | Healthy | 8 runs | 12 files (4 SE + 4 PE pairs) | 3.8 GB |
| SAMEA3879530 | Healthy | 8 runs | 12 files (4 SE + 4 PE pairs) | 4.5 GB |

Note: the project contains a mixture of single-end (`ERR*.fastq.gz`) and
paired-end (`ERR*_1.fastq.gz` / `ERR*_2.fastq.gz`) runs. Both mates of
paired-end runs were concatenated before profiling.

## Performance

| Sample | Merged reads (compressed) | Real time | Peak memory | Distinct sample k-mers |
|---|---|---|---|---|
| SAMEA3879525 | 3.8 GB | 48.6 s | 923 MB | 11,588,883 |
| SAMEA3879530 | 4.5 GB | 51.1 s | 1,345 MB | 15,409,820 |

Throughput: ~90 MB/s compressed input, or ~250,000–300,000 reads/s
(assuming ~100 bp reads).

## Strain calls

| Sample | Species | Strain cluster | Coverage | Support | Depth | Global abundance |
|---|---|---|---|---|---|---|
| SAMEA3879525 | Roseburia intestinalis | GCF_025151715.1_ASM2515171v1 | 0.57 | 100 | 1.28x | 0.872 |
| SAMEA3879525 | Bifidobacterium longum | GCF_005406285.1_ASM540628v1 | 0.12 | 614 | 0.19x | 0.128 |
| SAMEA3879530 | Roseburia intestinalis | GCF_025151715.1_ASM2515171v1 | 0.61 | 107 | 6.36x | 1.000 |

Both individuals were assigned the **same Roseburia intestinalis strain**
(GCF_025151715.1). This may reflect a common reference strain in the panel
(only 2 Roseburia genomes were available) rather than true biological
identity.

## Detected-but-not-resolved species

Species with markers seen in the sample but no strain cluster passing the
support gate:

| Sample | Akkermansia | Bacteroides | Bifidobacterium | Escherichia | Faecalibacterium | Prevotella | Roseburia |
|---|---|---|---|---|---|---|---|
| SAMEA3879525 | 5.54x | 6.16x | called | 0.01x | 0.64x | 12.31x | called |
| SAMEA3879530 | 0.01x | 4.45x | 0.03x | 0.03x | 3.70x | 0.03x | called |

Faecalibacterium prausnitzii, despite being a dominant gut species, was not
resolved to a strain in either sample. Its DB contains 25 genomes; at the
observed depths (0.64x and 3.70x) there is insufficient unique-marker coverage
to separate them.

## Coverage of the panel

Only **0.2–0.8%** of sample k-mer observations were assigned to the 7-species
panel. The remaining 99+% is host DNA, other microbes, sequencing error, or
species absent from the panel. This is expected for a tiny reference set.

## Conclusions and next steps

1. **Feasibility**: Strain2bScan k-mer mode works on WGS metagenome data.
   Runtime is fast (~50 s for 4 GB compressed WMS) and memory modest
   (~1 GB) for a 123-genome panel.

2. **Panel size is the bottleneck**: A 7-species panel misses >99% of reads.
   For GMHI application, a panel of hundreds of common gut species (e.g.
   GTDB representative genomes or the GMHI species set) is needed.

3. **Strain resolution requires depth**: Even within the 7 species, many were
   detected but not resolved. Faecalibacterium prausnitzii at 3.7x depth did
   not produce a strain call with 25 reference genomes, suggesting that deeper
   sequencing or a sparser representative panel may be needed.

4. **Recommended follow-up**:
   - Build a larger gut panel (50–200 species) using GTDB representatives.
   - Profile more PRJEB12449 samples across phenotypes (Healthy, CRC,
     Overweight, Obesity) to see if strain-level features improve disease
     classification.
   - Compare Strain2bScan strain calls to species-level abundance from a
     conventional profiler (e.g. MetaPhlAn, Kraken2/Bracken) for accuracy
     benchmarking.

## Files

| Path | Description |
|---|---|
| `genomes/flat/` | 123 reference genomes |
| `genomes/by_species/` | Genomes grouped by species prefix |
| `dbs/` | Per-species Strain2bScan k-mer DBs |
| `merged_reads/` | Merged FASTQ per individual |
| `preds/` | Strain-level predictions and stdout logs |
| `download_panel.py` | Genome download script |
| `build_dbs.py` | DB build script |
| `merge_samples.py` | Sample merge script |
