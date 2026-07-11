# figure_raw_data — raw sequencing data supporting the figures

Local store of the **raw sequencing reads** behind the real-data figures (git-ignored — the FASTQ are
large and not committed; the manifests and this README are tracked). Reference genomes are public
assemblies and are given as accession manifests, not stored here. Detailed legends:
`manuscript/figure_legends_detailed.md`.

| Figure / Table | Data type | Contents | Size | Accession |
|---|---|---|---|---|
| **Fig 6** (`Fig06_mock_hostcontam/`) | real native 2bRAD + real shotgun (WMS) | ATCC MSA-1002 mock, 0/90/99 % human; 4 × 2bRAD + 3 × WMS (R1+R2) = 14 FASTQ + `manifest.tsv` | 12 GB | SRA PRJNA1131785 |
| **Fig 7** (`Fig07_saliva_2bRAD/`) | real native BcgI 2bRAD | saliva, 8 subjects × 4 timepoints = 32 R1 FASTQ + `manifest.tsv` (S`<time>`-`<subject>` decode) | 12 GB | SRA PRJNA1131785 |
| **Fig 8** (`Fig08_saliva_shotgun/`) | real shotgun (WMS) | paired saliva WMS, 4 R1 FASTQ (truncated valid-prefix subsamples) + `manifest.tsv` | 22 GB | SRA PRJNA1131785 |
| **Table S3** (`TableS3_clinical_oral/`) | real native BcgI 2bRAD | 4 clinical oral samples (`S_#` cohort), R1 FASTQ + `manifest.tsv` | 0.6 GB | SRA PRJNA1131785 |
| **Fig 2** (`Fig02_16S_panel/`) | reference assemblies (no reads) | 15-species genome QC/accession manifest `genome_accessions_qc.tsv` | — | NCBI / ENA |
| **Fig S1** | real native 2bRAD | MSA-1002 DNA-input titration | — | Figshare 12272360 |

**Simulation-based figures (3, 4, 5, 9C, 10, S2)** use short reads simulated from pinned reference
genomes; they are not stored as reads — regenerate from `scripts/` + `data/accessions/`.

## Re-downloading the reads
All SRA reads are on ENA (reachable even when NCBI is blocked); each sample's `manifest.tsv` gives the
run accession and `fastq_ftp` URL. Example (one saliva sample):
```
curl -s "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR297/.../<SRR>_1.fastq.gz" -o <alias>_1.fastq.gz
```
The Fig 8 WMS files here are valid-prefix subsamples of the full ~14 GB R1 files (sufficient for the
concordance analysis); full-depth files are re-downloadable from the same accessions.

## Provenance
SRA BioProject **PRJNA1131785** (npj Biofilms Microbiomes 2025, s41522-025-00851-2): MSA-1002 mock
(+90/99 % human), HoC saliva (diurnal), oral cancer/ECC; 2bRAD-M / WMS. Figshare **12272360** (2B-RAD-M
unique-tags DB): MSA-1002 native 2bRAD titration.
