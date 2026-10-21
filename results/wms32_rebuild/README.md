# PRJNA1131785 paired saliva WMS rebuild

The complete deterministic dataset contains 32 samples and 64 mate files in canonical FASTQ gzip format:

- HPC directory: `/lustre1/g/aos_shihuang/Strain2b/data/saliva_data/WMS_data/complete_32_fqgz`
- Canonical names: `WMS_<SAMPLE>_1.fq.gz` and `WMS_<SAMPLE>_2.fq.gz`
- `MANIFEST.tsv` maps every mate file to its BioSample-style label and SRA run accession.

Four previously absent/partial libraries were rebuilt from SRA:

| Sample | SRA run |
|---|---|
| S11_1 | SRR29727008 |
| S11_2 | SRR29726982 |
| S11_6 | SRR29726981 |
| S11_11 | SRR29727007 |

`SRR29726985` was also fetched/converted during disambiguation; it corresponds to `S11_14`, whose gz files were already present in the source directory. Conversion logs report exact read/write mate-count equality for every rebuilt run. The canonical directory uses hard links to avoid duplicating 1.3 TB.

## Validation and cleanup

All 64 canonical gzipped mate files passed `gzip -t` (`gzip_validation.tsv`: 64 OK, 0 FAIL). Uncompressed FASTQ, temporary conversion directories, and superseded run logs were removed after conversion. The deleted-path inventory is in `cleanup_intermediates_2026-10-19.tsv`. Raw SRA archives, canonical gzipped FASTQ files, download/conversion scripts, and manifests remain on HPC.

## Full-depth profiling

All 32 paired WMS libraries and all 32 native 2bRAD libraries were profiled against the 19-species oral
BcgI panel with Strain2bScan commit `6ae760c2f3bd2f39efc0e6aa06c53a84366cb16e`. Both runs used 20
panel-specific markers, two detected species markers, support >= 2 and coverage >= 0.01. The WMS call
table is `wms32_strain_calls_all.tsv`; the gate-matched native call table is
`native20_matched/native_strain_calls.tsv`. `per_sample_summary.tsv` and
`full_wms_profile_summary.json` compare the two call sets. This is a modality and sequencing-depth
comparison, not a ground-truth accuracy validation.
