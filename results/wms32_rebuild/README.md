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
