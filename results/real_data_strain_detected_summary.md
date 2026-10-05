# Strain-level + detected-not-resolved depth ML comparison

Date: 2026-09-24

## What was tested

`Strain2bScan` strain-level calls (sample_fraction per strain cluster) were compared
against species-level abundance tables, with an additional feature set that also
includes the **depth of species that were detected but not resolved to a strain**.

Detected-not-resolved species come from `multi-profile` stdout and fall into two
categories:

1. `[detected, not strain-resolvable]` – species present in the panel but without a
   strain-resolved database, or with too few markers.
2. `[strain-resolved, no cluster above threshold]` – species with a strain DB, but
   no cluster passed the support/coverage gates.

These species are **not written to `.pred.tsv`**; they only appear in stdout.
The helper `scripts/parse_detected_depth.py` parses a saved stdout log into a
per-sample, per-species depth matrix.

## Datasets

| Dataset | Samples | Class balance | Source |
|---------|---------|---------------|--------|
| ECC saliva 2bRAD | 38 | 19 caries / 19 healthy | `/lustre1/g/aos_shihuang/Strain2b/data/saliva_data/ECC_saliva/2b/data/Clean_data` |
| Lim_ORPI denture 2bRAD-M | 97 | 42 unclean / 55 clean | `/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/lim_orpi_denture/fastq` |

## Method

- Strain features: each strain cluster (`species|cluster`) sample_fraction, with
  `min_prev=2` prevalence filter.
- Detected-species features: per-species depth from stdout.
- Classifiers: Random Forest (500 trees) and logistic regression (balanced,
  scaled for LR), evaluated with 5-fold stratified CV `cross_val_predict`.
- The combined matrix is row-normalized so that strain fractions + detected depths
  sum to 1 per sample.

## Results

### Lim_ORPI

| Feature set | RF accuracy | RF AUROC | RF F1 | LR AUROC |
|-------------|-------------|----------|-------|----------|
| species-level (1601 species) | 0.670 | 0.653 | 0.484 | 0.640 |
| strain-only (9 strains) | 0.670 | 0.604 | 0.407 | 0.598 |
| **strain + detected depth (9 + 22)** | **0.691** | **0.664** | **0.595** | 0.560 |

Adding detected-not-resolved depth improved AUROC over both strain-only and
species-level baselines in this dataset.

### ECC

| Feature set | RF accuracy | RF AUROC | RF F1 | LR AUROC |
|-------------|-------------|----------|-------|----------|
| species-level (970 species) | 0.816 | **0.936** | 0.821 | 0.873 |
| strain-only (18 strains) | 0.579 | 0.597 | 0.579 | 0.596 |
| **strain + detected depth (18 + 18)** | 0.763 | 0.864 | 0.769 | 0.626 |

Detected depth closed much of the gap between strain-only and species-level
performance, but species-level still outperformed strain+detected for ECC.

## Files

| File | Description |
|------|-------------|
| `results/lim_orpi/detected_depth.tsv` | Per-sample detected-species depth matrix |
| `results/lim_orpi/lim_orpi_ml_v3.out` | Full ML output |
| `results/ecc_strain/detected_depth.tsv` | Per-sample detected-species depth matrix |
| `results/ecc_strain/ecc_ml_v3.out` | Full ML output |
| `scripts/parse_detected_depth.py` | Parse `multi-profile` stdout to depth TSV |
| `scripts/lim_orpi/lim_orpi_ml_compare_v3.py` | Lim_ORPI ML comparison |
| `scripts/ecc_strain/ecc_ml_compare_v3.py` | ECC ML comparison |

## Caveats and next steps

- Current detected-depth matrices were extracted from runs using the default
  `--min-species-markers 50 --min-species-detect 3`. Loosening these thresholds
  to `--min-species-markers 1 --min-species-detect 1` would expose more
  low-coverage detected species, but would require re-profiling all samples.
- Depth and sample_fraction are on different scales; row-normalizing them
  together is a simple first attempt. Better performance might come from
  log-transforming depth or keeping the two feature blocks separately normalized.
- The improvement on Lim_ORPI is promising, but ECC still favors species-level
  features. This suggests strain-level resolution is most useful when the
  species-level signal is already noisy or lower-dimensional.
