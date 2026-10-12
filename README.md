# Strain2bScan — manuscript and reproducibility snapshot

This repository contains the manuscript, final figures, current result tables, and reproducibility
scripts for **Strain2bScan**. The software source is in
[HuangShiLab/Strain2bScan](https://github.com/HuangShiLab/Strain2bScan).

The current snapshot is tagged `paper-2026-10-12`.

## Start here

- `PAPER_OUTPUTS.md`: short map to the active manuscript, figures, evidence tables, and QA scripts.
- `manuscript/full_manuscript.md`: assembled paper.
- `manuscript/thesis_chapter.md`: extended thesis version.
- `figures/numbered/`: final Figures 1–11 and Supplementary Figures S1–S5, with checksums.
- `results/`: current evidence tables and frozen benchmark provenance.
- `archive/README.md`: policy for superseded outputs and local bulk intermediates.

## QA and regeneration

```bash
python3 scripts/assemble_manuscript.py
python3 scripts/check_manuscript_claims.py
make figures
```

`make figures` regenerates current source panels and the numbered publication set, then refreshes
`figures/numbered/MANIFEST.tsv`. Frozen benchmark tables and checksums should not be overwritten by
ad hoc runs.

## Data and code availability

- Strain2bScan source: <https://github.com/HuangShiLab/Strain2bScan>.
- This paper snapshot: <https://github.com/HuangShiLab/Strain2bScan-paper>, tag `paper-2026-10-12`.
- Public reads: PRJNA1131785, PRJNA288562, PRJNA1517970, PRJNA1191223 and PRJNA1191225.
- The 22 selected public-library accessions and SRA metadata are in
  `results/realworld_cohort_benchmark/public_run_manifest.tsv`.

Correspondence: Shi Huang, `shihuang@hku.hk`.
