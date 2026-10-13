# Current paper output map

This is the short index for the 2026-10-12 manuscript snapshot (`paper-2026-10-13-submission`).

## Manuscript

- Editable modular sources: `manuscript/abstract.md`, `introduction.md`, `results.md`,
  `discussion.md`, `methods.md`, `figures.md`, `tables.md`, `references.md`.
- Assembled paper: `manuscript/full_manuscript.md`.
- Thesis version: `manuscript/thesis_chapter.md`.
- Regenerate both assembled files with `python3 scripts/assemble_manuscript.py`.

## Final figures

- Publication set: `figures/numbered/Fig01...Fig11` and `FigS1...FigS5` (PNG + PDF).
- Checksums and source mapping: `figures/numbered/MANIFEST.tsv`.
- Current source panels: `figures/source_panels/`.
- Regenerate numbered images with `python3 scripts/number_figures.py`, then
  `python3 scripts/make_figure_manifest.py`.

## Current primary evidence tables

- Mock benchmark: `results/mock_benchmark_f26f234/` and `data/fig6_fig12_metrics.tsv`.
- StrainScan rerun: `results/strainscan_rerun/`.
- Current software regression: `results/strainscan_rerun/current_single_persample.tsv`, `current_uncertainty_summary.tsv` and `current_software_regression_provenance.json`.
- Cluster-resolution audit: `results/strainscan_rerun/database_cluster_resolution.tsv`.
- StrainGE rerun: `results/strainge_rerun/`.
- Saliva analyses: `results/saliva_*.tsv`.
- Public-cohort analyses: `results/realworld_cohort_benchmark/`.
- Frozen configuration and runtime metadata:
  `results/benchmark_configuration.json` and `results/runtime_environment.json`.

## QA entry points

- `python3 scripts/check_manuscript_claims.py`
- `python3 scripts/assemble_manuscript.py`
- `python3 scripts/number_figures.py`
- `python3 scripts/make_figure_manifest.py`

Historical diagnostics, planning notes, superseded figures, and bulk local reads are under
`archive/`. See `archive/README.md`.
