# Strain2bScan — manuscript & reproducibility

Scripts, accession lists, ground truth, result tables, and analysis notes for the
**Strain2bScan** manuscript. The software lives at
[HuangShiLab/Strain2bScan](https://github.com/HuangShiLab/Strain2bScan).

Large data (genomes, reads, databases) is **not** committed — everything is regenerated from
the committed accession lists + seeded scripts, so the repository stays small and fully
reproducible.

## Layout

```
scripts/    benchmark + simulation scripts (portable: set $WORKDIR and $STRAIN2BSCAN_BIN)
data/
  accessions/   NCBI accession lists  (cae_reference_545.txt, cae_truth_strains.txt,
                                       multispecies_40x4.tsv)
  truth/        per-sample ground-truth abundance tables (C. acnes mocks)
results/    measured result tables (TSV) — the source data for each figure
docs/       analysis write-ups (methods rationale + per-experiment results)
env/        conda env for the StrainScan comparison
figures/    figure-generation scripts + rendered figures (add as the manuscript develops)
Makefile    convenience targets for the multi-species benchmark
```

## Prerequisites
- **Strain2bScan** built: `cargo build --release` in the software repo; then either put
  `target/release/strain2bscan` on `PATH` or `export STRAIN2BSCAN_BIN=/abs/path/strain2bscan`.
- Python 3 (stdlib only for simulation/eval); network access to NCBI Datasets.
- For the StrainScan head-to-head: `conda env create -f env/strainscan-linux.yml` (Linux).

## Data sources
- **Reference genomes** — fetched from NCBI by the accession lists in `data/accessions/`
  (`scripts/dl_multispecies.py`; the *C. acnes* 545-genome panel is in `cae_reference_545.txt`).
- **C. acnes mock reads** — from
  [HuangShiLab/MockMetagenomes4Benchmark](https://github.com/HuangShiLab/MockMetagenomes4Benchmark)
  (`Cutibacterium_acnes_0.01/`).
- **Multi-species samples** — regenerated deterministically (fixed seeds) by
  `scripts/sim_samples.py`; no reads are stored.

## Reproduce each result

| result / figure | table | how |
|---|---|---|
| Strain2bScan vs StrainScan, per-sample time & memory | `results/headtohead_performance.tsv` | `scripts/run_acnes.sh` (Strain2bScan) + StrainScan on the 5 *C. acnes* mocks |
| Parallel speedup + DB-build scaling | `results/parallel_and_build_scaling.tsv` | `STRAIN2BSCAN_THREADS=1\|16 strain2bscan cluster …`; `scripts/gen_scale.py` for synthetic panels |
| Detection vs depth (0.1–5×) | `results/depth_sensitivity.tsv` | `scripts/sim_depth.py` then profile with each tool |
| Multi-species scaling (species & sample gradients) | `results/multispecies_scaling.tsv` | `make genomes samples build scaling` |
| Multi-species accuracy vs species gate | `results/multispecies_accuracy_gate.tsv` | `make gate-sweep` |
| Reference-genome quality vs accuracy (supp.) | `results/refqual_degradation.tsv` + `figures/refqual_figure.*` | `scripts/run_refqual.py` (Strain2bScan) + `run_refqual_strainscan.py` (Linux); `plot_refqual.py` |
| Enzyme count vs performance (the 2bRAD knob) | `results/enzyme_sweep.tsv` + `figures/enzyme_sweep.*` | `make enzyme-sweep` |
| Cross-species generalization (3 species) | `results/cross_species.tsv` + `figures/cross_species.*` | `make cross-species` |

### Quick multi-species run
```bash
export STRAIN2BSCAN_BIN=/path/to/strain2bscan
make genomes      # download 40 species × 4 strains (NCBI)
make samples      # simulate 20 mock metagenomes (seeded) + truth
make build scaling
make gate-sweep   # accuracy vs Layer-1 species gate
```

## Headline results (see `results/` and `docs/` for detail)
- **Speed/memory:** ~14× faster, ~7× less memory per sample than StrainScan (real *C. acnes*).
- **Scalability:** per-sample cost independent of species count (digest once, match many);
  linear in samples; ~9× parallel build speedup; MinHash clustering identical to exact.
- **Accuracy (40-species panel):** with the Layer-1 species gate, species precision 0.98 /
  recall 1.0 and strain recall 1.0; without it, precision collapses (cross-species sharing) —
  evidence that the two-layer (species → strain) design is necessary.
- **Depth:** below ~1× per strain, full-k-mer StrainScan is more sensitive; Strain2bScan
  matches it at sufficient depth while being far faster/lighter.
- **Reference quality (supp.):** degrading the reference DB (completeness↓, contamination↑,
  fragmentation↑) lowers precision (0.90→0.64) and roughly doubles abundance error
  (Bray–Curtis 0.25→0.46) by ≤50% completeness / 10% contamination — Jaccard clustering is
  completeness-sensitive, motivating the `quality.rs` filters.
- **Enzyme count (the 2bRAD knob):** ≥4 enzymes reach full strain resolution; **~8 is the
  accuracy sweet spot** (precision 1.0, Bray–Curtis 0.15) — more enzymes add cost and noise
  with no accuracy gain (14: precision 0.90). Enzyme count trades efficiency ↔ resolution.
- **Cross-species:** accuracy tracks intra-species **resolvability** (clusters/genomes): high
  for diverse *C. acnes* (0.94 → precision 0.90), lower for similar-strain *S. epidermidis*
  (0.20 → precision 0.59). Occurrence-based uniqueness (a marker unique iff present at any copy
  in one cluster) lifts similar-strain precision/abundance; a residual gap needs overlap-aware
  deconvolution.

## Notes
- Committed scripts are the exact ones used (paths parameterized via `$WORKDIR` /
  `$STRAIN2BSCAN_BIN`). The HaeIV/Hin4I enzymes are excluded from the multi-enzyme set
  (degenerate recognition in the Fast2bRAD-M table → over-matching).
- Reads are simulated error-free; an error model is a planned addition (see `docs/`).
