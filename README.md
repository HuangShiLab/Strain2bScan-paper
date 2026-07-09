# Strain2bScan — manuscript & reproducibility

Scripts, accession lists, ground truth, result tables, and analysis notes for the
**Strain2bScan** manuscript. The software lives at
[HuangShiLab/Strain2bScan](https://github.com/HuangShiLab/Strain2bScan).

> **Note — strand-invariance digestion fix (all benchmarks re-run).** A core correctness bug was
> found and fixed mid-project: tag extraction was not reverse-complement invariant, so reference
> genomes (digested forward-only) and reads (both strands) lived in different marker spaces,
> causing spurious over-detection of similar strains and corrupting clustering. The fix
> (both-strand digestion) makes accuracy near-perfect (precision 1.0 across the cross-species,
> species-expansion and panel-size benchmarks; low-depth detection now matches StrainScan) at the
> cost of ~2× more markers (so the per-sample speed edge is ~6× rather than the pre-fix ~14×).
> **Every benchmark below has been re-run with the fixed binary**; the numbers here are the
> fixed-binary results.

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
| Multi-species scaling, earlier 40-species panel (species & sample gradients) | `results/multispecies_scaling.tsv` | superseded by the 55-species run below; kept for reference |
| Multi-species scaling, **55-species** panel + community throughput vs StrainScan | `results/multispecies_species_gradient55.tsv`, `..._sample_gradient55.tsv`, `..._accuracy_gate55.tsv`, `community_throughput.tsv` + `figures/community_throughput.*` | `make genomes samples build scaling gate-sweep` (55-species pinned accessions), then `make figures`; see `docs/multispecies.md` |
| Multi-species accuracy vs species gate | `results/multispecies_accuracy_gate.tsv`, `..._accuracy_gate55.tsv` | `make gate-sweep` |
| Reference-genome quality vs accuracy (supp.) | `results/refqual_degradation.tsv` + `figures/refqual_figure.*` | `scripts/run_refqual.py` (Strain2bScan) + `run_refqual_strainscan.py` (Linux); `plot_refqual.py` |
| Enzyme count vs performance (the 2bRAD knob) | `results/enzyme_sweep.tsv` + `figures/enzyme_sweep.*` | `make enzyme-sweep` |
| Cross-species generalization (3 species, our own panels) | `results/cross_species.tsv` + `figures/cross_species.*` | `make cross-species` |
| Species expansion (3 more, StrainScan's own pre-built DBs) | `results/species_expansion.tsv` + `figures/species_expansion.*` | `make species-expansion`; see `docs/species_expansion.md` |
| P. copri panel-size test (does a bigger DB fix it? — no) | `results/panelsize_prevotella.tsv` + `figures/panelsize_prevotella.*` | `make panelsize`; see `docs/species_expansion.md` |

### Quick multi-species run
```bash
export STRAIN2BSCAN_BIN=/path/to/strain2bscan
make genomes      # download 40 species × 4 strains (NCBI)
make samples      # simulate 20 mock metagenomes (seeded) + truth
make build scaling
make gate-sweep   # accuracy vs Layer-1 species gate
```

## Headline results (see `results/` and `docs/` for detail; all strand-fixed)
- **Speed/memory:** ~6× faster, ~7× less memory per sample than StrainScan (*C. acnes*:
  1.25 s / 115 MB vs 7.06 s / 828 MB). The ~6× (vs a pre-fix ~14×) reflects the 2× markers from
  both-strand digestion — a correctness-over-speed tradeoff, recoverable by framing one canonical
  tag per site.
- **Scalability:** per-sample cost independent of species count (digest once, match many),
  flat from 10 to **55 real species** (~4.6 s); linear in samples (~4.1 s/sample, 10→200
  measured, extrapolated to 500); ~5× parallel build speedup.
- **Complex-community efficiency (55 species):** **~91× faster at 100 samples** (6.8 min vs a
  projected 10.4 h) and **~92× at 500 samples** (34 min vs ~52 h) than running StrainScan once
  per species per sample — the only mode available to a tool with no multi-species support.
  See `docs/multispecies.md`.
- **Accuracy (55-species panel):** with the Layer-1 species gate, species precision **0.96→1.0**
  (gate 10→400), species recall 0.98, strain recall 0.98.
- **Cross-species (C. acnes, S. aureus, S. epidermidis):** **precision 1.0 for all three**,
  recall 0.81–1.0, Bray–Curtis 0.01–0.24. The earlier "resolvability gradient" (precision
  0.90→0.48) was a strand-bug artifact; clustering now matches StrainScan's own granularity.
- **Depth:** Strain2bScan now **detects down to 0.5×, matching StrainScan** (both miss at 0.1×).
  The earlier "less sensitive below ~1×" was the strand bug.
- **Enzyme count (the 2bRAD knob):** precision 1.0 from ≥2 enzymes; recall climbs to 0.81 and
  plateaus by **~8 enzymes (the sweet spot)** — 14 is no better. 1 enzyme is too sparse to resolve.
- **Reference quality (supp.):** precision stays **1.0 down to 70% completeness / 8%
  contamination**; recall and abundance error degrade progressively and hit a hard cliff (total
  failure) at 50% completeness / 10% contamination / 400 contigs.
- **Species expansion (StrainScan's own pre-built DBs — *A. muciniphila*, *P. copri*,
  *M. tuberculosis*):** Strain2bScan **precision 1.0** on all three (~10–13× faster, ~7–15×
  lighter where StrainScan completes), matching StrainScan's precision while recall-first. The
  one low recall, *M. tuberculosis* (0.46), reflects **genuine near-clonality** (40 genomes → 8
  clusters), not a bug. **StrainScan itself did not complete** on *M. tuberculosis* (>3.3 h,
  >25.5 GB RAM, 1 sample) — a real scalability ceiling of its full-k-mer regression Layer-2 on
  near-clonal species. See `docs/species_expansion.md`.

## Notes
- Committed scripts are the exact ones used (paths parameterized via `$WORKDIR` /
  `$STRAIN2BSCAN_BIN`). The HaeIV/Hin4I enzymes are excluded from the multi-enzyme set
  (degenerate recognition in the Fast2bRAD-M table → over-matching).
- Reads are simulated error-free; an error model is a planned addition (see `docs/`).
