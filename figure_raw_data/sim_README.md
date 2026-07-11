# Simulated benchmark datasets (consistent, StrainScan-style)

One consistent simulated dataset built from the **15-species complete/near-complete genome pool**
(`data/genome_qc_16s_panel.tsv`; pinned cluster membership in `sim_pool_manifest.tsv`), used across all
simulated figures (3–5, 9–10). Reads simulated with **ART** (`art_illumina -p -l 250 -f <cov> -m 600
-s 150`, paired-end Illumina error model — StrainScan's method, PMID 37587527). Design rationale and
audit: `docs/simulation_audit.md`. Generators: `scripts/simulate_single_species.py`,
`scripts/simulate_multi_species.py`, `scripts/build_sim_pool.py`. FASTQ are git-ignored; truth tables
and manifests are tracked.

## A. Single-species, multi-strain — `sim_single_species/<Species>/`
For each of the 15 species: multi-strain communities with **uneven StrainScan coverage profiles**
(2-strain 100:10; 3-strain 100:50:10; 5-strain 100:70:50:20:10), under two difficulty strategies —
strains from **different clusters** (`diff`) and from the **same cluster** (`same`, the hard
near-clonal case) — with **5 reps**, each simulated across the **depth ladder 0.5/1/3/5/10×** (mean
per-strain coverage). The *same* strain community is simulated at every depth.
- `reads/<Species>__<strat>_k<k>_rep<r>_d<D>_R{1,2}.fastq.gz`
- `truth/<...>.truth.tsv` — columns: strain, cluster, relative_abundance, coverage, mean_depth,
  n_strains, strategy.
- Infeasible combinations are skipped and simply absent (e.g. *M. tuberculosis* = 1 cluster → only
  `same`; *P. copri* = all singletons → only `diff`).

## B. Multi-species, multi-strain — `sim_multi_species/depth_{low,med,high}/`
Mixed community of **all 15 species** per sample (1–2 strains/species at log-normal relative
abundance), simulated at **3 depth gradients**: `low` = 1×, `med` = 5×, `high` = 10× mean per-strain
coverage, **20 samples each** (60 total).
- `depth_<grad>/reads/sample<i>_R{1,2}.fastq.gz`
- `depth_<grad>/truth/sample<i>.truth.tsv` — columns: species, strain, cluster, relative_abundance,
  coverage, mean_depth.

## Figure mapping (all from this one dataset)
- **Fig 3A** accuracy → single-species at a fixed depth (5×); **Fig 3B** depth → the 0.5–10× ladder.
- **Fig 4** robustness → panel-size / degradation on one species (P. copri) from the same pool.
- **Fig 5** enzyme knob → one species (C. acnes) from the same pool.
- **Fig 9** scale → multi-species community + its 3 depth gradients.
- **Fig 10** head-to-head → the 3 StrainScan species (A. muciniphila, P. copri, M. tuberculosis),
  which are in the 15-species pool, at the ladder depths.

Consistency: the **same genome pool, cluster definitions, and truth format** underlie every figure.
