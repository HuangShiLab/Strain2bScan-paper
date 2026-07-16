# Fig 11 — Systematic Strain2bScan vs StrainScan head-to-head (15-species simulated benchmark)

**Question.** On a single controlled benchmark — same reference pool, same reads, each tool building its
own database — how do Strain2bScan and StrainScan compare on strain-identification accuracy (precision,
recall, F1 per sample) and on run time and memory, across a range of sequencing depths, for both
single-species and multi-species communities?

## Design
- **Pool:** 15 species, 15–50 complete/near-complete NCBI genomes each (`figure_raw_data/sim_pool_manifest.tsv`).
- **Single-species:** per species, 2/3/5 co-present strains (same- or different-cluster) × depth
  {0.5,1,3,5,10}× × 5 reps = 2 025 samples. Reads: ART `art_illumina`, 150 bp PE.
- **Multi-species:** ~18 co-present species across 3 community depth gradients = 60 samples.
- **Truth** (`figure_raw_data/sim_{single,multi}_species/*/truth/*.tsv`) records species, genome
  accession and 0.95-cluster per strain.

## Tools and fairness
- **Strain2bScan** (native arm64): `cluster --enzyme all --similarity 0.95`; `profile` / `multi-profile
  --enzyme all` (R1+R2 decompressed and concatenated).
- **StrainScan** v1.0.14 (Linux-x86 only) in Docker `linux/amd64` (QEMU on Apple Silicon):
  `strainscan_build`; `strainscan -i R1 -j R2 -d DB`. VM raised to 56 GB (build peaks ~28 GB).
- **Own-cluster-space scoring:** predicted clusters vs truth strains mapped into each tool's clusters
  (Strain2bScan: truth `cluster` column; StrainScan: `Cluster_Result/hclsMap_95_recls.txt`, report
  `Cluster_ID` = `C`+id). P/R/F1 over cluster sets per sample — the honest short-read resolution unit.
- **Coverage:** Strain2bScan profiled all 2 085 samples; StrainScan a matched subset (different-cluster,
  k=2/3/5, rep1, all depths; near-clonal *M. tuberculosis* via same-cluster samples) → 204 paired
  single-species samples over 14 species + 4 multi-species/depth. StrainScan has no multi-species mode, so
  a community's cost = Σ over the 14 species databases (peak RSS = max).
- **Same-environment timing:** a `linux/amd64` Strain2bScan binary was run in the same container on the
  same subset to remove the emulation confound from the profile-speed ratio (Fig 11E, Table 2).

## Results (medians)
- **Accuracy:** both precision 1.0 (Strain2bScan every species; StrainScan 0.80–0.83 on *A. muciniphila*,
  *C. difficile*, *P. dorei*, *S. enterica*). Strain2bScan full recall by 3×, StrainScan only at 10×;
  overall R/F1 0.80/0.89 vs 0.67/0.75. Multi-species accuracy comparable.
- **Build:** Strain2bScan 0.7–5.1 s / 0.1–0.4 GB vs StrainScan 5–43 min / 8–28 GB = **249–614× faster,
  43–138× lighter**. **StrainScan failed to build *K. pneumoniae*** (killed >100 min; stalled >1 h 40 in
  the k-mer-matrix step on retry); Strain2bScan built it in 5.1 s. That species is dropped from the paired
  accuracy set.
- **Profile:** same-env 4–33× faster, 5–39× lighter (0.16–1.7 s / 20–160 MB vs 3.5–6.9 s / ~830 MB).
- **Multi-species:** 1–9 s (one pass) vs 100–398 s (Σ 14 runs) = **50–105× faster**.

## Files
- Figure: `figures/sim_headtohead.{png,pdf}` (script `scripts/plot_sim_headtohead.py`).
- Aggregates: `results/sim_headtohead_{single_by_depth,by_species,multi,build_cost,profile_sameenv}.tsv`.
- Tables: `manuscript/tables.md` (Table 1–3). Per-sample raw + README: `figure_raw_data/sim_headtohead/`.
- Drivers (this session): `scratchpad/eval/run_s2b_{single,multi}.py`,
  `run_strainscan_{single,multi}.py`, `run_s2b_emulated_single.py`, `analyze_headtohead.py`.
