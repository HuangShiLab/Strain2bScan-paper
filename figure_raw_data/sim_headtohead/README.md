# Strain2bScan vs StrainScan — head-to-head on simulated metagenomes

Both tools run on the same 15-species reference pool and the same simulated reads
(`../sim_single_species/`, `../sim_multi_species/`), scored **each in its own cluster space**
(the honest unit for short-read strain calling). Strain2bScan ran natively (Apple Silicon);
StrainScan is Linux-x86 only and ran under Docker `linux/amd64` QEMU emulation, so its wall-clock
is inflated — a same-environment ratio is provided from an in-container linux/amd64 Strain2bScan
build (`profile_timing_sameenv.tsv`), and DB-build wall-clock is reported in the emulated
environment for both where noted.

## Datasets
- **Single-species:** 15 species × {same,diff cluster} × k∈{2,3,5} strains × 5 reps × depth
  {0.5,1,3,5,10}× = 2,025 samples. Strain2bScan profiled all 2,025; StrainScan profiled a matched
  subset (diff, k2/3/5, rep1, all depths; near-clonal *M. tuberculosis* via its `same` samples) →
  **204 paired samples across 14 species**.
- **Multi-species:** 60 community samples over 3 depth gradients. Strain2bScan profiled all 60 in a
  single `multi-profile` pass; StrainScan (no multi-species mode) profiled 4 samples/depth by running
  each species DB separately (per-sample cost = Σ over species).

## Headline results (medians)

| metric | Strain2bScan | StrainScan |
|---|---|---|
| single-species precision | **1.00** (all 14 species) | 1.00 (10/14; 0.80–0.83 on Akkermansia, C. difficile, P. dorei, Salmonella) |
| single-species recall | **0.80** (full recall by 3×) | 0.67 (full recall only at 10×) |
| single-species F1 | **0.889** | 0.75 |
| DB build time / species | **0.7–5.1 s** | 5–43 min (**249–614× slower**) |
| DB build memory / species | **0.1–0.4 GB** | 8–28 GB (**43–138× heavier**) |
| profile time / sample (same env) | **0.16–1.7 s** | 3.5–6.9 s (**4–33× slower**) |
| multi-species time / community sample | **1–9 s** (1 pass) | 100–398 s (Σ 14 runs; **50–105× slower**) |

- **Klebsiella pneumoniae** (47 genomes × 5.5 Mb, the heaviest species): StrainScan build did **not
  complete** — killed at the 100-min cap, and on retry it stalled >1h40 in the final matrix step
  (C20 cluster). Strain2bScan built it in **5.1 s / 397 MB**. Excluded from the paired accuracy set.
- **Multi-species accuracy** is comparable (StrainScan marginally higher recall on deep communities);
  the decisive multi-species difference is throughput (one digest-once pass vs one run per species).

## Files
- `headtohead_accuracy.png` — 6-panel summary (single-species P/R/F1 vs depth; build time/species;
  same-env profile time; multi-species time).
- `single_headtohead_by_depth.tsv`, `single_headtohead_by_species.tsv`, `multi_headtohead.tsv`,
  `build_cost.tsv`, `profile_timing_sameenv.tsv` — aggregated head-to-head tables.
- `*_persample.tsv` — raw per-sample P/R/F1 + time/memory for both tools (single, multi, emulated).
- `strain2bscan_build_times.tsv`, `strainscan_build_times.tsv` — per-species DB build cost.

## Reproduce
Scripts in `scratchpad/eval/` (this session): `run_s2b_single.py`, `run_s2b_multi.py`,
`run_strainscan_single.py`, `run_strainscan_multi.py`, `run_s2b_emulated_single.py`,
`analyze_headtohead.py`. StrainScan runs inside Docker container `ss` (strainscan 1.0.14,
linux/amd64); DBs built by `ss_build_driver.sh`. Docker VM raised to 56 GB (StrainScan build peaks
28 GB) by writing `MemoryMiB` to `settings-store.json` while Docker is fully stopped.
