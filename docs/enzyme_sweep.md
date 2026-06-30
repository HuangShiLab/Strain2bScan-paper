# Enzyme count vs strain-level performance (the 2bRAD tuning knob)

**Question.** Does digesting with more type-IIB enzymes help strain-level analysis? More
enzymes → more 2b tags → more strain-specific markers, but also more compute (eroding the
sparsity advantage). Is there an optimum?

**Design.** Cumulative enzyme ladders on the real *C. acnes* 64-genome panel (samples fixed):
1 (`BcgI`), 2 (`+CspCI`), 4 (`+BaeI,AlfI`), 8 (`+AloI,BsaXI,CjeI,PpiI`), 14 (full curated set,
HaeIV/Hin4I excluded as degenerate). For each: strain-specific marker count, cluster count
(resolution), accuracy on the 5 samples, and build/profile time + memory.
(`scripts/run_enzyme_sweep.py`, `results/enzyme_sweep.tsv`, `figures/enzyme_sweep.*`.)

| #enzymes | strain-specific markers | clusters | build (s) | profile (s) | precision | recall | Bray–Curtis |
|---|---|---|---|---|---|---|---|
| 1 (BcgI) | 891 | **16** | 0.24 | 0.12 | 1.00\* | 1.00\* | 0.00\* |
| 2 | 1,093 | 51 | 0.40 | 0.15 | 0.00 | 0.00 | 1.00 |
| 4 | 6,118 | 60 | 0.93 | 0.22 | 1.00 | 0.31 | 0.38 |
| **8** | 10,475 | 60 | 1.66 | 0.34 | **1.00** | **0.56** | **0.15** |
| 14 | 18,566 | 60 | 2.81 | 0.55 | 0.90 | 0.56 | 0.25 |

\* The 1-enzyme "perfect" scores are a **coarse-clustering artifact**: BcgI yields only 16
clusters from 64 genomes, so distinct strains are lumped and the cluster-level metric is
trivially satisfied — it is not strain resolution (cf. `docs/acnes_results.md`).

## Findings

1. **≥4 enzymes are needed for full strain resolution.** Clusters saturate at 60 (≈ full
   genome-level resolution) by 4 enzymes; 1–2 enzymes are degenerate — either too coarse
   (1: 16 clusters) or finely split but **under-marked** (2: 51 clusters with too few
   cluster-specific markers each → detection collapses to P=R=0).
2. **~8 enzymes is the accuracy sweet spot** (precision 1.0, recall 0.56, lowest Bray–Curtis
   0.15). Beyond that, markers and compute keep rising but accuracy **plateaus and slightly
   degrades**: at 14 enzymes precision drops to 0.90 and Bray–Curtis rises to 0.25 (extra
   markers add some cross-talk/noise). *Practical note: the default 14-enzyme set is not the
   accuracy optimum here — ~8 is both more accurate and ~1.7× cheaper.*
3. **Cost scales ~linearly** with enzyme count (markers 891→18,566; build 0.24→2.81 s;
   profile 0.12→0.55 s; peak RSS 127→201 MB).

## The tuning knob (ties the low-depth result together)

Enzyme count is the lever that trades **sparsity (efficiency)** against **marker density
(resolution + low-depth sensitivity)**. The low-depth limitation (Strain2bScan needs ≥~5×
because sparse markers run out of statistical power, `docs/depth_sensitivity.md`) and this
result are two views of the same axis: when you need more sensitivity/resolution you add
enzymes (more markers, more compute); when you need maximum throughput you use fewer. The
recommended default is a **moderate set (~8 enzymes)** — near-optimal accuracy at lower cost
than the full 16/14.

## Caveats
- The 2-enzyme P=R=0 point is a genuine "fine-split but under-marked" regime, not a bug.
- A clean "low-depth detection vs #enzymes" curve is confounded here because adding enzymes
  changes the clustering granularity at the same time; the sweep TSV records it but the figure
  omits it. A controlled version (fixed clusters, vary only marker density) is future work.
- Single species (*C. acnes*); the cross-species generalization is in `docs/cross_species.md`.
