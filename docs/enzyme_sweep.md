# Enzyme count vs strain-level performance (the 2bRAD tuning knob)

**Question.** Does digesting with more type-IIB enzymes help strain-level analysis? More
enzymes → more 2b tags → more strain-specific markers, but also more compute (eroding the
sparsity advantage). Is there an optimum?

**Design.** Cumulative enzyme ladders on the *C. acnes* 64-genome panel (5 simulated
strain-mixture mocks, fixed across levels): 1 (`BcgI`), 2 (`+CspCI`), 4 (`+BaeI,AlfI`),
8 (`+AloI,BsaXI,CjeI,PpiI`), 14 (full curated set). **Strand-fixed binary.**
(`scripts/run_enzyme_sweep.py`, `results/enzyme_sweep.tsv`, `figures/enzyme_sweep.*`.)

| #enzymes | strain-specific markers | clusters | build (s) | precision | recall | Bray–Curtis |
|---|---|---|---|---|---|---|
| 1 (BcgI) | 0 | 1 | 1.05 | 0.00 | 0.00 | 1.00 |
| 2 | 314 | 15 | 1.72 | 1.00 | 0.38 | 0.48 |
| 4 | 7,866 | 16 | 3.37 | 1.00 | 0.75 | 0.25 |
| **8** | 14,332 | 16 | 6.30 | **1.00** | **0.81** | **0.24** |
| 14 | 26,460 | 16 | 10.82 | 1.00 | 0.81 | 0.24 |

## Findings

1. **Precision is 1.0 from 2 enzymes up.** After the strand fix, there is no over-detection —
   the earlier "14 enzymes degrade precision to 0.90" was a strand artifact; now all enzyme
   counts ≥2 give precision 1.0.
2. **Recall is the axis that improves with more enzymes, and it plateaus at ~8.** Recall climbs
   0.38 → 0.75 → **0.81** (2 → 4 → 8 enzymes) and does **not** improve further at 14 (0.81).
   So **~8 enzymes is the sweet spot**: maximum accuracy at ~1.7× less compute than the full 14.
3. **One enzyme is too sparse to resolve strains.** BcgI alone yields so few markers that all
   64 genomes collapse into a single cluster (no strain-specific markers) → detection is
   meaningless (P=R=0). At least ~4 enzymes are needed for useful resolution.
4. **Cost scales ~linearly** with enzyme count (markers 0→26,460; build 1.05→10.82 s). The
   both-strand digestion roughly doubles marker counts vs the pre-fix numbers.

Note the corrected cluster count: 64 *C. acnes* genomes resolve into ~16 clusters (not the
pre-fix 60), consistent with StrainScan's own C. acnes granularity — the 60 was a strand artifact.

## The tuning knob

Enzyme count trades **sparsity (efficiency)** against **marker density (recall + low-depth
sensitivity)**. Precision is robust across the range; recall and low-depth detection improve
with more markers up to ~8 enzymes. Recommended default: a **moderate set (~8 enzymes)** —
near-optimal accuracy at lower cost than the full 14. (With the strand fix, even the depth
onset now matches StrainScan at 0.5×, `docs/depth_sensitivity.md`.)

## Caveats
- Single species (*C. acnes*); cross-species generalization in `docs/cross_species.md`.
- Simulated error-free reads; an error model would slightly favor more enzymes (redundancy).
