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
| 1 (BcgI) | 882 | 16 | 0.58 | 1.00 | 0.56 | 0.25 |
| 2 | 1,029 | 16 | 1.02 | 1.00 | 0.50 | 0.35 |
| **4** | 1,524 | 16 | 1.77 | **1.00** | **0.75** | **0.23** |
| 8 | 4,036 | 16 | 3.97 | 1.00 | 0.75 | 0.24 |
| 14 | 7,400 | 16 | 6.38 | 1.00 | 0.75 | 0.24 |

## Findings

1. **Precision is 1.0 at every enzyme count, including BcgI alone.** With the correct tag lengths
   there is no over-detection at any point on the ladder.
2. **BcgI alone already resolves and profiles** (16 clusters, precision 1.0, recall 0.56) — the
   canonical single-enzyme 2bRAD workflow works out of the box.
3. **Recall improves with more enzymes and plateaus by ~4.** Recall rises 0.56 → 0.75 (1 → 4
   enzymes) and does **not** improve further at 8 or 14 (0.75). So **~4 enzymes is the sweet
   spot**: maximum accuracy at a fraction of the 14-enzyme cost. (Cluster resolution is already
   saturated at 16 from one enzyme — adding enzymes adds *marker density* / recall, not clusters.)
4. **Cost scales ~linearly** with enzyme count (strain-specific markers 882→7,400; build
   0.58→6.38 s). The correct (Fast2bRAD-M) tag lengths give a sparser, correct marker set —
   ~3.5× fewer markers than the earlier both-strand workaround, so builds are faster too.

The cluster count (64 *C. acnes* genomes → 16 clusters) is consistent with StrainScan's own
C. acnes granularity.

## The tuning knob

Enzyme count trades **sparsity (efficiency)** against **marker density (recall + low-depth
sensitivity)**. Precision is robust across the range; recall and low-depth detection improve
with more markers up to ~8 enzymes. Recommended default: a **moderate set (~8 enzymes)** —
near-optimal accuracy at lower cost than the full 14. (With the strand fix, even the depth
onset now matches StrainScan at 0.5×, `docs/depth_sensitivity.md`.)

## Caveats
- Single species (*C. acnes*); cross-species generalization in `docs/cross_species.md`.
- Simulated error-free reads; an error model would slightly favor more enzymes (redundancy).
