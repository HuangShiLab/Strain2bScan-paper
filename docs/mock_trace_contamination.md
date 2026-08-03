# Mock trace contamination and the trace-gap filter

## Observation

On the WMS mock communities MSA-1005 and MSA-1007, Strain2bScan's default output
contains 11–21 false-positive clusters per sample while retaining all true
strains (recall = 1.0). Almost all of these false positives are other mocks'
ATCC type strains detected at 0.01–0.02% relative abundance.

Example breakdown per sample (WMS, 164-panel, Strain2bScan default):

| Mock    | FP/sample | other-mock ATCC | same-species decoy | other |
|---------|-----------|-----------------|--------------------|-------|
| MSA1005 | 11–13     | 11–12           | 0–1                | 0     |
| MSA1007 | 15–21     | 7–12            | 7–9                | 0     |
| MSA1003 | 9–13      | 8               | 1–5                | 0     |

The "other-mock ATCC" calls have coverage/depth/consistency ratios consistent
with genuinely present genomes at ~0.1–1× depth, not with shared-marker shadows.
Because the reference database is built from exactly the four mock communities,
these trace reads map to bona fide ATCC reference genomes rather than random
contaminants.

Interpretation: the signal is cross-library contamination (index hopping or
carry-over during multiplexed library preparation), not an algorithmic
false-positive. The software is sensitive enough to detect it; the default
output therefore reports it.

## Trace-gap filtering

A fixed abundance floor (e.g. `--min-abundance 0.001`) removes most trace false
positives, but it also deletes the rare tail of staggered mocks such as
MSA-1003, whose design includes true members below 0.02%.

A sample-adaptive alternative is implemented as `--trace-gap R --trace-floor F`:

1. Sort detected strains by relative abundance.
2. Find the largest ratio between consecutive abundances.
3. If that ratio exceeds `R`, keep only the strains above the gap.
4. Regardless of the gap, drop strains below `F`.

On the four WMS mocks, `--trace-gap 10 --trace-floor 1e-4` gives:

| Mock    | default P/R/F1 | trace-gap P/R/F1 | FP removed |
|---------|----------------|------------------|------------|
| MSA1002 | 1.00/1.00/1.00 | 1.00/1.00/1.00   | 0/0        |
| MSA1003 | 0.65/1.00/0.79 | 0.94/1.00/0.97   | 28/32      |
| MSA1005 | 0.33/1.00/0.50 | 1.00/1.00/1.00   | 36/36      |
| MSA1007 | 0.25/1.00/0.40 | 1.00/1.00/1.00   | 54/54      |

Recall is preserved because the true community and the trace tail are separated
by more than 10×, while MSA-1003's staggered true members are separated only by
their design ratio.

## Why this is not a default

The trace-gap filter should **not** be enabled by default on real metagenomes.
Real communities typically have long-tailed abundance distributions; there is no
reliable gap between "true rare members" and "contamination". Cutting at the
largest abundance gap would arbitrarily remove real low-abundance organisms.

Recommended practice:

- Report the default output (all detected calls) as the primary result.
- Use `--trace-gap` / `--trace-floor` only when the sample is a defined
  community, or when negative/blank controls and technical replicates justify a
  contamination threshold.
- For open communities, rely on experimental controls and cross-sample
deccntamination rather than a per-sample abundance knee.

## Metrics for mock figures

The main mock figures therefore report:

- **AUPR** as the primary summary metric (threshold-independent).
- **Precision/recall/F1 at abundance ≥ 1e-4** as a fixed-threshold operating
  point.
- The trace-gap-filtered result is shown separately as a supplementary analysis,
  not as the headline performance of the algorithm.
