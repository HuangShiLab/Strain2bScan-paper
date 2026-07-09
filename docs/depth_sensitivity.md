# Depth-sensitivity: Strain2bScan vs StrainScan (the <1× question)

**Question:** at low per-strain depth, is Strain2bScan less sensitive than StrainScan?

**Answer (strand-fixed): No — they match.** Both detect down to 0.5× and both miss at 0.1×.

> This doc previously reported that StrainScan was *more* sensitive (Strain2bScan onset ~1–5×
> vs StrainScan 0.5×) and framed it as a fundamental 2bRAD marker-sparsity limit. That was an
> artifact of the **strand bug**: forward-only digestion recovered only ~half the strain's tags,
> so at low depth too few were observed. With both-strand digestion (the fix), the observed
> marker count per unit depth roughly doubles, and the low-depth gap disappears.

## Experiment
Single C. acnes strain **GCF_009737125** (present in both tools' DBs), simulated shotgun reads
(150 bp, error-free) at 0.1× / 0.5× / 1× / 5× coverage. Strain2bScan profiled against the
64-genome panel (detection floors 10 and 3); StrainScan against its 275-strain DB with its
low-depth modes. `results/depth_sensitivity.tsv`, `figures/depth_sensitivity.*`.

| depth | Strain2bScan (floor 10) | Strain2bScan (floor 3) | StrainScan |
|---|---|---|---|
| 0.1× | miss | miss | miss |
| **0.5×** | **DETECT** | **DETECT** | **DETECT** |
| 1× | DETECT | DETECT | DETECT |
| 5× | DETECT | DETECT | DETECT |

**Detection onset is 0.5× for both tools** (was 1–5× for Strain2bScan pre-fix). Even at the
default floor of 10, Strain2bScan now detects at 0.5×.

## For the manuscript
- **Low-depth sensitivity now matches StrainScan** — the earlier "less sensitive below ~1×"
  caveat is removed; it was the strand bug, not marker sparsity.
- Combined with the ~6× speed / ~7× memory advantage, Strain2bScan is now equal-sensitivity and
  far lighter across the tested depth range (≥0.5×).
- Both tools miss at 0.1× — below ~0.5× per strain, neither has enough signal; this is a shared
  limit of low-coverage strain detection, not a Strain2bScan-specific one.

## To strengthen for publication
- Replicate across strains/species and add ≥3 replicates per depth; report precision/recall.
- Add a realistic error model.
- Same-panel StrainScan build (Linux/dashing) for a fully matched head-to-head.
