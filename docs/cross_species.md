# Cross-species generalization (C. acnes, S. aureus, S. epidermidis)

**Goal.** Test whether strain profiling generalizes across species. Real NCBI/ENA panels
(*C. acnes* 64 genomes, *S. aureus* 60, *S. epidermidis* 60), identical 14-enzyme pipeline,
5 simulated strain-mixture mocks per species (2–5 strains, log-normal depth ≥1×).
`results/cross_species.tsv`, `figures/cross_species.*`, `scripts/run_cross_species.py`.

> **This experiment was re-run after the strand-invariance digestion fix** (see
> `docs/species_expansion.md` and the software repo). The earlier version of this doc reported a
> "resolvability gradient" in which precision fell from *C. acnes* (0.90) to *S. epidermidis*
> (0.48), and attributed it to clustering collapse in low-diversity species. **That gradient was
> an artifact of the strand bug**, not biology — see below.

## Result (strand-fixed)

| species | genomes | clusters | precision | recall | Bray–Curtis |
|---|---|---|---|---|---|
| *C. acnes* | 64 | 16 | **1.00** | 0.75 | 0.24 |
| *S. aureus* | 60 | 17 | **1.00** | 0.79 | 0.33 |
| *S. epidermidis* | 60 | 11 | **1.00** | **1.00** | **0.01** |

**Precision is 1.0 for all three species** — no false positives. Recall is high (0.75–1.0) and
abundance error (Bray–Curtis) is low (0.02–0.33). *S. epidermidis*, which in the buggy run
looked like the worst case (precision 0.48, Bray–Curtis 0.80), is now the *best* (1.00 / 0.01).

Before → after the fix:

| species | precision | recall | Bray–Curtis |
|---|---|---|---|
| *C. acnes* | 0.90 → **1.00** | 0.56 → 0.75 | 0.25 → 0.24 |
| *S. aureus* | 0.50 → **1.00** | 0.73 → 0.79 | 0.59 → 0.33 |
| *S. epidermidis* | 0.48 → **1.00** | 0.91 → 1.00 | 0.80 → **0.01** |

## What changed and why

The earlier "accuracy tracks intra-species **resolvability** (clusters/genomes)" story rested on
a cluster count that was itself corrupted: forward-only digestion made genome assemblies in
opposite strand orientations look nearly disjoint (Jaccard 0.64 for two ~identical genomes),
so near-identical genomes were split into separate clusters (*C. acnes* 64→**60** clusters) and
each split cluster's spurious "unique" markers cross-detected. With the correct (Fast2bRAD-M) tag
lengths — which restore the forward/reverse patterns as reverse-complement pairs, making a
single-strand scan strand-invariant — the
clustering is correct — *C. acnes* 64→**16**, *S. aureus* 60→**17**, *S. epidermidis* 60→**11**
— now consistent with StrainScan's own granularity (its *C. acnes* DB clusters 275 genomes into
28; its *S. epidermidis* 995 into 378). The clusters/genomes "gradient" (0.94/0.47/0.20) is gone
(now 0.25/0.28/0.18), and precision is uniformly 1.0.

**Occurrence-based uniqueness and the coverage-fraction gate**, which the earlier doc credited
with partially mitigating the over-detection, were treating a symptom; the real cause was the
strand bug. They remain in place and are harmless, but are no longer load-bearing for these
species.

## For the manuscript

1. **Strain profiling generalizes cleanly across species**: precision 1.0 on all three, good
   recall, low abundance error — a much stronger and simpler claim than the earlier gradient.
2. **Recall, not precision, is now the species-dependent axis**, and it reflects *genuine*
   intra-species diversity: where strains are near-clonal (see *M. tuberculosis* in
   `docs/species_expansion.md`, recall ~0.25), cluster resolution is intrinsically coarse; where
   strains are distinct, recall is high. This is honest resolvability, no longer confounded by a
   digestion artifact.

## Caveats

- Simulated error-free reads, closed-world (truth strains in the DB). An error model and
  held-out strains would harden it.
- Genomes fetched from EBI/ENA (NCBI was IP-blocked during the re-run); GCA vs GCF assemblies
  are near-identical for these prokaryotic WGS genomes.
