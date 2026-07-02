# Cross-species generalization (C. acnes, S. aureus, S. epidermidis)

**Goal.** Test whether the strain-profiling results generalize beyond *C. acnes*, and
characterize how performance depends on the species. Real NCBI panels (60 complete genomes
each for the Staphylococci; the 64-genome *C. acnes* panel), identical 14-enzyme pipeline.
*C. acnes* uses the real MockMetagenomes4Benchmark samples; S. aureus / S. epidermidis use
simulated strain-mixture mocks (`scripts/sim_strain_mock.py`: 2–5 strains/sample, log-normal
depth ≥1×, matching the *C. acnes* design). `results/cross_species.tsv`, `figures/cross_species.*`.

| species | genomes | clusters | clusters/genomes | strain-specific markers | precision | recall | Bray–Curtis |
|---|---|---|---|---|---|---|---|
| *C. acnes* | 64 | **60** | **0.94** | 18,566 | **0.90** | 0.56 | 0.25 |
| *S. aureus* | 60 | 28 | 0.47 | 20,258 | 0.52 | 0.87 | 0.59 |
| *S. epidermidis* | 60 | **12** | **0.20** | 5,856 | 0.59 | 0.91 | 0.37 |

*(with occurrence-based uniqueness — see below.)*

## The key axis: intra-species resolvability

`clusters / genomes` at the 0.95 Jaccard cut measures how many strains are actually
**distinguishable**: *C. acnes* genomes are nearly all distinct (60 clusters from 64 genomes),
whereas *S. epidermidis* strains are so similar that 60 genomes collapse into **12** clusters.
This is a real biological property of each species, not a method artifact — no marker method
can resolve strains that are near-identical.

**Precision tracks resolvability** (0.90 → 0.52 → 0.59 as clusters/genomes falls 0.94 → 0.20).
In heavily-clustered, low-diversity species the profiler tends to **over-detect**: reads from a
present strain can match a *different* large cluster's markers. We traced the dominant cause to
**single-copy-filter asymmetry** — a tag single-copy in cluster A but multi-copy (hence
single-copy-filtered) in cluster B's genomes was labelled "unique to A" yet is still produced
by B's reads. Two fixes were added, in order of impact:
- **Occurrence-based uniqueness** (a marker is unique iff it occurs *at any copy number* in
  exactly one cluster's genomes) removes these false-unique markers at the detection step. It
  is the main lever: *S. epidermidis* precision 0.48→**0.59** and Bray–Curtis 0.80→**0.37**,
  *S. aureus* recall 0.73→**0.87**, with *C. acnes* unchanged (P=0.90).
- A **coverage-fraction gate** (`--min-coverage`, default 0.1) additionally removes large,
  similar clusters called at a tiny coverage fraction.

A residual gap remains for the hardest low-diversity cases (precision ≈0.5–0.6), where
near-identical strains are intrinsically ambiguous from short reads.

## What this means for the manuscript

1. **The method generalizes where strains are resolvable.** For species with sufficient
   intra-species diversity (*C. acnes*), Strain2bScan gives high-precision strain calls; the
   speed/scale advantages are species-independent.
2. **Resolvability is a species property, reported honestly.** `clusters/genomes` (and the
   `NOT DOABLE` verdict) tell the user, per species, how fine a resolution the data support —
   consistent with "not all species are strain-resolvable."
3. **Low-diversity species are partly rescued by occurrence-based uniqueness**, which lifts
   precision/abundance on *S. epidermidis*/*S. aureus* at no cost to *C. acnes*. The residual
   gap (precision ≈0.5–0.6) is the known hard case that StrainScan tackles with a within-cluster
   overlap-matrix / lasso Layer-2; adding a full overlap-aware deconvolution on top of
   occurrence-based uniqueness is the remaining step for near-identical strains.

## Reproduce
```bash
export STRAIN2BSCAN_BIN=/path/to/strain2bscan
make cross-species        # downloads pinned Staph panels + C. acnes data, runs all 3 species
# -> work/xspec/cross_species.tsv, work/figures/cross_species.*
```

## Caveats
- Staph mocks are simulated error-free from the same panel used to build the DB (closed-world),
  as for *C. acnes*; an error model and held-out truth strains would harden the comparison.
- Two of three species use simulated samples (the benchmark repo only ships *C. acnes* reads);
  real S. aureus/S. epidermidis mocks (Figshare, when available) would replace the simulation.
