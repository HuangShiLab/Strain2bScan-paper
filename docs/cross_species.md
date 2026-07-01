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
| *S. aureus* | 60 | 28 | 0.47 | 20,258 | 0.50 | 0.73 | 0.58 |
| *S. epidermidis* | 60 | **12** | **0.20** | 5,856 | 0.48 | 0.91 | 0.80 |

## The key axis: intra-species resolvability

`clusters / genomes` at the 0.95 Jaccard cut measures how many strains are actually
**distinguishable**: *C. acnes* genomes are nearly all distinct (60 clusters from 64 genomes),
whereas *S. epidermidis* strains are so similar that 60 genomes collapse into **12** clusters.
This is a real biological property of each species, not a method artifact — no marker method
can resolve strains that are near-identical.

**Precision tracks resolvability** (0.90 → 0.50 → 0.48 as clusters/genomes falls 0.94 → 0.20).
In heavily-clustered, low-diversity species the simple Layer-2 **over-detects**: because unique
markers are scarce and clusters are large, reads from a present strain can match a small
fraction of a *different* large cluster's "unique" markers (single-copy-filter asymmetry under
high similarity), producing false positives and inflated abundance error. A coverage-fraction
gate (`--min-coverage`, default 0.1) removes the most egregious cases (S. aureus precision
0.40→0.50, S. epidermidis 0.37→0.48; *C. acnes* unchanged) but does not fully solve it.

## What this means for the manuscript

1. **The method generalizes where strains are resolvable.** For species with sufficient
   intra-species diversity (*C. acnes*), Strain2bScan gives high-precision strain calls; the
   speed/scale advantages are species-independent.
2. **Resolvability is a species property, reported honestly.** `clusters/genomes` (and the
   `NOT DOABLE` verdict) tell the user, per species, how fine a resolution the data support —
   consistent with "not all species are strain-resolvable."
3. **Low-diversity species need overlap-aware deconvolution.** The over-detection on
   *S. epidermidis*/*S. aureus* is the known hard case that StrainScan addresses with its
   within-cluster overlap-matrix / lasso Layer-2. Porting that (beyond the current
   detect-and-depth Layer-2) is the clear next step to lift precision on similar-strain species.

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
