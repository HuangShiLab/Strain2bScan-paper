# Strain2bScan — Results (draft)

> **⚠ Narrative + numbers below PREDATE the strand-invariance digestion fix and need revision.**
> A core bug (forward-only tag extraction) was found and fixed after this draft. It changes the
> story substantially — for the better: accuracy is now precision 1.0 across species (the
> "resolvability gradient" and "low-diversity over-detection" here were bug artifacts), low-depth
> detection now matches StrainScan at 0.5×, and the per-sample speed edge is ~6× (not ~14×; the
> DB doubled). The authoritative fixed numbers/figures are in `results/*.tsv`, `figures/*`, and
> the updated `docs/*` (cross_species, multispecies, depth_sensitivity, enzyme_sweep,
> refgenome_quality, species_expansion). This prose should be rewritten around: *accurate +
> fast/light + community-scale + 2bRAD-capable*, dropping the "operating envelope" framing.

## Faster and lighter than StrainScan per sample

On five real *C. acnes* mock samples profiled against a 64-genome panel, Strain2bScan and
StrainScan detected the same dominant strains, but Strain2bScan profiled each sample in
**0.50 s using ~120 MB** of memory versus **~7.0 s and 828 MB** for StrainScan — **~14× faster
and ~7× lighter** (Fig 1a). The memory gap is structural: StrainScan counts the full k-mer set
with jellyfish (~800 MB hash), whereas Strain2bScan streams reads and keeps only the sparse
2b-tag markers. Database construction and profiling parallelise cleanly (dependency-free `std`
threads): building the 64-genome database sped up 8.5× (22 → 2.6 s) and profiling 6.5× (3.4 →
0.5 s) from 1 to 16 threads (Fig 1b).

## Scaling to many species and many samples

The per-sample cost is **independent of the number of species**: profiling one sample against
10, 20, 30 or 40 per-species databases took ~1.0–1.1 s throughout (Fig 2a), because the sample
is digested once and each additional species is a cheap hash-set match rather than a re-count.
Throughput is **linear in samples** (~1.2 s/sample; Fig 2b), and database build scales
near-linearly with genome count once MinHash-sketch clustering replaces exact all-pairs Jaccard
above 96 genomes (Fig 2c), with MinHash yielding partitions identical to exact on real data.

This is decisive for multi-species communities. Because full k-mer tools have no multi-species
mode, resolving *S* species requires running them once per species database, re-counting the
sample's k-mers each time; Strain2bScan pays the read-digestion cost once. For a 40-species
panel, profiling *N* samples costs Strain2bScan ~20 s (build) + *N* × 1.2 s, versus a projected
*N* × 40 × ~6.8 s for StrainScan run per species. At 100 samples this is **~2 minutes versus
~7.5 hours (~200×)**; the advantage grows linearly with the number of species (Discussion). A
Layer-1 species gate — requiring a minimum number of species-specific markers before a species
is profiled — keeps community accuracy high: across 20 simulated 8-species samples, species
recall stayed at 1.0 for every gate setting while precision rose to 0.94–0.99 as the gate
tightened (Fig 2d), confirming that a species→strain two-layer design is both necessary and
sufficient.

## Accuracy envelope: sequencing depth

Sensitivity depends on per-strain depth. For a single *C. acnes* strain present in both tools'
databases, StrainScan detected it from 0.5× coverage, whereas Strain2bScan required ~1×
(permissive threshold) to ~5× (default) (Fig 3). This is expected: the number of observed
unique markers scales with depth × the number of marker loci, and the 2bRAD set has ~50–100×
fewer loci, so at very low depth it runs out of statistical power first. Below ~1× per-strain
depth, therefore, full k-mer StrainScan is more sensitive; Strain2bScan matches it at
sufficient depth (≈≥5×) while being far faster and lighter.

## Enzyme count tunes resolution against cost

Because more enzymes yield more tags, the enzyme set is a tunable knob (Fig 4). On *C. acnes*,
≥4 enzymes were needed to reach full strain resolution (60 clusters from 64 genomes); one or
two enzymes were degenerate (too coarse, or finely split but under-marked). Accuracy peaked at
~8 enzymes (precision 1.0, Bray–Curtis 0.15) and did not improve — indeed slightly declined —
with 14 enzymes (precision 0.90), while marker count and compute continued to rise ~linearly.
A moderate set (~8 enzymes) is therefore both more accurate and cheaper than the full 16.

## Cross-species generalization and resolvability

Across *C. acnes*, *Staphylococcus aureus* and *S. epidermidis* (real panels; simulated strain
mixtures), accuracy tracked a species-intrinsic property we call **resolvability** — the
fraction of genomes that form distinct clusters at the 0.95 similarity cut (Fig 5). *C. acnes*
genomes are nearly all distinct (60/64 = 0.94; precision 0.90), whereas *S. epidermidis*
strains are so similar that 60 genomes collapse into 12 clusters (0.20). We found that
low-diversity species were prone to over-detection caused by single-copy-filter asymmetry — a
tag single-copy in one cluster but multi-copy (and thus filtered) in another was mislabelled
unique yet reachable from the latter's reads. Defining uniqueness by **full-genome occurrence**
(a marker unique iff present, at any copy number, in exactly one cluster's genomes) removed
these false-unique markers and improved the hard cases without affecting *C. acnes*: *S.
epidermidis* precision rose 0.48→0.59 and abundance error (Bray–Curtis) fell 0.80→0.37, and
*S. aureus* recall rose 0.73→0.87. A residual gap (precision ≈0.5–0.6) reflects the intrinsic
ambiguity of near-identical strains from short reads.

## Robustness to reference-genome quality

Because Jaccard clustering is completeness-sensitive, we degraded the reference genomes of the
truth strains (completeness 100→50 %, contamination 0→10 %, fragmentation 1→400 contigs) while
holding samples fixed (Supp Fig). Precision fell from 0.90 to 0.64 and abundance error roughly
doubled (Bray–Curtis 0.25→0.46) by ≤50 % completeness / 10 % contamination, motivating the
built-in assembly-quality filters (contig count and single-copy tag count relative to the
conspecific median) that flag or remove low-quality references before clustering.
