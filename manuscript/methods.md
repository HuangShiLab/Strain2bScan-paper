# Strain2bScan — Methods (draft)

## Overview

Strain2bScan reimplements the two-layer StrainScan strategy — cluster near-identical strains,
then score samples on markers unique to a strain or cluster — but replaces the full k-mer set
with **2bRAD tags**, and is written in dependency-free Rust for speed and parallelism. The
pipeline is: (i) digest reference genomes and sample reads into 2bRAD-tag markers; (ii) build,
per species, a within-species cluster database annotated with unique markers; (iii) profile a
sample by detecting present clusters from their unique markers and estimating their abundance;
(iv) at the community level, digest each sample once and match it against every per-species
database, gated by a species-level Layer-1 check.

## 2bRAD tag extraction

Type-IIB restriction enzymes cut on both sides of their recognition site, releasing a
fixed-length fragment (the 2bRAD tag, 32–38 bp). Each of the 16 enzymes in the Fast2bRAD-M
table is modelled as a set of anchored sequence patterns; scanning every offset of a sequence
and testing the anchors reproduces the enzyme's digestion sites (forward and reverse patterns
allow scanning a single strand). Each tag is canonicalised (the lexicographically smaller of
the tag and its reverse complement) and hashed to a 64-bit integer marker (FNV-1a; genome and
sample tags use the same hash, so marker values are internally consistent). Two input modes
are supported: **BcgI only** for 2bRAD experimental libraries whose reads already are tags, or
a user-chosen enzyme set (`--enzyme all` for all 16) for in-silico digestion of conventional
shotgun reads, which enriches the marker set ~*n*-fold for *n* enzymes. For reference genomes
we retain only **single-copy** tags (occurring exactly once in the genome), following
StrainScan's and Fast2bRAD-M's use of single-copy markers for unbiased quantification.

## Reference database construction

**Within-species clustering.** For each species, genomes are grouped by single-linkage
hierarchical clustering at 0.95 marker-set similarity (0.05 distance), matching StrainScan's
`hclsMap_95`. Single-linkage at threshold τ is exactly the connected components of the graph
whose edges join genome pairs with Jaccard ≥ τ, computed with union-find. For panels of ≤96
genomes we use exact all-pairs Jaccard on the tag sets; above that we estimate Jaccard from
bottom-*k* MinHash sketches (*k* = 2000) of each genome's markers, which reduces the pairwise
cost from O(n²·m) to O(n²·k) with *k* ≪ *m* and yields partitions identical to exact on real
data (Results). Clusters are the finest reliable resolution unit: strains within one cluster
are too similar to separate from short reads.

**Marker classification.** Within a species, each tag is labelled by its within-species
incidence — present in all clusters (*species-core*; detects the species, not strains), in one
cluster with ≥2 genomes (*cluster-specific*), in a single genome (*strain-specific*), or in
several but not all clusters (*shared-partial*). Cluster- and strain-specific tags are the
Layer-2 markers. Crucially these are derived from **all** tags of the species' genomes, not
from a pre-built species-unique database: species-unique markers (a genome compared against
genomes of *other* species) are computed for species detection and are orthogonal to
within-species strain structure. Each cluster's database is the union of its member genomes'
single-copy tags; a marker is *unique* to a cluster iff it occurs in exactly one cluster.

**Assembly-quality filtering.** Variable reference completeness biases Jaccard clustering
toward spurious splits: an incomplete genome's marker set is approximately a subset of its
complete twin's, so their Jaccard falls below 1 and they fail to cluster. Because CheckM is
not run in-line, two dependency-free proxies computed from data already at hand are used —
contig count (`--max-contigs`), and single-copy tag count relative to the conspecific median
(`--min-tag-fraction`, a completeness proxy). Genomes far below the median are always flagged;
they are removed only when a threshold is set.

## Layer-2 strain profiling

Sample reads are digested with the database's enzyme set (recorded in the database header) to
give per-marker counts. A cluster is called **present** iff at least *N* of its unique markers
are observed at count ≥ 2 (default *N* = 10, in tag units; the full-k-mer StrainScan floor
of ~1240 k-mers is inappropriate for the ~50–100× sparser tag set). Using only unique markers
makes detection immune to the shared-marker cross-talk that otherwise lets a greedy set-cover
over a large conspecific panel strip shared markers and starve true strains. Each present
cluster's **relative abundance** is estimated from the median sample count over its detected
unique markers — robust to repeat and contamination outliers, and less prone than a joint
regression over shared markers to mis-attributing signal between very similar co-present
strains (a non-negative Elastic Net solver is also provided). Calls are then filtered by a
minimum coverage fraction of their unique markers (`--min-coverage`, default 0.1; suppresses
spurious detection of large, similar clusters whose absolute unique-marker count clears the
floor at a tiny coverage fraction) and a minimum relative abundance (0.02), and renormalised.
When no cluster passes, Strain2bScan reports the species as detectable but not
strain-resolvable with the given enzyme set.

## Multi-species profiling

For community samples, the reads are digested **once** into a shared set of tag counts, which
is then matched against every per-species cluster database in parallel; the per-species
marginal cost is a hash-set lookup rather than a re-count, so the total cost is independent of
the number of species. A **Layer-1 species gate** requires that at least *G* species-specific
markers (markers unique to one species across the panel — the Fast2bRAD-M species layer) be
present before a species is profiled at the strain level (default *G* = 200), suppressing the
cross-species false positives that arise because strain markers are unique only within a
species. In practice the species set to resolve is supplied by an upstream Fast2bRAD-M
species-level profiling step.

## Implementation

Strain2bScan is written in Rust with no third-party dependencies. Data-parallelism (genome
digestion, sketch construction, pairwise clustering, read digestion) uses scoped `std` threads
(`STRAIN2BSCAN_THREADS`; default = all cores). Databases are stored as sparse strain×marker
tables with an inverted unique-marker index and the enzyme set in the header.

## Benchmarking

**Datasets.** (i) A real *C. acnes* benchmark: a 64-genome reference panel (14 ground-truth
strains + 50 background, NCBI accessions pinned) and five paired-end mock samples from
MockMetagenomes4Benchmark (~100k read pairs each, ~12× total). (ii) A simulated multi-species
benchmark: 40 real species × 4 strains (160 NCBI genomes) and 20 samples, each mixing strains
from eight species at log-normal depth ≥1×. (iii) Cross-species mocks for *Staphylococcus
aureus* and *S. epidermidis* (60-genome panels each; 2–5 strains/sample, log-normal ≥1×,
matching the *C. acnes* design). (iv) A reference-degradation gradient in which the truth
strains' database genomes are degraded to completeness 100→50 % (with co-varying contamination
0→10 % and fragmentation), samples held fixed. All simulated reads are error-free 150 bp.

**Comparison to StrainScan.** StrainScan (v1.0) was run on the same *C. acnes* samples with
its low-depth modes for the depth series. StrainScan's DB build requires Linux-only
dependencies (dashing); the same-panel build-and-profile comparison is provided as a Linux
script.

**Metrics.** Detection precision, recall and F1 at a 0.01 presence threshold; abundance error
by L1 distance and Bray–Curtis dissimilarity over the union of predicted and true labels,
evaluated at cluster resolution (ground-truth strains mapped to their clusters). Wall-clock
time and peak resident set size were measured with `/usr/bin/time -l` on a 16-core Apple
silicon machine. All scripts, pinned accession lists, result tables and figure code are in the
Strain2bScan-paper repository; every figure is regenerable with `make figures`.
