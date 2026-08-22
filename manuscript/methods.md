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
fixed-length fragment (the 2bRAD tag, 25–33 bp depending on enzyme). Each of the 16 enzymes in
the Fast2bRAD-M table is modelled as a set of anchored sequence patterns — literal motifs at
fixed offsets within the tag window, plus, for the three IUPAC-degenerate enzymes (BaeI, HaeIV,
Hin4I), single positions restricted to a base class; unanchored positions are unconstrained but
must be A/C/G/T, which excludes tags spanning ambiguity codes. Scanning every offset and testing
the anchors reproduces the enzyme's digestion sites. Only the forward strand is scanned: each
enzyme carries forward and reverse patterns that are exact reverse-complement pairs at its tag
length (the palindromic enzymes BplI, FalI and AlfI carry a single self-complementary pattern),
so one pass finds sites in either orientation. Digestion is therefore strand-invariant —
digest(*S*) = digest(revcomp(*S*)) — without scanning both strands and without doubling the
marker set; this is asserted for every enzyme by a regression test. Each tag is canonicalised (the lexicographically smaller of
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

**Containment clustering for uneven-completeness panels (`--containment`).** Jaccard penalises
incompleteness: an incomplete genome's markers are approximately a *subset* of a complete relative's,
so |A∩B|/|A∪B| falls below τ and the two spuriously split. The optional `--containment` mode instead
links on **max-containment**, |A∩B| / min(|A|,|B|), which stays ≈ 1 when one marker set is contained in
the other — the containment estimator used by Mash-screen and sourmash for uneven-completeness genomes.
It is exact for small panels; for large panels the intersection is estimated from the MinHash-sketch
Jaccard and the exact set sizes (|A∩B| = J·(|A|+|B|)/(1+J)), then divided by min(|A|,|B|). Because
max-containment ≥ Jaccard it merges at least as much, so it is opt-in (for reference sets of mixed
completeness) while the default stays Jaccard; the assembly-quality filter below is the
complementary first line of defence.

**Marker classification.** Within a species, each tag is labelled by its within-species
incidence — present in all clusters (*species-core*; detects the species, not strains), in one
cluster with ≥2 genomes (*cluster-specific*), in a single genome (*strain-specific*), or in
several but not all clusters (*shared-partial*). Cluster- and strain-specific tags are the
Layer-2 markers. Crucially these are derived from **all** tags of the species' genomes, not
from a pre-built species-unique database: species-unique markers (a genome compared against
genomes of *other* species) are computed for species detection and are orthogonal to
within-species strain structure. Each cluster's database is the union of its member genomes'
single-copy tags. A marker is *unique* to a cluster iff it is absent — **at any copy number** —
from every other cluster's genomes. The weaker test (degree 1 over the single-copy sets alone)
mislabels a tag as unique when it is multi-copy, and therefore filtered, in another cluster while
still being reachable from that cluster's reads.

**Assembly-quality filtering.** Variable reference completeness biases Jaccard clustering
toward spurious splits: an incomplete genome's marker set is approximately a subset of its
complete twin's, so their Jaccard falls below 1 and they fail to cluster. Because CheckM is
not run in-line, two dependency-free proxies computed from data already at hand are used —
contig count (`--max-contigs`), and single-copy tag count relative to the conspecific median
(`--min-tag-fraction`, a completeness proxy). Genomes far below the median are always flagged;
they are removed only when a threshold is set.

## Layer-2: detection and abundance

Sample reads are digested with the database's enzyme set (recorded in the database header) to
give per-marker counts *c*<sub>*m*</sub>. For cluster *j* with discriminating panel *U*<sub>*j*</sub>:

&nbsp;&nbsp;&nbsp;&nbsp;*N*<sub>*j*</sub> = |*U*<sub>*j*</sub>| &nbsp;(panel size), &nbsp;
*D*<sub>*j*</sub> = |{*m* ∈ *U*<sub>*j*</sub> : *c*<sub>*m*</sub> ≥ 1}|, &nbsp;
coverage<sub>*j*</sub> = *D*<sub>*j*</sub> / *N*<sub>*j*</sub>

**Depth-adaptive singleton policy.** Evidence is counted at a threshold *t*<sub>*j*</sub> that
depends on the estimated depth: *t* = 2 at or above 3 reads/tag, *t* = 1 below it. At high depth a
genuine marker is essentially never observed exactly once, so *c* = 1 is dominated by sequencing
error and is filtered, as in StrainScan. At low depth the reverse holds — under Poisson(λ) the
share of *detected* markers seen exactly once is λ/(e<sup>λ</sup> − 1), 78 % at λ = 0.5 — so a fixed
*c* ≥ 2 rule discards most of the signal precisely where signal is scarce. Sequencing errors
generate essentially random tags, which almost never coincide with one specific cluster's panel,
so admitting singletons there costs little specificity while the support floor still requires many
independent hits on that one panel.

**Detection.** A cluster is called present when support<sub>*j*</sub> = |{*m* ∈ *U*<sub>*j*</sub> :
*c*<sub>*m*</sub> ≥ *t*<sub>*j*</sub>}| ≥ 8 (`--min-support`) and coverage<sub>*j*</sub> ≥ 0.1
(`--min-coverage`). The support floor follows from the arithmetic of the marker space rather than
being chosen round: support tracks *N* · (1 − e<sup>−λ</sup>), and on 2bRAD both factors are small —
a discriminating panel is a few dozen tags (median 53 across a 419-cluster *C. acnes* panel), and a
strain at 5 % of a sample sequenced to ~5× per tag sits at λ ≈ 0.27, where only ~24 % of any panel
is observable, giving ~8 expected observations.

**Depth–breadth consistency.** A cluster is rejected when

&nbsp;&nbsp;&nbsp;&nbsp;coverage<sub>*j*</sub> / (1 − e<sup>−depth<sub>*j*</sub></sup>) &lt; 0.5 &nbsp;(`--min-consistency`)

Under Poisson sampling a genuinely present cluster at depth λ must show breadth 1 − e<sup>−λ</sup>,
so this ratio is ≈ 1 for a real cluster at any depth. It is ≈ *f* for a **shadow** — a cluster
called because the strain in the sample happens to carry a fraction *f* of its distinguishing loci,
so those markers appear at the sample strain's full depth across only *f* of the panel. No coverage
floor can separate the two, because a shadow and a genuinely rare strain have the same breadth and
differ only in depth: on a constructed shadow the spurious cluster showed breadth 0.350 against
0.392 for a genuinely present cluster at 0.4×, while their depths were 7.68× and 0.44×. Swept
synthetically, genuine clusters scored 0.949–1.018 across 0.3×–20× and shadows
0.200/0.300/0.495/0.691/0.897 at *f* = 0.2/0.3/0.5/0.7/0.9.

**Abundance.** Each called cluster's depth is the **zero-inclusive** mean count over its whole
discriminating panel, with the top 1 % of non-zero observations winsorized to the 99th percentile:

&nbsp;&nbsp;&nbsp;&nbsp;depth<sub>*j*</sub> = (1/*N*<sub>*j*</sub>) Σ<sub>*m* ∈ *U*<sub>*j*</sub></sub> min(*c*<sub>*m*</sub>, κ<sub>*j*</sub>)

Both halves are load-bearing. Averaging over the whole panel — zeros included — is what keeps the
estimate proportional to true depth; an estimator restricted to *detected* markers (for example
their median) pins a rare cluster near 1 read/tag however rare it is, compressing the ratio between
an abundant and a rare cluster and flattening the whole composition. Winsorizing rather than
discarding, and taking the fraction of the *non-zero* observations rather than of the panel,
prevents the guard against collapsed repeats from deleting real signal when few markers are
detected. Because single-copy tags are one per genome copy, reads-per-tag cancels genome size, so
depth is proportional to cell (taxonomic) abundance and depth × *G*<sub>*j*</sub> — where
*G*<sub>*j*</sub> is the cluster's tag count — is proportional to DNA mass.

**Three abundance scopes** are reported, because per-species fractions cannot be concatenated into
a community composition:

| column | definition | denominator | interpretation |
|---|---|---|---|
| `abundance` | depth<sub>*j*</sub> / Σ<sub>*k* ∈ species</sub> depth<sub>*k*</sub> | this species | within-species split (primary) |
| `global_abundance` | depth<sub>*j*</sub> / Σ<sub>*k*</sub> depth<sub>*k*</sub> | clusters this run resolved | community composition, **cell** fraction |
| `sample_fraction` | depth<sub>*j*</sub> · *G*<sub>*j*</sub> / Σ<sub>*m*</sub> *c*<sub>*m*</sub> | all tag observations | share of the sequencing, **DNA** fraction |

The first two compose exactly — global_abundance<sub>*j*</sub> = species_abundance<sub>*s*(*j*)</sub>
× abundance<sub>*j*</sub> — so an externally computed species layer (for example Fast2bRAD-M's) can
be substituted for the species term. `sample_fraction` is the only column whose denominator is
fixed by the sequencing rather than by how well profiling went, and is therefore the only one
comparable *between* samples; the unclassified remainder is reported rather than hidden. Ground
truth stated as taxonomic abundance (as in the ATCC MSA standards) should be compared against
`global_abundance`, and truth stated as genomic DNA against `sample_fraction`; the two differ by
genome size, which spans >6× within a single mock community.

When no cluster passes, Strain2bScan reports the species as detectable but not strain-resolvable
with the given enzyme set.

## Multi-species profiling and species selection

For community samples, the reads are digested **once** into a shared set of tag counts, matched
against every per-species cluster database in parallel; the per-species marginal cost is a
hash-set lookup rather than a re-count, so the total cost is independent of the number of species.

**Which species to strain-profile — the Layer-1 gate.** Strain markers are unique only *within* a
species, so a species absent from a sample can be spuriously hit by a present relative's shared
tags. Strain2bScan therefore decides per species from **absolute species-specific marker
evidence**, never relative abundance (which conflates community composition with sequencing
depth). Let *total* be the species-specific markers a species carries — tags unique to a single
species across the panel, the same tag space as the Fast2bRAD-M species layer — and *present* the
subset observed in the sample at count ≥ 2. The gate is

&nbsp;&nbsp;&nbsp;&nbsp;*r* = max(1 − e<sup>−λ<sub>*s*</sub></sup>, 0.25) &nbsp;(the reachable fraction of the panel)

&nbsp;&nbsp;&nbsp;&nbsp;*resolve_gate* = max(⌈*G*·*r*⌉, ⌈*f*·*total*·*r*⌉, *d*, 1), &nbsp; *detect_gate* = min(*d*, *resolve_gate*)

where λ<sub>*s*</sub> is the species' estimated per-tag depth, taken as the zero-inclusive mean
count over its species-specific markers — a quantity that does not presuppose the species passed
any gate, so the rule is not circular. Scaling by *r* is what keeps a fixed 200-marker bar from
being **unreachable by construction** in a low-input or high-host sample: at 0.05× depth only ~5 %
of any panel is observable, so an unscaled floor files a genuinely present species as absent
however clean the data is. The scaling is clamped at 25 % of the configured floor because an
unbounded version is self-cancelling — the observed count is itself proportional to *r*, so
*present* ≥ *G*·*r* reduces to *total* ≥ *G* at every depth, leaving *d* as the only real
threshold.

with an absolute floor *G* (`--min-species-markers`, default 200), a breadth fraction *f*
(`--min-species-marker-frac`, default 0) that scales the bar to each species' panel size, and a
low detection floor *d* (`--min-species-detect`, default 10). This yields three outcomes per
species: **strain-resolved** (*present* ≥ *resolve_gate*; Layer-2 runs), **detected but not
strain-resolvable** (*detect_gate* ≤ *present* < *resolve_gate*; reported at species level with its
observed marker breadth, no strain claim), or **absent**. The middle tier is the honest treatment
of a low-abundance species — present but too faint to support strain calls — rather than a binary
drop or an over-call. All inputs are computed by Strain2bScan from its own databases and a single
digest of the reads, so the gate needs no external abundance input; for open-world samples the
species presence call can instead be taken from an upstream Fast2bRAD-M step whose species
database is far broader than the strain panel.

**Gate calibration.** On the 55-species panel across normal and low (median 0.62×) depth, the
default floor gives species precision 1.0 at both depths, with leakage species correctly held in
the middle tier; at this panel the breadth term only trades recall, so *f* = 0 is optimal and is
the shipped default. The breadth term is scale insurance: when the floor is relaxed — or the panel
grows large enough for a fixed floor to be outrun by leakage — a small *f* (≈0.02) restores
precision to 1.0 at negligible recall cost, because it raises the bar in proportion to panel size,
where large-panel leakage concentrates (Results; `docs/gate_calibration.md`).

## Cross-species restriction on quantification

Cluster-uniqueness is defined only *within* one species database, so a tag can be unique to a
cluster there and still occur in a congener's genomes. When that congener is co-present, its reads
land on the tag and inflate the cluster's depth. Panels routinely contain such pairs —
*S. aureus*/*S. epidermidis*, three streptococci and two lactobacilli in ATCC MSA-1002, and most
oral communities — so this is a systematic abundance error rather than a rare accident.

Under `multi-profile`, detection and depth are therefore restricted to markers that are specific
to their species **across the whole panel**, using the same species-degree index the Layer-1 gate
is built from. On a two-congener mock the affected cluster's depth was overstated 3× (29.9×
against a true 10×); with the restriction it is 10.2×. The proportion of markers excluded is
reported per run, and `--no-cross-species-filter` disables it for comparison. Single-species
`profile` carries no such information and cannot apply it — a database on its own knows nothing
about the rest of the panel.

## Ported StrainScan layers, and why neither is the default

Both stages of StrainScan's resolution framework are implemented and selectable, so the
architectural choice can be tested rather than asserted. Neither is default, on measurement.

**Layer-1 — Cluster Search Tree (`--layer1 cst`).** A strictly binary hierarchy is built above the
clusters; each node stores the markers core to its subtree and absent from every genome outside it.
The descent prunes a whole subtree on one test and, at a leaf, pools the markers of every ancestor
whose sibling branch was never entered — which is what lets a leaf with too few markers of its own
be called at all. `--layer1` defaults to **`auto`**, which reads off the database whether a tree can
help *here*: descend only if some cluster falls below the support floor (the only case pooling can
change an outcome) **and** some internal node carries enough markers to pool. On a dense
conspecific panel the second condition fails — of 542 internal nodes on 543 *C. acnes* genomes, 373
carry zero group-specific markers, because clustering at τ has already merged anything similar
enough for a clade to have a distinct core. Whether a tree helps is a property of the panel, not of
the software, so the decision is made per database and printed with the counts behind it.

**Layer-2 — joint non-negative ElasticNet (`--layer2 enet`).** A design matrix over the *shared*
markers, which the unique-only estimator discards, fitted jointly across co-present clusters. It
can in principle resolve a cluster whose tag set is contained in a relative's — one with no unique
markers at all, invisible to the flat path. Measured on identical detections it is worse:
Bray–Curtis 0.035 → 0.127 and mean absolute relative error 0.183 → 0.794. The cause is structural
collinearity rather than tuning — each cluster carries ~33 100 markers of which only 29–115 are
unique, so design columns are ~99.7 % identical and the shared rows constrain the *sum* of two
near-identical clusters while saying almost nothing about the split. Penalising makes it
monotonically worse (Bray–Curtis 0.127 / 0.145 / 0.223 / 0.330 / 0.360 at α = 0 / 0.001 / 0.01 /
0.1 / 1.0).

These are results about the marker space, not about the implementation: both mechanisms assume a
dense k-mer set, and their preconditions do not hold on a ~1–2 % genomic subsample. The few unique
markers carry the split directly; the many shared ones do not.

## Implementation

Strain2bScan is written in Rust with **no third-party dependencies**, which keeps the build
reproducible and the binary self-contained. Data-parallelism (genome digestion, sketch
construction, the pairwise similarity scan, read digestion) uses scoped `std` threads
(`STRAIN2BSCAN_THREADS`; default = all cores).

**Database.** Sparse: per cluster, the set of marker hashes it carries, plus an inverted
marker → cluster-degree index and the enzyme set in the header. A dense strain × marker matrix
reaches tens of GB at real panel sizes.

**Digestion hot path.** Allocation-free end to end. Enzyme scanning is case-insensitive, so
sequences are read directly from the input buffer with no upper-cased copy per read (this also
handles soft-masked reference genomes correctly); canonicalisation chooses the orientation by
comparing the forward strand against its reverse complement one base at a time and hashes the
winner in place, with no reverse-complement buffer; counts land directly in a hash map with no
intermediate vector per sequence. Because marker keys are `u64`, the maps use an inlined FxHash
rather than the default SipHash — marker *values* are unchanged (FNV-1a of the canonical tag), so
only in-memory bucket assignment differs and databases remain readable across versions.

**I/O.** FASTA and FASTQ are streamed, plain or gzipped, with decompression piped through `gzip`
so the zero-dependency property is preserved. Peak memory is one batch rather than the whole file,
which is what makes multi-GB samples tractable.

**Measured.** On 4 M reads (289 MB) against a two-species panel with all 16 enzymes, wall time is
0.51 s and peak resident memory 11 MB (1.71 s / 282 MB before these optimisations); gzipped input
costs no additional wall-clock. Tree construction scales as ~O(*n*<sup>1.6</sup>) after replacing a
per-merge recomputation of max-linkage with an incrementally updated cluster-level similarity
matrix, and a per-node set-subtraction with a single carrier-set index; on 543 genomes the
pairwise similarity scan that dominates a large build is parallel, taking the build from 71.4 s to
58.4 s. Every optimisation was verified to leave output unchanged — for the tree, node for node
against the serial build at *n* = 543.

Two intuitive optimisations were **rejected on measurement** and are recorded here because both
are commonly assumed to help: replacing the hash-set Jaccard with a sorted-vector two-pointer
intersection is 0.82× — slower, because a hash lookup on `u64` is cheap while the two-pointer must
traverse both arrays — although it would halve the memory; and fusing the multi-enzyme scan into a
single pass is 1.43× on a 2.5 Mb contig but 0.93× on 150 bp reads, so it would have to be
dispatched on sequence length rather than applied globally.

## Benchmarking

**Datasets.** (i) A real *C. acnes* benchmark: a 64-genome reference panel (14 ground-truth
strains + 50 background, NCBI accessions pinned) and five paired-end mock samples from
MockMetagenomes4Benchmark (~100k read pairs each, ~12× total). (ii) A simulated multi-species
benchmark: 55 real species × ~4 strains (218 NCBI genomes) and 30 samples, each mixing strains
from twelve species at log-normal depth ≥1× (plus a low-depth variant, median 0.62×, used for
gate calibration). (iii) Cross-species mocks for *Staphylococcus
aureus* and *S. epidermidis* (60-genome panels each; 2–5 strains/sample, log-normal ≥1×,
matching the *C. acnes* design). (iv) A reference-degradation gradient in which the truth
strains' database genomes are degraded to completeness 100→50 % (with co-varying contamination
0→10 % and fragmentation), samples held fixed. All simulated reads are error-free 150 bp.

**Real-data and motivation datasets.** (v) *2bRAD-vs-16S motivation* (Fig 2): 15
pathogenic/commensal species, ~50 genomes each from NCBI accession lists (ENA FASTA), **restricted to
complete/near-complete assemblies** (CheckM completeness ≥ 97 %, contamination ≤ 5 %, assembly level
Complete Genome/Chromosome; `data/genome_qc_16s_panel.tsv`). Between-strain distance was computed in
three spaces — whole-genome (bottom-3000 canonical 21-mer MinHash), 2bRAD (Strain2bScan `build` BcgI
tags) and 16S (longest gene per genome via barrnap 0.9 + HMMER, 21-mer Jaccard) — all with the Mash
transform D(J) = −ln(2J/(1+J)); per species the 2bRAD and 16S pairwise vectors were correlated (Spearman)
against the whole-genome vector, with 95 % CIs from 500 genome subsamples. (vi) *ATCC DNA mocks,
strain-level (Fig 6, Fig 12, Fig S3, Fig S4)*: four whole-cell mocks — MSA-1002 (20 strains,
even; native BcgI 2bRAD and shotgun WMS across a 0/90/95/99/99.9 % human-DNA ladder and a 1→0.001 ng
low-biomass ladder, SRA PRJNA1131785), MSA-1003 (20 strains, staggered), MSA-1005 and MSA-1007 (6 strains
each). A single unified combined tree was built from **28 species × up to 6 genomes = 164 genomes** (each
mock species = its ATCC genome + up to 5 high-quality conspecific decoys, CheckM completeness ≥ 90 %,
contamination ≤ 5 %, within-species ANI 95–99.9 % to the ATCC reference by skani), clustered at 0.95
similarity with `--containment`; native 2bRAD used the BcgI tree and shotgun used the all-enzyme tree.
Strain2bScan was run with `--min-abundance 0 --min-coverage 0.2`. On the shotgun samples it was compared
against **StrainScan** 1.0.14 (per-species databases, `linux/amd64` container) and **inStrain** 1.10.0
(bowtie2 → `inStrain profile` against a 98 %-ANI dereplicated reference; the non-dereplicated reference is
shown as a control in Fig S4). Each tool was scored in its own 0.95-similarity cluster space
against the mock ground truth (`Ground_truth/*`, sequence abundance), reporting precision, recall, F1,
AUPR (abundance-threshold sweep, Ye et al. 2019), and Bray–Curtis and L2 similarity to the truth profile
(2bRAD-M, 2021); scorer `scripts/score_all.py`, figures `scripts/plot_figs_h.py`. (vii) *Real saliva* (Fig 7, Fig 8): native BcgI
2bRAD (and paired shotgun WMS) saliva from PRJNA1131785, 8 subjects × 4 within-day timepoints, profiled
against a 19-species oral-commensal panel (up to 25 genomes/species). Strain- and species-level relative
abundances → Bray–Curtis → PERMANOVA (adonis, subject/timepoint factors) and leave-one-timepoint-out
1-NN host classification; shotgun R1 (in-silico BcgI) compared to native 2bRAD calls per sample. Full
per-dataset procedures and accessions are in `docs/` (`motivation_16s.md`,
`saliva_individual_discrimination.md`, `saliva_temporal_ml.md`, `saliva_concordance.md`).

**Systematic head-to-head on a 15-species simulated benchmark (Fig 11, Table 1–3).** A common
benchmark was built from a fixed pool of 15 pathogenic/commensal species (15–50 complete/near-complete
NCBI genomes each; `figure_raw_data/sim_pool_manifest.tsv`). *Single-species* samples were generated for
every species as 2/3/5 co-present strains drawn either from the same or from different 0.95 clusters, at
per-strain coverages 0.5/1/3/5/10× with uneven abundance ratios (following StrainScan's simulation
design), 5 replicates per cell — 2 025 samples. *Multi-species* samples mixed ~18 co-present species
(one to a few strains each) across three community depth gradients — 60 samples. Reads were simulated
with ART (`art_illumina`, 150 bp paired-end) from the truth genomes; truth tables record each strain's
species, genome accession and 0.95-cluster assignment.

Both tools built their databases from the **same genome pool** and profiled the **same reads**.
Strain2bScan databases were built with `cluster --enzyme all --similarity 0.95` and profiled with
`profile` / `multi-profile --enzyme all` (reads decompressed, R1+R2 concatenated). StrainScan (v1.0.14,
bioconda) is Linux-x86-only — it ships `dashing_s128` and `jellyfish-linux` ELF binaries, a Python-3.7
`.so`, and an R reclustering step — so it was run inside a Docker `linux/amd64` container (QEMU emulation
on Apple Silicon; `strainscan_build`, then `strainscan -i R1 -j R2 -d DB`). Because the two tools cluster
genomes independently, **each tool was scored in its own cluster space**: predicted clusters were compared
against the truth strains mapped into that tool's clusters — for Strain2bScan via the truth `cluster`
column, for StrainScan via its `Cluster_Result/hclsMap_95_recls.txt` (report `Cluster_ID` = `C`+cluster
id) — and precision/recall/F1 computed over the cluster sets per sample. Strain2bScan profiled all 2 025 +
60 samples; StrainScan profiled a matched subset (different-cluster mixtures, k = 2/3/5, one replicate, all
depths; near-clonal *M. tuberculosis* via its same-cluster samples) — 204 depth-matched paired
single-species samples across 14 species, plus 4 multi-species samples per depth. StrainScan has no
multi-species mode, so each community sample was profiled once per species database and the per-sample
cost taken as the sum of wall-clock over species (peak RSS as the maximum).

*Timing.* Strain2bScan build and profile times are native (arm64). To compare profiling speed free of the
emulation confound, a `linux/amd64` Strain2bScan binary (zero-dependency `cargo build`) was run **inside
the same container** on the same subset, giving the same-environment ratio of Fig 11E/Table 2; StrainScan
build times are reported in the emulated environment (an upper bound). DB build for *K. pneumoniae*
(47 genomes × 5.5 Mb) did not complete under StrainScan (killed at a 100-min cap; still in the k-mer-matrix
step past 1 h 40 min on a longer retry), and that species is omitted from the paired accuracy set. The
container's VM memory was raised to 56 GB because StrainScan's build peaks at ~28 GB (vs ≤0.4 GB for
Strain2bScan). Scripts: `scripts/plot_sim_headtohead.py` and the drivers under `scratchpad/eval/`
(`run_s2b_{single,multi}.py`, `run_strainscan_{single,multi}.py`, `run_s2b_emulated_single.py`,
`analyze_headtohead.py`); raw per-sample tables in `figure_raw_data/sim_headtohead/`.

**Comparison to StrainScan (curated-DB and per-sample benchmarks).** In addition to the common benchmark
above, StrainScan (v1.0) was run on its **own** reference databases (Fig 10) and on the same *C. acnes*
per-sample profiling comparison (Fig 9A), using its low-depth modes for the depth series.

**Metrics.** Detection precision, recall and F1 at a 0.01 presence threshold; abundance error
by L1 distance and Bray–Curtis dissimilarity over the union of predicted and true labels,
evaluated at cluster resolution (ground-truth strains mapped to their clusters). Wall-clock
time and peak resident set size were measured with `/usr/bin/time -l` on a 16-core Apple
silicon machine. All scripts, pinned accession lists, result tables and figure code are in the
Strain2bScan-paper repository; every figure is regenerable with `make figures`.
