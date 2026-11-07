# Methods

## Overview

Strain2bScan reimplements the two-layer StrainScan strategy [1]; cluster near-identical strains,
then score samples on markers unique to a strain or cluster; but replaces the full k-mer set
with **2bRAD tags**, and implemented in Rust for speed and parallelism. The
pipeline is: (i) digest reference genomes and sample reads into 2bRAD-tag markers; (ii) build,
per species, a within-species cluster database annotated with unique markers; (iii) profile a
sample by detecting present clusters from their unique markers and estimating their abundance;
(iv) at the community level, digest each sample once and match it against every per-species
database, gated by a species-level Layer-1 check.

## 2bRAD tag extraction

Type-IIB restriction enzymes cut on both sides of their recognition site, releasing a
fixed-length fragment (the 2bRAD tag, 25–33 bp depending on enzyme). Each of the 16 enzymes in
the Fast2bRAD-M table [7] is modelled as a set of anchored sequence patterns; literal motifs at
fixed offsets within the tag window, plus, for the three IUPAC-degenerate enzymes (BaeI, HaeIV,
Hin4I), single positions restricted to a base class; unanchored positions are unconstrained but
must be A/C/G/T, which excludes tags spanning ambiguity codes. Scanning every offset and testing
the anchors reproduces the enzyme's digestion sites. Only the forward strand is scanned: each
enzyme carries forward and reverse patterns that are exact reverse-complement pairs at its tag
length (the palindromic enzymes BplI, FalI and AlfI carry a single self-complementary pattern),
so one pass finds sites in either orientation. Digestion is therefore strand-invariant (digest(*S*) = digest(revcomp(*S*))) without scanning both strands and without doubling the
marker set; this is asserted for every enzyme by a regression test. Each tag is canonicalised (the lexicographically smaller of
the tag and its reverse complement) and hashed to a 64-bit integer marker (FNV-1a; genome and
sample tags use the same hash, so marker values are internally consistent). Two input modes
are supported: **BcgI only** for 2bRAD experimental libraries whose reads already are tags, or
a user-chosen enzyme set (`--enzyme all` for all 16) for in-silico digestion of conventional
shotgun reads, which enriches the marker set ~*n*-fold for *n* enzymes. Pooling several enzymes
requires care: ten of the sixteen emit 27 bp tags and their recognition sites genuinely coincide
(every AloI site with a C at offset 16 is also a BsaXI site, and a PpiI site as well if offset 19
is C; on 2 Mb of sequence 64 of 257 AloI sites are all three). A locus is therefore counted once
however many enzymes recognise it, deduplicating on (position, tag length) rather than on the
marker, since two enzymes of different tag lengths at the same offset cut genuinely different
loci. Counting per enzyme instead gives such a locus an apparent copy number of 2–3, and the
single-copy filter then discards it; a systematic 6–7% hole in multi-enzyme panels. For reference genomes
we retain only **single-copy** tags (occurring exactly once in the genome), following
StrainScan's and Fast2bRAD-M's use of single-copy markers for unbiased quantification.

## An alternative marker source: FracMinHash sketching

Everything downstream of extraction consumes only a set of 64-bit marker identities and their
per-sample counts; nothing in clustering, the search tree, detection or abundance refers to
restriction sites. Strain2bScan therefore accepts a second, purely computational marker source
(`--marker-source kmer`): every canonical *k*-mer (default *k* = 31) is hashed with the same
function used for tags, and kept iff its hash falls below 2^64 / *S* for a user-chosen scale *S*,
so approximately 1 / *S* of the *k*-mers survive. This is the FracMinHash rule, and like
restriction digestion (but unlike minimizers) it is *context-free*: whether a locus is selected
depends on that locus alone, never on its neighbours, so marker identity is stable across
genomes. The single-copy restriction, clustering, tree construction and both identification
layers are applied unchanged.

The two sources are not interchangeable in practice and are not offered as competitors. 2bRAD is
a wet-lab protocol that delivers markers directly from a mixed sample without assembly;
FracMinHash is an in-silico rule that needs assembled genomes for the reference side. The sketch
is therefore for shotgun analysis, not for native 2bRAD libraries.

What it buys is **density**. An enzyme panel offers 1, 2, 4 or at most 16 discrete steps; the
sketch scale is continuous. This matters because the Cluster Search Tree's internal nodes are
defined by intersection over a subtree minus the union of everything outside it, and are squeezed
from both sides as a panel grows: a marker absent from any one member of the subtree is lost from
the intersection, and the union to subtract only grows. On a 28-genome panel, where every cluster
comfortably clears the support floor on its own markers, the internal nodes do not:

| marker source | markers per cluster (median) | internal-node group-specific markers (median) | tree usable? |
|---|---|---|---|
| BcgI (1 enzyme) | 650 | 2 | no |
| `recommended` (14) | 11,711 | 32 | marginal |
| `all` (16; the ceiling) | 16,008 | 45 | marginal |
| sketch *S* = 100 | 18,991 | 40 | marginal |
| sketch *S* = 30 | 63,554 | 134 | yes |
| sketch *S* = 10 | 189,986 | 372 | yes |

Two things follow. At matched density the two sources are equivalent (16,008 markers → 45 versus
18,991 → 40), which is the expected result if the tree mathematics is a property of the metric
rather than of restriction enzymes. But the enzyme path has a ceiling; `all` is all sixteen
enzymes and there is nothing above it; and that ceiling sits below the density at which the tree
becomes descendable. The sketch has no ceiling, and reaches it at *S* = 30 for ~3.5× the database
size (76 MB versus 22 MB for `all`, on 28 genomes). This is the mechanism behind the tree being
inert on dense same-species panels, reported in Results.

The sketch was verified to behave as FracMinHash [10] requires. Marker count tracks 1 / *S* to within
sampling noise across three decades of scale (ratios 0.995–1.008 relative to the unsketched set,
on *E. coli* K-12), and the FNV-1a hash produces no collisions on the 989,962 distinct canonical
31-mers of a 1 Mb window, against 0.027 expected for a uniform 64-bit hash. The apparent
shortfall against genome length (4,523,995 markers from 4,641,622 windows) is entirely the
single-copy filter removing the 30,274 multi-copy 31-mers, matching an independent count exactly.

## Reference database construction

**Within-species clustering.** For each species, genomes are grouped by single-linkage
hierarchical clustering at 0.95 marker-set similarity (0.05 distance), matching StrainScan's [1]
`hclsMap_95`. Single-linkage at threshold τ is exactly the connected components of the graph
whose edges join genome pairs with Jaccard ≥ τ, computed with union-find. For panels of ≤96
genomes we use exact all-pairs Jaccard on the tag sets; above that we estimate Jaccard from
bottom-*k* MinHash [9] sketches (*k* = 2000) of each genome's markers, which reduces the pairwise
cost from O(n²·m) to O(n²·k) with *k* ≪ *m* and yields partitions identical to exact on real
data (Results), without the multiple whole-genome alignments used by tools such as SibeliaZ [13]. Clusters are the finest reliable resolution unit: strains within one cluster
are too similar to separate from short reads.

**Containment clustering for uneven-completeness panels (`--containment`).** Jaccard penalises
incompleteness: an incomplete genome's markers are approximately a *subset* of a complete relative's,
so |A∩B|/|A∪B| falls below τ and the two spuriously split. The optional `--containment` mode instead
links on **max-containment**, |A∩B| / min(|A|,|B|), which stays ≈ 1 when one marker set is contained in
the other; the containment estimator used by Mash-screen [14] and sourmash [10] for uneven-completeness genomes.
It is exact for small panels; for large panels the intersection is estimated from the MinHash-sketch
Jaccard and the exact set sizes (|A∩B| = J·(|A|+|B|)/(1+J)), then divided by min(|A|,|B|). Because
max-containment ≥ Jaccard it merges at least as much, so it is opt-in (for reference sets of mixed
completeness) while the default stays Jaccard; the assembly-quality filter below is the
complementary first line of defence.

**Marker classification.** Within a species, each tag is labelled by its within-species
incidence; present in all clusters (*species-core*; detects the species, not strains), in one
cluster with ≥2 genomes (*cluster-specific*), in a single genome (*strain-specific*), or in
several but not all clusters (*shared-partial*). Cluster- and strain-specific tags are the
Layer-2 markers. Crucially these are derived from **all** tags of the species' genomes, not
from a pre-built species-unique database: species-unique markers (a genome compared against
genomes of *other* species) are computed for species detection and are orthogonal to
within-species strain structure. Each cluster's database is the union of its member genomes'
single-copy tags. A marker is *unique* to a cluster iff it is absent (**at any copy number**) from every other cluster's genomes. The weaker test (degree 1 over the single-copy sets alone)
mislabels a tag as unique when it is multi-copy, and therefore filtered, in another cluster while
still being reachable from that cluster's reads.

**Assembly-quality filtering.** Variable reference completeness biases Jaccard clustering
toward spurious splits: an incomplete genome's marker set is approximately a subset of its
complete twin's, so their Jaccard falls below 1 and they fail to cluster. Because CheckM [15] is
not run in-line, two dependency-free proxies computed from data already at hand are used;
contig count (`--max-contigs`), and single-copy tag count relative to the conspecific median
(`--min-tag-fraction`, a completeness proxy). Genomes far below the median are always flagged;
they are removed only when a threshold is set.

## Layer-2: detection and abundance

Sample reads are digested with the database's enzyme set (recorded in the database header) to
give per-marker counts *c*<sub>*m*</sub>. For cluster *j* with discriminating panel *U*<sub>*j*</sub>:

&nbsp;&nbsp;&nbsp;&nbsp; *N*<sub>*j*</sub> = |*U*<sub>*j*</sub>| &nbsp;(panel size), &nbsp;
*D*<sub>*j*</sub> = |{*m* ∈ *U*<sub>*j*</sub> : *c*<sub>*m*</sub> ≥ 1}|, &nbsp;
coverage<sub>*j*</sub> = *D*<sub>*j*</sub> / *N*<sub>*j*</sub>

**Depth-adaptive singleton policy.** Evidence is counted at a threshold *t*<sub>*j*</sub> that
depends on the estimated depth: *t* = 2 at or above 3 reads/tag, *t* = 1 below it. At high depth a
genuine marker is essentially never observed exactly once, so *c* = 1 is dominated by sequencing
error and is filtered, as in StrainScan. At low depth the reverse holds; under Poisson(λ) the
share of *detected* markers seen exactly once is λ/(e<sup>λ</sup> − 1), 78 % at λ = 0.5; so a fixed
*c* ≥ 2 rule discards most of the signal precisely where signal is scarce. Sequencing errors
generate essentially random tags, which almost never coincide with one specific cluster's panel,
so admitting singletons there costs little specificity while the support floor still requires many
independent hits on that one panel.

**Detection.** A cluster is called present when support<sub>*j*</sub> = |{*m* ∈ *U*<sub>*j*</sub> :
*c*<sub>*m*</sub> ≥ *t*<sub>*j*</sub>}| ≥ 8 (`--min-support`) and coverage<sub>*j*</sub> ≥ 0.1
(`--min-coverage`). The support floor follows from the arithmetic of the marker space rather than
being chosen round: support tracks *N* · (1 − e<sup>−λ</sup>), and on 2bRAD both factors are small;
a discriminating panel is a few dozen tags (median 53 across a 419-cluster *C. acnes* panel), and a
strain at 5 % of a sample sequenced to ~5× per tag sits at λ ≈ 0.27, where only ~24 % of any panel
is observable, giving ~8 expected observations.

**Depth–breadth consistency.** A cluster is rejected when

&nbsp;&nbsp;&nbsp;&nbsp; coverage<sub>*j*</sub> / (1 − e<sup>−depth<sub>*j*</sub></sup>) &lt; 0.5 &nbsp;(`--min-consistency`)

Under Poisson sampling a genuinely present cluster at depth λ must show breadth 1 − e<sup>−λ</sup>,
so this ratio is ≈ 1 for a real cluster at any depth. It is ≈ *f* for a **shadow**. a cluster
called because the strain in the sample happens to carry a fraction *f* of its distinguishing loci,
so those markers appear at the sample strain's full depth across only *f* of the panel. No coverage
floor can separate the two, because a shadow and a genuinely rare strain have the same breadth and
differ only in depth: on a constructed shadow the spurious cluster showed breadth 0.350 against
0.392 for a genuinely present cluster at 0.4×, while their depths were 7.68× and 0.44×. Swept
synthetically, genuine clusters scored 0.949–1.018 across 0.3×–20× and shadows
0.200/0.300/0.495/0.691/0.897 at *f* = 0.2/0.3/0.5/0.7/0.9.

**Abundance.** Each called cluster's depth is the **zero-inclusive** mean count over its whole
discriminating panel, with the top 1 % of non-zero observations winsorized to the 99th percentile:

&nbsp;&nbsp;&nbsp;&nbsp; depth<sub>*j*</sub> = (1/*N*<sub>*j*</sub>) Σ<sub>*m* ∈ *U*<sub>*j*</sub></sub> min(*c*<sub>*m*</sub>, κ<sub>*j*</sub>)

Both halves are load-bearing. Averaging over the whole panel (zeros included) is what keeps the
estimate proportional to true depth; an estimator restricted to *detected* markers (for example
their median) pins a rare cluster near 1 read/tag however rare it is, compressing the ratio between
an abundant and a rare cluster and flattening the whole composition. Winsorizing rather than
discarding, and taking the fraction of the *non-zero* observations rather than of the panel,
prevents the guard against collapsed repeats from deleting real signal when few markers are
detected. Because single-copy tags are one per genome copy, reads-per-tag cancels genome size, so
depth is proportional to cell (taxonomic) abundance and depth × *G*<sub>*j*</sub>; where
*G*<sub>*j*</sub> is the cluster's tag count; is proportional to DNA mass.

**Three abundance scopes** are reported, because per-species fractions cannot be concatenated into
a community composition:

| column | definition | denominator | interpretation |
|---|---|---|---|
| `abundance` | depth<sub>*j*</sub> / Σ<sub>*k* ∈ species</sub> depth<sub>*k*</sub> | this species | within-species split (primary) |
| `global_abundance` | depth<sub>*j*</sub> / Σ<sub>*k*</sub> depth<sub>*k*</sub> | clusters this run resolved | community composition, **cell** fraction |
| `sample_fraction` | depth<sub>*j*</sub> · *G*<sub>*j*</sub> / Σ<sub>*m*</sub> *c*<sub>*m*</sub> | all tag observations | share of the sequencing, **DNA** fraction |

The first two compose exactly; global_abundance<sub>*j*</sub> = species_abundance<sub>*s*(*j*)</sub>
× abundance<sub>*j*</sub>; so an externally computed species layer (for example Fast2bRAD-M's) can
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
hash-set lookup rather than a re-count. The marginal work of an additional species is one compact database
lookup/search, while total cost remains weakly dependent on panel content; throughput was approximately flat
over the tested 20–55-species range.

**Which species to strain-profile; the Layer-1 gate.** Strain markers are unique only *within* a
species, so a species absent from a sample can be spuriously hit by a present relative's shared
tags. Strain2bScan therefore decides per species from **absolute species-specific marker
evidence**, never relative abundance (which conflates community composition with sequencing
depth). Let *total* be the species-specific markers a species carries; tags unique to a single
species across the panel, the same tag space as the Fast2bRAD-M species layer; and *present* the
subset observed in the sample at count ≥ 2. The gate is

&nbsp;&nbsp;&nbsp;&nbsp; *r* = max(1 − e<sup>−λ<sub>*s*</sub></sup>, 0.25) &nbsp;(the reachable fraction of the panel)

&nbsp;&nbsp;&nbsp;&nbsp; *resolve_gate* = max(⌈*G*·*r*⌉, ⌈*f*·*total*·*r*⌉, *d*, 1), &nbsp; *detect_gate* = min(*d*, *resolve_gate*)

where λ<sub>*s*</sub> is the species' estimated per-tag depth, taken as the zero-inclusive mean
count over its species-specific markers; a quantity that does not presuppose the species passed
any gate, so the rule is not circular. Scaling by *r* is what keeps a fixed 200-marker bar from
being **unreachable by construction** in a low-input or high-host sample: at 0.05× depth only ~5 %
of any panel is observable, so an unscaled floor files a genuinely present species as absent
however clean the data is. The scaling is clamped at 25 % of the configured floor because an
unbounded version is self-cancelling; the observed count is itself proportional to *r*, so
*present* ≥ *G*·*r* reduces to *total* ≥ *G* at every depth, leaving *d* as the only real
threshold.

with an absolute floor *G* (`--min-species-markers`, default 200), a breadth fraction *f*
(`--min-species-marker-frac`, default 0) that scales the bar to each species' panel size, and a
low detection floor *d* (`--min-species-detect`, default 10). This yields three outcomes per
species: **strain-resolved** (*present* ≥ *resolve_gate*; Layer-2 runs), **detected but not
strain-resolvable** (*detect_gate* ≤ *present* < *resolve_gate*; reported at species level with its
observed marker breadth, no strain claim), or **absent**. The middle tier is the honest treatment
of a low-abundance species (present but too faint to support strain calls) rather than a binary
drop or an over-call. All inputs are computed by Strain2bScan from its own databases and a single
digest of the reads, so the gate needs no external abundance input; for open-world samples the
species presence call can instead be taken from an upstream Fast2bRAD-M step whose species
database is far broader than the strain panel.

**Gate calibration.** On the 55-species panel across normal and low (median 0.62×) depth, the
default floor gives species precision 1.0 at both depths, with leakage species correctly held in
the middle tier; at this panel the breadth term only trades recall, so *f* = 0 is optimal and is
the shipped default. The breadth term is scale insurance: when the floor is relaxed; or the panel
grows large enough for a fixed floor to be outrun by leakage; a small *f* (≈0.02) restores
precision to 1.0 at negligible recall cost, because it raises the bar in proportion to panel size,
where large-panel leakage concentrates (Results; `docs/gate_calibration.md`).

## Cross-species restriction on quantification

Cluster-uniqueness is defined only *within* one species database, so a tag can be unique to a
cluster there and still occur in a congener's genomes. When that congener is co-present, its reads
land on the tag and inflate the cluster's depth. Panels routinely contain such pairs;
*S. aureus*/*S. epidermidis*, three streptococci and two lactobacilli in ATCC MSA-1002, and most
oral communities; so this is a systematic abundance error rather than a rare accident.

Under `multi-profile`, detection and depth are therefore restricted to markers that are specific
to their species **across the whole panel**, using the same species-degree index the Layer-1 gate
is built from. On a two-congener mock the affected cluster's depth was overstated 3× (29.9×
against a true 10×); with the restriction it is 10.2×. The proportion of markers excluded is
reported per run, and `--no-cross-species-filter` disables it for comparison. Single-species
`profile` carries no such information and cannot apply it; a database on its own knows nothing
about the rest of the panel.

## Ported StrainScan layers, and why neither is the default

Both stages of StrainScan's [1] resolution framework are implemented and selectable, so the
architectural choice can be tested rather than asserted. Neither is default, on measurement.

**Layer-1; Cluster Search Tree (`--layer1 cst`).** A strictly binary hierarchy is built above the
clusters; each node stores the markers core to its subtree and absent from every genome outside it.
The descent prunes a whole subtree on one test and, at a leaf, pools the markers of every ancestor
whose sibling branch was never entered; which is what lets a leaf with too few markers of its own
be called at all. `--layer1` defaults to **`auto`**, which reads off the database whether a tree can
help *here*: descend only if some cluster falls below the support floor (the only case pooling can
change an outcome) **and** some internal node carries enough markers to pool. On a dense
conspecific panel the second condition fails; of 542 internal nodes on 543 *C. acnes* genomes, 373
carry zero group-specific markers, because clustering at τ has already merged anything similar
enough for a clade to have a distinct core. Whether a tree helps is a property of the panel, not of
the software, so the decision is made per database and printed with the counts behind it.

**Layer-2; joint non-negative ElasticNet (`--layer2 enet`).** A design matrix over the *shared*
markers, which the unique-only estimator discards, fitted jointly across co-present clusters. It
can in principle resolve a cluster whose tag set is contained in a relative's; one with no unique
markers at all, invisible to the flat path. Measured on identical detections it is worse:
Bray–Curtis 0.035 → 0.127 and mean absolute relative error 0.183 → 0.794. The cause is structural
collinearity rather than tuning; each cluster carries ~33 100 markers of which only 29–115 are
unique, so design columns are ~99.7 % identical and the shared rows constrain the *sum* of two
near-identical clusters while saying almost nothing about the split. Penalising makes it
monotonically worse (Bray–Curtis 0.127 / 0.145 / 0.223 / 0.330 / 0.360 at α = 0 / 0.001 / 0.01 /
0.1 / 1.0).

These are results about the marker space, not about the implementation: both mechanisms assume a
dense k-mer set, and their preconditions do not hold on a ~1–2 % genomic subsample. The few unique
markers carry the split directly; the many shared ones do not.

The sketch marker source lets that claim be tested rather than argued, by sweeping density with
the panel, the clustering and both layers held fixed. It holds for Layer-1: internal-node marker
counts rise monotonically with density and the tree crosses from unusable to usable between
sketch scales 100 and 30 (median 40 → 134 group-specific markers), a density the enzyme path
cannot reach at all. The Layer-2 result was not re-measured at higher density and no claim is
made about it here; the collinearity that defeats the ElasticNet is a property of how much of
each cluster's marker set is *shared*, which denser sampling need not change.

## Implementation

Strain2bScan is written in Rust with no required third-party runtime dependencies. Data-parallelism (genome digestion, sketch
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
rather than the default SipHash; marker *values* are unchanged (FNV-1a of the canonical tag), so
only in-memory bucket assignment differs and databases remain readable across versions.

**I/O.** FASTA and FASTQ are streamed, plain or gzipped, with decompression piped through `gzip`
so the zero-dependency property is preserved. Peak memory is one batch rather than the whole file,
which is what makes multi-GB samples tractable.

## Benchmarking

### Overview

All simulated reads were generated with ART [18] using `art_illumina -p -l 250 -m 600 -s 150`; this
250 bp paired-end configuration includes the ART Illumina quality-error model and was used for every
primary simulation. Reference panels were restricted to complete or near-complete genomes unless stated
otherwise. Detection metrics use cluster-set precision, recall and F1 in each tool's own reference or
call space; abundance error uses Bray–Curtis dissimilarity or L2 distance against sequence-abundance
truth. Timing used the resource manager or `/usr/bin/time` and reports wall-clock time and maximum
resident set size.

### 2bRAD versus 16S reference panels

For the 15-species motivation analysis, ~50 genomes per species were obtained from NCBI through ENA and
restricted to CheckM [15] completeness ≥ 97%, contamination ≤ 5% and assembly level Complete
Genome/Chromosome. Whole-genome distance used bottom-3000 canonical 21-mer MinHash; 2bRAD distance used
BcgI tags built by Strain2bScan; 16S distance used the longest barrnap [20] gene, HMMER [21] alignment
and 21-mer Jaccard. Mash [9] distances were converted with
D(J) = −ln(2J/(1+J)), and per-species Spearman correlations against whole-genome distance used 95% CIs
from 500 genome subsamples.

### Native and shotgun ATCC DNA mock communities

The four ATCC whole-cell mocks were MSA-1002 (20 even strains), MSA-1003 (20 staggered strains),
MSA-1005 and MSA-1007 (six strains each). A combined reference tree contained 28 species and up to six
genomes per species (164 genomes), including the ATCC strain and up to five high-quality conspecific
decoys per species; genomes were clustered at 0.95 similarity with `--containment`. Native BcgI used the
BcgI tree; shotgun used the all-enzyme tree. Strain2bScan used default gates (`--min-support 8`,
`--min-coverage 0.1`, `--min-abundance 0`, `--min-consistency 0.5`). Shotgun StrainScan 1.0.14 [1] used
per-species databases in a `linux/amd64` container. inStrain 1.9.1 [8] used Bowtie2 [23] alignments,
sorted BAM files and skani [22]-dereplicated references at 98% ANI. Genome-level calls required breadth
≥ 0.5 and reference popANI ≥ 0.99999; abundance was normalized and detection thresholded at 10⁻⁴.
Supplementary Fig S4 shows the non-dereplicated control. Mock species-to-ATCC-strain mappings are in
Supplementary Table S14.

### Real saliva time-course

Native BcgI 2bRAD and paired shotgun WMS saliva from PRJNA1131785 comprised eight subjects and four
within-day timepoints. Libraries were profiled against a 19-species oral-commensal panel (up to 25
genomes/species) with species gates of 20 panel-specific markers and two detected species markers and call
gates of two supporting markers and 0.01 coverage. Strain- and species-level cross-species matrices used
`sample_fraction`; within-species analyses used within-species relative abundance. Bray–Curtis
dissimilarities were tested by PERMANOVA with 4,999 label permutations (seed 42) [25]. Host assignment used
leave-one-timepoint-out 1-nearest-neighbour classification, and subject-aware cluster bootstrap with
20,000 resamples provided uncertainty. Full-depth paired comparisons used all 32 WMS/native pairs.

### Simulation pool and comparator workflows

A fixed 15-species pool containing 15–50 complete or near-complete NCBI genomes per species was used for
single-species and multi-species simulations. Single-species samples contained 2, 3 or 5 co-present strains
from the same or different 0.95 clusters at per-strain depths of 0.5, 1, 3, 5 and 10×, with five replicates
per cell (2,025 samples). Multi-species samples mixed approximately 18 species across three community depth
gradients (60 samples). Truth tables recorded genome accession and 0.95-cluster assignment. Supplementary
Table S15 enumerates the single-species design.

Strain2bScan databases were built with `--enzyme all --similarity 0.95` (16 enzymes) and samples were
profiled with `profile` or `multi-profile --enzyme all` using default depth-adaptive gates. StrainScan
v1.0.14 [1] clustered the same genome pool independently and was run in a `linux/amd64` Docker container.
Because clustering differs, each tool was scored against truth genomes mapped to its own clusters. The
matched subset contained 225 intended single-species entries; 204 produced StrainScan final reports and 21
were scored as absent final reports in the failure-aware sensitivity. Multi-species analyses used 12
samples, and absent StrainScan species reports were treated as no detection.

StrainGST 1.3.9 [2], the reference-search component of StrainGE, used default k = 23, StrainGE's
near-subset removal and 0.90 Jaccard reference clustering. One species-level pan-genome database was built
for each panel. Each sample was k-merized once and searched against every relevant species database; calls
used StrainGST score ≥ 0.02. Accuracy was mapped through StrainGE's 0.90-reference clusters. Wall time
summed sample k-merization and all species searches; peak RSS was the maximum across those stages.

inStrain 1.11.0 [8] aligned reads to full species genome pools with Bowtie2 and scored genomeCoverage
detections at coverage ≥ 0.1 and breadth ≥ 0.05 after mapping references to truth 0.95 clusters. Profiling
runtime used the same matched simulation set.

### FracMinHash marker mode

The optional k-mer marker mode was benchmarked with k = 31 and FracMinHash scale 100. Species databases
were clustered from the same genome pool at 0.95 similarity, and all 225 matched single-species samples were
profiled with `--marker-source kmer`. Calls were scored against truth strains mapped to the k-mer-mode
0.95 clusters.

### Computational profiling

The 225-sample matched simulation set was used for runtime and memory summaries. Strain2bScan/StrainScan
same-container comparisons report median wall time and peak RSS by depth for the 190 completed matched pairs.
StrainGST reports 225 samples as the sum of sample k-merization and species searches. inStrain reports 225
quick-profile runs from sorted BAM files; Bowtie2 alignment was not included in this stage. Database-build cost
was summarized across the 15 species panels. Multi-species throughput used 55 real species panels and the
projection described in Results.

### Metrics and mock scoring

For primary mock comparisons, precision, recall and F1 use abundance ≥ 10⁻⁴; AUPR summarizes the
threshold-free abundance ranking. Abundance error is Bray–Curtis dissimilarity or L2 distance, and
corresponding similarities are reported as one minus dissimilarity. All mock metrics are evaluated at
0.95-cluster resolution for Strain2bScan/StrainScan and 98%-ANI dereplicated genome references for inStrain.
StrainGST uses its own 0.90-reference space. Metrics were therefore not pooled across tools.

### Rhesus macaque GI spatial cohort

The final real-cohort application used an author-generated rhesus macaque GI spatial cohort of 148 libraries:
128 GI samples (66 luminal and 62 mucosal) and 20 body controls. Four additional vendor throat/trachea
files were excluded. Clean native BcgI libraries contained 32-bp tags. We selected 18 prevalent gut species
and 315 study MAGs from the cohort assembly catalogue, clustered genomes within species at 0.95 Jaccard
similarity, and built per-species BcgI databases representing 303 strain-cluster units. Cross-species marker
filtering retained 29,162 of 30,008 markers. All 148 libraries were profiled with `--enzyme BcgI
--min-species-markers 20 --min-species-detect 2 --min-support 2 --min-coverage 0.01
--min-consistency 0.5`. Sample-level, species-level and cluster-level summaries were calculated for GI,
luminal, mucosal and body-control groups. LI-enriched clusters required peak LI prevalence >= 0.3 and
region specificity >= 0.8.

For mucosal/luminal tests, luminal and mucosal samples were matched within monkey, GI region and GI
position. Missing strain-cluster values were represented as zero, and within-species global-abundance
compositions were compared with two-sided Wilcoxon signed-rank tests followed by Bonferroni correction.
The 66 luminal libraries were also profiled from processed shotgun reads in the same BcgI marker space.
Species and strain-call concordance used Jaccard distances; dominant clusters were compared only where both
modalities detected the species.

For metabolic coupling, targeted bile-acid SRM chromatograms were extracted from mzML files, peak areas
were normalized to internal standards, and negative-mode precursor classes were grouped by precursor mass
into LCA-like, dihydroxy, cholate-like and conjugated/high-mass categories. HUMAnN `pathabundance-cpm`
tables were aggregated and screened for butyrate/SCFA, propionate, lactate, mixed-acid, sialic-acid,
N-acetyl-sugar, fucose/rhamnose and glucarate pathways. In LI samples, strain composition and BA or
functional composition were residualized against monkey identity and LI position before Spearman
correlation; Benjamini–Hochberg FDR was estimated separately for BA and functional tests. Because BA names
were not available from the vendor quantification table, BA features were interpreted at precursor-class
rather than individual-compound resolution. Associations are ecological and do not assign pathway genes to
individual strain clusters.

### Strain-sharing rates

Pairwise strain sharing was computed only among samples with at least one species-level call and used the
same gate-filtered 0.95 cluster calls as the preceding analyses. For each sample pair, species detected in
both samples formed the denominator, following the strain-sharing-rate framework used in epidemiological
microbiome studies [29,30]. A species contributed to the primary dominant-cluster rate when its
highest-within-species-abundance cluster label was identical in both samples. A species contributed to the
secondary any-cluster rate when at least one cluster label was shared; this secondary metric is more
sensitive to multi-strain coexistence but less conservative. Exact merged labels (for example,
`C0|C1|C2`) were considered shared only when the identical label occurred in both samples.

Saliva comparisons comprised 496 pairs from 32 native BcgI libraries (eight subjects and four within-day
timepoints) and contrasted same-subject between-timepoint pairs with between-subject pairs. Rhesus GI
comparisons comprised 8,128 pairs from 128 spatial libraries. Groups distinguished within-monkey
adjacent-position, within-monkey mucosal/luminal, between-monkey same-age same-site and old/young
corresponding-site comparisons. O1/O2 were 20-year-old animals and Y1/Y2 were 8-year-old animals. Differences
between prespecified groups used 10,000 pair-level bootstrap resamples with seed 20261105. Because samples
are nested within subjects or monkeys, these intervals describe pairwise-rate variability rather than
independent-host inference. Sharing of a 0.95 marker-cluster label is not, without SNV or popANI
validation, evidence of strain identity at nucleotide resolution.

### Public real-metagenome cohorts and isolate-derived panels

To supplement the controlled analyses, we profiled representative paired-end WGS subsets from four
public BioProjects: PRJNA288562 (pregnancy saliva/vagina/distal-gut time series), PRJNA1517970 (preterm
vaginal and meconium metagenomes), PRJNA1191223 (preterm-infant stool time series) and PRJNA1191225
(preterm-infant isolate WGS). The 22 run accessions, SRA spot/base counts, mean read lengths, instruments, sample roles and generic-panel outcomes are enumerated in `results/realworld_cohort_benchmark/public_run_manifest.tsv`. These were exploratory application subsets selected to test multi-site
stability, low-biomass specificity and isolate-derived strain recovery; they were not complete
epidemiological cohorts. Reads were supplied directly to Strain2bScan as paired FASTQ without an external host-removal or trimming step; the tool's own streaming I/O and marker-count gates were used as described above. For the generic-panel screen we used the existing 20-species MSA database in
`all`-enzyme mode with `--min-species-markers 20 --min-species-detect 2 --min-support 2
--min-coverage 0.01 --min-abundance 0`. PRJNA288562 analyses used subject T23 at gestational days 84 and
273 for all three body sites. PRJNA1517970 analyses used three vaginal, three meconium and one negative
extraction-control library.

For cohort-specific compatibility testing, six PRJNA1191225 isolate read sets were assembled with SPAdes 4.3.0 [28]
(`--isolate`, 8 threads). Three *E. coli* assemblies were clustered at 0.95 similarity into three
resolvable units; one assembly each from *B. longum*, *B. breve* and *B. bifidum* was built as a
single-genome species database. Panels used `--enzyme recommended`. All six isolate read sets were then
profiled with the cohort panel using `--min-species-markers 50 --min-species-detect 3 --min-support 2
--min-coverage 0.01 --min-abundance 0`. Because these reads also generated the assemblies, the 6/6
self-recovery is a panel/read compatibility check, not an independent classification benchmark. P08
W1–W3 metagenomes were profiled with the same cohort panel using `--min-species-markers 20
--min-species-detect 2` and otherwise the gates above. Timings were measured on a local Apple Silicon
workstation; full command lines and per-sample outputs are provided in
`results/realworld_cohort_benchmark/`.
To avoid reference mismatch in PRJNA1517970, we built a body-site-oriented panel from 13 species commonly
representing vaginal, neonatal-gut and meconium communities: *Bifidobacterium adolescentis*,
*Cutibacterium acnes*, *Enterococcus faecalis*, *Escherichia coli*, *Gardnerella vaginalis*,
*Lactobacillus gasseri*, *Lactobacillus jensenii*, *Prevotella bivia*, *Staphylococcus aureus*,
*Staphylococcus epidermidis*, *Streptococcus agalactiae*, *Streptococcus mitis* and *Streptococcus
mutans*. We used up to six RefSeq or ATCC genomes per species (74 genomes across the 13 species; two *S. mitis* genomes were
available). Databases were built with `--enzyme recommended` and profiled with the same thresholds as the
generic-panel screen. To test behaviour when a true strain is absent from the panel, we also built three
leave-one-isolate-out *E. coli* databases, each containing only the other two *E. coli* assemblies plus
the three Bifidobacterium controls, and profiled the held-out isolate reads.
