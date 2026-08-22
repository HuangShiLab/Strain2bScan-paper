# Strain2bScan — algorithm and computational design (Methods source)

Written against `Strain2bScan` branch `strainscan-port`, HEAD `f404a83`. Every constant and
formula below was read out of the source; file:line references are given so the manuscript text
can be re-verified against the code it describes.

---

## 1. Markers: what a 2bRAD tag becomes

**Digestion.** A type-IIB restriction enzyme cuts on both sides of its recognition site and
releases a short fixed-length fragment — the 2bRAD tag. All 16 enzymes from the Fast2bRAD-M
table are implemented (`enzymes.rs:142-249`), with tag lengths 25–33 bp taken from the `@site`
regexes of `2bRADExtraction.pl`, so extracted tags are identical to the reference Perl and
Fast2bRAD-M implementations. An enzyme is a set of *patterns*; a pattern is a list of literal
anchors at fixed offsets within the tag window, plus, for the three IUPAC-degenerate enzymes
(BaeI, HaeIV, Hin4I), single positions restricted to a base class. Positions with no anchor are
unconstrained but must be A/C/G/T (`is_pure_atcg`, `enzymes.rs:80`), which excludes tags spanning
ambiguity codes.

**Strand handling.** Only the forward strand is scanned. Each enzyme carries forward and reverse
patterns that are exact reverse-complement pairs at its tag length, so one pass finds sites in
either orientation; the three palindromic enzymes (BplI, FalI, AlfI) carry a single self-RC
pattern. Digestion is therefore strand-invariant — `digest(S) = digest(revcomp(S))` — without
scanning both strands and without doubling the marker set. This is asserted for every enzyme by
a regression test (`markers.rs`, `digestion_is_strand_invariant`).

**Marker definition.** A marker is the 64-bit FNV-1a hash of the *canonical* tag, i.e. the
lexicographically smaller of the tag and its reverse complement, upper-cased
(`markers.rs::marker_from_tag`). Hashing makes the representation length-agnostic, which matters
because tag length varies by enzyme; a fixed-width 2-bit packing would need a different code path
per enzyme and could not hold 33 bp in a `u64`. The hash defines on-disk marker identity, so it
is fixed: databases built by any version remain readable.

**Single-copy restriction.** Only tags occurring exactly once in a genome are retained
(`markers.rs::single_copy_markers`). Multi-copy tags inflate depth estimates and blur strain
identity; the same restriction is applied by StrainScan and by Fast2bRAD-M's `remove_redundant`.
This is what makes per-tag depth proportional to genome copy number, and it is the assumption
every abundance quantity below rests on.

---

## 2. Clustering and the Cluster Search Tree

**Clusters.** Genomes are grouped by single-linkage clustering on tag-set Jaccard at a default
threshold τ = 0.95 (`cst.rs:36`), matching StrainScan's `hclsMap_95`. Single linkage at τ is
exactly the connected components of the τ-similarity graph, computed by union–find. Above 96
genomes (`MINHASH_ABOVE`, `cst.rs:159`) the pairwise similarities are estimated from bottom-k
MinHash sketches with k = 2000 (`SKETCH_K`) rather than computed exactly; the partition is
unchanged on the panels tested. A `--containment` mode substitutes max-containment for Jaccard,
for panels of uneven assembly completeness, where an incomplete genome otherwise looks distant
from its complete twin.

**Marker taxonomy.** Within a species each tag is classified by how many genomes and clusters
carry it: *species-core* (all clusters), *shared-partial*, *cluster-specific* (one cluster, ≥2
genomes) and *strain-specific* (one genome). Cluster- and strain-specific tags are the
discriminating markers; they are found by within-species incidence, not by reusing a
species-specific database, which would be dominated by species-core tags.

**Uniqueness is occurrence-based.** A marker counts as unique to a cluster only if it is absent —
*at any copy number* — from every other cluster's genomes (`cst.rs::cluster_db`). The weaker test
(`degree == 1` over single-copy sets) mislabels a tag as unique when it is multi-copy, and hence
filtered, in another cluster, while still being reachable from that cluster's reads.

**The tree.** A strictly binary hierarchy is built above the clusters by agglomerative merging on
a cluster-level similarity matrix, updated by the single-linkage rule (a merged node's similarity
to any other is the maximum of its children's). Each node stores

    K(v) = markers core to v's subtree and absent from every genome outside it

computed by indexing each marker once by its carrier genome set: an internal node takes the
markers whose carrier set equals its genome set, and a leaf takes those whose carriers all lie
inside it. Leaves deliberately take the **union** of their members' markers rather than the
intersection, matching what the flat path scores the same cluster on; the intersection makes the
tree a weaker competitor rather than an extension, and empties entirely for clusters whose
members disagree (on real data, a 5-genome *C. acnes* cluster went from 115 markers to 0).

---

## 3. Detection (Layer-1)

For cluster *j* with discriminating panel `U_j` and observed per-marker counts `c_m`:

    N_j  = |U_j|                                   panel size
    D_j  = |{ m in U_j : c_m >= 1 }|               detected markers
    coverage_j = D_j / N_j                         breadth
    support_j  = |{ m in U_j : c_m >= t_j }|       evidence at the singleton policy

**Singleton policy.** `t_j = 2` when the estimated depth is at or above 3 reads/tag,
`t_j = 1` below it (`min_count_for`, `SINGLETON_SAFE_DEPTH = 3.0`). At high depth a genuine
marker is essentially never seen exactly once, so `c = 1` is dominated by sequencing error; at
low depth the reverse holds — under Poisson(λ) the share of *detected* markers seen exactly once
is λ/(e^λ − 1), which is 78 % at λ = 0.5 — so a fixed `c >= 2` rule discards most of the signal
precisely where signal is scarce. Errors generate essentially random tags, which almost never
coincide with one specific cluster's panel, so admitting singletons there costs little
specificity while the support floor still demands many independent hits.

A cluster is called when `support_j >= 8` (`min_support_markers`) and `coverage_j >= 0.1`
(`min_coverage`). The support floor is set from the arithmetic of the marker space rather than
chosen round: support tracks `N × (1 − e^(−λ))`, and on 2bRAD both factors are small — a
discriminating panel is a few dozen tags (median 53 on a 419-cluster *C. acnes* panel) and a
strain at 5 % of a sample sequenced to ~5× per tag sits at λ ≈ 0.27, where only ~24 % of any
panel is observable, giving ~8 expected observations.

**Depth–breadth consistency.** A cluster is rejected when

    coverage_j / (1 − e^(−depth_j))  <  0.5      (`min_consistency`)

Under Poisson sampling a genuinely present cluster at depth λ must show breadth `1 − e^(−λ)`, so
this ratio is ≈1 for a real cluster at any depth. It is ≈f for a *shadow* — a cluster called
because the strain in the sample happens to carry a fraction f of its distinguishing loci, so
those markers appear at the sample strain's full depth across only f of the panel. No coverage
floor can separate the two: measured on a constructed shadow, the spurious cluster showed breadth
0.350 against 0.392 for a genuinely present cluster at 0.4×, while their depths were 7.68× and
0.44×. Swept synthetically, genuine clusters scored 0.949–1.018 across 0.3×–20× and shadows
0.200/0.300/0.495/0.691/0.897 at f = 0.2/0.3/0.5/0.7/0.9.

---

## 4. Abundance (Layer-2)

**Depth estimator.** Per-cluster depth is the **zero-inclusive** mean count over the whole
discriminating panel, with the top 1 % of non-zero observations winsorized down to the 99th
percentile (`panel_stats`, `TRIM_FRACTION = 100`):

    depth_j = (1/N_j) * sum_{m in U_j} min(c_m, kappa_j)

Both halves matter. Averaging over the whole panel — zeros included — is what keeps the estimate
proportional to true depth; an estimator restricted to *detected* markers (for example their
median) pins a rare cluster near 1 read/tag however rare it is, compressing the ratio between an
abundant and a rare cluster and flattening the composition. Winsorizing rather than discarding,
and taking the fraction of the non-zero observations rather than of the panel, prevents the guard
against collapsed repeats from deleting real signal when few markers are detected.

**Three abundance scopes.** Because single-copy tags are one per genome copy, reads-per-tag
cancels genome size, so depth is a *cell* (taxonomic) fraction and `depth × n_markers` is a *DNA*
fraction. Three columns are reported:

| column | definition | denominator | interpretation |
|---|---|---|---|
| `abundance` | `depth_j / Σ_{k in species} depth_k` | this species | within-species split (**primary**) |
| `global_abundance` | `depth_j / Σ_k depth_k` | clusters this run resolved | community composition, **cell** fraction |
| `sample_fraction` | `depth_j × G_j / Σ_m c_m` | all tag observations | share of the sequencing, **DNA** fraction |

`abundance` and `global_abundance` compose exactly:
`global_abundance_j = species_abundance_{s(j)} × abundance_j`, verified to 1e-6, so an externally
computed species layer (for example Fast2bRAD-M's) can be substituted for the species term.
`sample_fraction` is the only column whose denominator is fixed by the sequencing rather than by
how well profiling went, and is therefore the only one comparable *between* samples; the
unclassified remainder is printed rather than hidden.

**Ported alternatives, and why they are not the default.** StrainScan's Cluster Search Tree
descent (`--layer1 cst`) and its non-negative ElasticNet over the shared-marker design matrix
(`--layer2 enet`) are both implemented and selectable. Neither is default, on measurement:

- `--layer1` defaults to **`auto`**, which reads off the database whether a tree can help at all —
  descend only if some cluster falls below the support floor (the only case ancestor pooling can
  change an outcome) *and* some internal node carries enough markers to pool. On a dense
  conspecific panel the second condition fails: of 542 internal nodes on 543 *C. acnes* genomes,
  373 carry zero group-specific markers, because clustering at τ has already merged anything
  similar enough for a clade to have a distinct core. The decision and both counts are printed.
- `--layer2 enet` is measurably worse here. Scored on abundance accuracy over identical
  detections, Bray–Curtis rose 0.035 → 0.127 and mean absolute relative error 0.183 → 0.794.
  The cause is structural collinearity, not tuning: each cluster carries ~33 100 markers of which
  only 29–115 are unique, so design columns are ~99.7 % identical and the shared rows constrain
  the *sum* of two near-identical clusters while saying almost nothing about the split. Sweeping
  the penalty makes it monotonically worse (BC 0.127 / 0.145 / 0.223 / 0.330 / 0.360 at
  α = 0 / 0.001 / 0.01 / 0.1 / 1.0).

This is a substantive result about the marker space rather than about the implementation: the
mechanisms StrainScan relies on assume a dense k-mer set, and their preconditions do not hold on
a ~1–2 % genomic subsample.

---

## 5. Multi-species architecture

`multi-profile` digests a sample **once** and matches it against every per-species database in
parallel, which is what makes community-scale profiling cheap — a per-species tool must re-read
and re-process the reads for each species.

**Species gate.** Cluster-uniqueness is defined only *within* one species database, so an absent
species can be hit by a present relative's shared tags. Species-specific markers — tags carried
by exactly one species across the whole panel — are derived from the loaded databases, and each
species is placed in one of three tiers from **absolute** marker evidence, never relative
abundance:

    r            = max(1 − e^(−λ_s), 0.25)          reachable fraction, floored
    resolve_gate = max(ceil(floor × r), ceil(frac × total × r), detect, 1)
    detect_gate  = min(detect, resolve_gate)

with `floor = 200`, `detect = 10`, `frac = 0` by default. Scaling by the reachable fraction is
what keeps a fixed 200-marker bar from being unreachable by construction in a low-input or
high-host sample, where only a few percent of any panel is observable at all. The floor is
clamped at 25 % of its configured value (`MIN_FLOOR_FRACTION`) because an unbounded scaling is
self-cancelling — the observed count is itself proportional to r, so `present ≥ floor·r` reduces
to `total ≥ floor` at every depth.

**Cross-species quantification filter.** Detection and depth are restricted to markers specific
to their species across the whole panel. Without it, a tag that merely looks cluster-specific
inside one database but also occurs in a congener's genomes receives that congener's reads. Panels
routinely contain such pairs — *S. aureus*/*S. epidermidis*, three streptococci and two
lactobacilli in ATCC MSA-1002, most oral communities — so this is a systematic abundance error
rather than a rare accident. On a two-congener mock the affected cluster's depth was overstated
**3×** (29.9× against a true 10×); with the filter, 10.2×. The fraction of markers excluded is
printed per run.

---

## 6. Computational design

The crate is **dependency-free** (`std` only), which keeps the build reproducible and the binary
self-contained. Optimizations below were each verified to leave output unchanged.

**Database.** Sparse: per cluster, the set of marker hashes it carries, plus an inverted
`marker → #clusters` degree index. A dense strain × marker matrix, as `strainscan-rust` uses,
reaches tens of GB at real panel sizes.

**Digestion hot path.** Allocation-free end to end. Enzyme scanning is case-insensitive so
sequences are read directly from the input buffer with no upper-cased copy per read;
canonicalization chooses the orientation by comparing the forward strand against its reverse
complement one base at a time and hashes the winner in place, with no `revcomp` buffer; counts
land directly in a hash map with no intermediate vector per sequence.

**Hashing.** Marker keys are `u64`, so the map uses an inlined FxHash rather than the default
SipHash. Marker values are unchanged — FNV-1a of the canonical tag — so only in-memory bucket
assignment differs and existing databases still load.

**I/O.** FASTA/FASTQ are streamed, plain or gzipped (decompression pipes through `gzip -dc`,
preserving the zero-dependency property). Peak memory is one batch rather than the whole file.

**Measured.** On 4 M reads (289 MB) against a two-species panel, wall time fell from 1.71 s to
0.51 s and peak RSS from 282 MB to 11 MB with all 16 enzymes; gzipped input costs no extra
wall-clock. Tree construction was reduced from ~O(n^3.6) to ~O(n^1.6) by updating a cluster-level
similarity matrix instead of re-deriving max-linkage per merge, and by the carrier-set index for
K(v); on 543 genomes the pairwise similarity scan, which dominates a large build, is parallel and
takes the build from 71.4 s to 58.4 s.

**Two optimizations were rejected on measurement**, and are recorded because both are intuitive:
replacing the hash-set Jaccard with a sorted-vector two-pointer intersection is 0.82× — slower,
because a hash lookup on `u64` is cheap and the two-pointer must traverse both arrays — though it
would halve memory; and fusing the multi-enzyme scan into a single pass is 1.43× on a 2.5 Mb
contig but 0.93× on 150 bp reads, so it would have to be dispatched on sequence length rather
than applied globally.

---

## 7. What this design makes possible

Three capabilities follow from operating on 2bRAD tags rather than a full k-mer set, and none is
available to a tool built on the latter:

1. **Native 2bRAD libraries can be strain-profiled at all.** Wet-lab 2bRAD-M produces reads that
   *are* tags. Every existing strain profiler consumes shotgun reads, so a laboratory running
   2bRAD-M has no route to strain resolution — including on the low-biomass, degraded and
   host-dominated samples where 2bRAD-M was adopted precisely because shotgun fails.
2. **A community is digested once, not once per species**, which is what makes hundreds of
   samples against hundreds of species tractable.
3. **Absolute, cross-comparable depth**, and with it an honest unclassified fraction — a
   regression coefficient is meaningful only inside its own design, whereas reads-per-tag on
   single-copy markers is a physical quantity.

The cost is stated plainly: markers are ~1–2 % of the genome, so a cluster's discriminating panel
is tens rather than thousands of tags, and the detection limit is correspondingly set by how many
of those few can be observed. That constraint is what motivates the depth-adaptive gating in §3,
and it is why the two StrainScan mechanisms in §4 do not transfer.
