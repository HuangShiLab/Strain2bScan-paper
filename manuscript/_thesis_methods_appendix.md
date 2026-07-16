The subsections above summarise the pipeline; the following give the algorithmic detail, data
structures and complexity in full, as implemented in the Rust source (`src/`).

### Type-IIB tag model and enzyme patterns

Each type-IIB restriction enzyme recognises a short, partially degenerate motif and cleaves a fixed
distance to either side, releasing a tag of constant length. An enzyme is modelled as a triple
(*upstream gap*, *anchored pattern set*, *downstream gap*): for BcgI the excised fragment is 32 bp with a
central `CGA…TGC`-type recognition anchor and *N*-runs on either flank; the sixteen enzymes of the
Fast2bRAD-M table differ in tag length (32–38 bp) and anchor. Digestion scans every offset of a sequence
and tests the anchor set; because type-IIB sites are palindromically constrained, providing both the
forward and reverse anchor patterns lets a single left-to-right pass over one strand recover the tags that
would be produced from both strands, which halves the work and — importantly — yields exactly one canonical
marker per site rather than the two-fold-inflated set an explicit both-strand scan produces (the earlier
both-strand workaround was retired for this reason, ~3.5× fewer, correct-length markers). For conventional
shotgun input the union of a chosen enzyme set is applied (`--enzyme all` uses all sixteen), enriching the
marker yield ~*n*-fold for *n* enzymes; for native BcgI 2bRAD input the single BcgI pattern is used, because
the reads already are BcgI tags.

### Canonicalisation and hashing

Every tag is reduced to a single 64-bit integer *marker* by (i) taking the lexicographically smaller of
the tag and its reverse complement (the *canonical* form, so a tag and its complement collide), and (ii)
hashing that canonical string with FNV-1a. Genome tags and sample-read tags pass through the identical
canonicalisation and hash, so marker values are internally consistent across the reference database and any
sample; comparison and counting are then exact integer-set operations. Using a 64-bit hash rather than the
raw 2·*L*-bit packed sequence keeps the marker width constant across enzymes (tag lengths differ) and lets
every downstream structure be a hash set or hash map keyed by `u64`; the collision probability at the
marker counts used here (10⁴–10⁶ per species) is negligible.

### Single-copy filtering

For reference genomes only *single-copy* tags — those occurring exactly once in the genome — are retained,
following the practice of StrainScan and Fast2bRAD-M. Multi-copy tags carry copy-number rather than
presence/absence information and would bias both clustering (a repeat expansion inflates set overlap) and
abundance estimation (a repeat contributes disproportionate counts); dropping them makes every retained
marker a clean presence/absence locus. Sample reads are *not* single-copy-filtered — a marker's observed
count in a sample is the quantity of interest — but detection thresholds are stated in tag (marker) units
throughout.

### Within-species clustering: single-linkage, union-find, and MinHash acceleration

Genomes of a species are grouped by their marker sets. Two genomes are joined if their similarity meets a
threshold τ = 0.95 (distance ≤ 0.05); the clusters are the connected components of the resulting graph,
i.e. **single-linkage** clustering, computed with a union-find (disjoint-set) structure in near-linear
time in the number of edges. Single-linkage at threshold τ is exactly the transitive closure of the
"similar-enough" relation, which is the correct semantics here: a chain of pairwise-near-identical genomes
should share a cluster even if its endpoints are not themselves within τ.

The default similarity is the **Jaccard index** of the two marker sets, *J*(A,B) = |A∩B| / |A∪B|. For
panels of ≤ 96 genomes the tool computes exact all-pairs Jaccard directly on the tag sets, at cost
O(*n*²·*m̄*) for *n* genomes of mean marker-set size *m̄*. Above that size it estimates Jaccard from
**bottom-*k* MinHash sketches** (*k* = 2000): each genome's marker set is reduced to its *k* smallest hash
values, and *J* is estimated as the fraction of shared values in the union of the two sketches' bottom-*k*.
This lowers the pairwise cost to O(*n*²·*k*) with *k* ≪ *m̄*, and on the real panels used here yields
partitions *identical* to exact Jaccard (verified on the *P. copri* panel, which reproduces StrainScan's own
112→51 clustering). Sketch construction is O(*n*·*m̄*) once, parallel across genomes.

The resulting clusters are the finest resolution the data support: strains within one 0.95 cluster share
almost all of their markers and cannot be separated from short reads. All accuracy is therefore evaluated
at cluster resolution (ground-truth strains mapped to their clusters), which is the honest unit of claim
for any short-read strain caller.

### Containment clustering for uneven-completeness panels

Jaccard penalises incompleteness. If genome B is an incomplete assembly of the same strain as complete
genome A, its marker set is approximately a subset of A's, so |A∩B| ≈ |B| but |A∪B| ≈ |A|, giving
*J* ≈ |B|/|A| — which falls below τ as soon as B is materially smaller than A, spuriously splitting the two
into different clusters. The consequence is not merely coarser clustering: the shared markers, now present
in *two* clusters, are demoted from *cluster-specific* (discriminating) to *shared-partial*
(non-discriminating), so reads from the complete strain match *both* fragments, injecting false positives
and missing calls (Fig 4A).

The optional `--containment` mode replaces Jaccard with **max-containment**,

&nbsp;&nbsp;&nbsp;&nbsp;*C*(A,B) = |A∩B| / min(|A|,|B|),

which stays ≈ 1 when one marker set is contained in the other and so keeps the incomplete genome clustered
with its complete relative; the merged cluster's marker set is the *union* of its members, so the complete
member supplies the markers the incomplete one lacks and the discriminating tags are preserved (Fig 4A).
This is the containment estimator used by Mash-screen and sourmash for genomes of uneven completeness. It
is exact for small panels; for large panels the intersection is recovered from the sketch-estimated Jaccard
*Ĵ* and the exact set sizes via |A∩B| = *Ĵ*·(|A|+|B|)/(1+*Ĵ*), then divided by min(|A|,|B|). Because
*C* ≥ *J* always, containment merges at least as aggressively as Jaccard and can coarsen clusters on
already-complete panels; it is therefore opt-in (recommended for reference sets of mixed completeness),
with the default remaining Jaccard complemented by the assembly-quality filter.

### Marker classification and the database

Within a species, each tag is labelled by its incidence across the species' clusters: present in **all**
clusters (*species-core*; identifies the species, not a strain), in exactly **one** cluster of ≥ 2 genomes
(*cluster-specific*), in a **single genome** (*strain-specific*), or in **several but not all** clusters
(*shared-partial*). Cluster- and strain-specific tags are the Layer-2 markers; formally, a marker is
*unique* to a cluster iff it occurs in exactly one cluster. Crucially these labels are derived from the
incidence of *all* of the species' single-copy tags, not from any pre-built "species-unique" database:
species-unique markers (a genome compared against genomes of *other* species) are computed separately, for
species detection (Layer-1), and are orthogonal to within-species strain structure.

The database is stored as a **sparse strain × marker table** with the enzyme set in the header and an
**inverted index** from each unique marker to the single cluster that owns it. Profiling therefore reduces
to streaming a sample's markers through the inverted index and incrementing per-cluster counters — the
per-marker work is a single hash lookup.

### Layer-1: which species to strain-profile

Strain markers are unique only *within* a species, so a species that is absent from a sample can be hit
spuriously by shared tags of a present relative. Species selection is therefore made on **absolute
species-specific marker evidence**, never on relative abundance (which conflates community composition with
depth). Let *total* be the number of species-specific markers a species carries (tags unique to that
species across the panel — the same tag space as the Fast2bRAD-M species layer) and *present* the subset
observed in the sample at count ≥ 2. The gate is

&nbsp;&nbsp;&nbsp;&nbsp;*resolve_gate* = max(*G*, ⌈*f* · *total*⌉),&nbsp;&nbsp;&nbsp;*detect_gate* = min(*d*, *resolve_gate*),

with an absolute floor *G* (default 200), a breadth fraction *f* (default 0) that scales the bar to each
species' panel size, and a low detection floor *d* (default 10). This produces three outcomes per species:
**strain-resolved** (*present* ≥ *resolve_gate*; Layer-2 runs), **detected but not strain-resolvable**
(*detect_gate* ≤ *present* < *resolve_gate*; reported at species level with its observed marker breadth,
no strain claim), or **absent**. The middle tier is the honest treatment of a low-abundance species —
present but too faint for a strain claim — rather than a binary drop or an over-call. The breadth term *f*
is scale insurance: on the panels used here *f* = 0 (the shipped default) already gives species precision
1.0, but as a panel grows large enough for a fixed floor to be outrun by cross-species leakage, a small
*f* (≈ 0.02) restores precision by raising the bar in proportion to panel size, where large-panel leakage
concentrates.

### Layer-2: strain detection and abundance

Within each strain-resolved species, a cluster is **called present** iff at least *N* of its unique markers
are observed at count ≥ 2 (default *N* = 10, in tag units — the full-k-mer StrainScan floor of ~1240 k-mers
is inappropriate for the ~50–100× sparser tag set). Using *only* unique markers makes detection immune to
the shared-marker cross-talk that would otherwise let a greedy set-cover over a large conspecific panel
strip shared markers and starve true strains. Each present cluster's **relative abundance** is estimated
from the *median* sample count over its detected unique markers — robust to repeat and contamination
outliers, and less prone than a joint regression over shared markers to mis-attributing signal between very
similar co-present strains; a non-negative Elastic-Net solver over the marker×cluster incidence matrix is
also provided for users who prefer a regression estimate. Calls are then filtered by a minimum coverage
fraction of their unique markers (`--min-coverage`, default 0.1 — suppresses spurious detection of large,
similar clusters whose absolute unique-marker count clears the floor at a tiny coverage fraction) and a
minimum relative abundance (0.02), and renormalised. When no cluster passes, the species is reported as
detectable but not strain-resolvable at the given enzyme set.

### Complexity and the community-scale argument

Let *n* be genomes per species, *m̄* the mean single-copy marker count, *S* species, *N* samples, and *R*
the reads per sample. **Building** a species database costs O(*n*·(genome length)) to digest, O(*n*·*m̄*)
to sketch, and O(*n*²·*k*) to cluster — dominated in practice by digestion, and independent across species
(embarrassingly parallel). **Profiling** one sample against one species costs O(*R*·*L*) to digest the
reads into markers once, plus O(#markers) hash lookups against the inverted index. The decisive point is
the *community* cost. A full-k-mer tool has no shared per-sample representation across species, so it pays
the sample-side count once *per species*: total ≈ *N × S ×* (k-mer count + search). Strain2bScan digests
each sample **once** into a shared marker multiset and matches it against every species' inverted index at a
marginal cost *ε* of a hash-set intersection: total ≈ *N × (digest + S·ε)*. Because *ε* ≪ (k-mer count),
the ratio grows with *S* — an *S*-fold structural advantage confirmed empirically at ~132× on 100 samples
of a 55-species community (Fig 9C) and reproduced in the head-to-head, where StrainScan, lacking a
multi-species mode, must run each community sample once per species and so pays 100–398 s per sample versus
1–9 s for a single Strain2bScan pass (Fig 11F).

### Implementation and determinism

Strain2bScan is written in Rust with **no third-party dependencies**; data-parallelism (genome digestion,
sketch construction, pairwise clustering, read digestion) uses scoped `std` threads
(`STRAIN2BSCAN_THREADS`, default = all cores). All hashing is deterministic (fixed FNV-1a seed) and
clustering is order-independent (union-find over a symmetric edge set), so a given database and reads
produce identical output across runs, thread counts and platforms — verified here by cross-checking the
native arm64 binary against a `linux/amd64` build of the same source, which produced identical cluster
calls on the benchmark samples.
