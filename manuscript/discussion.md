# Strain2bScan — Discussion (draft)

Strain2bScan reframes strain-level profiling around a sparse, reproducible marker set. By porting
StrainScan's clustering-plus-unique-marker framework onto 2bRAD tags in Rust, it turns a
~50–100× smaller marker database into large, structural gains in speed, memory and scalability
**without sacrificing accuracy**: across every species tested it held precision 1.0, matched
StrainScan's detection onset at 0.5× coverage, and matched or exceeded StrainScan's recall on
StrainScan's own databases. The contribution is therefore a strain profiler that is
simultaneously accurate, an order of magnitude faster and lighter per sample, uniquely scalable
across many species, and — because its tags are identical to Fast2bRAD-M / `2bRADExtraction.pl` —
interoperable with a 2bRAD species layer and able to profile native BcgI 2bRAD experimental data.

**Cluster resolution is the honest unit of strain analysis.** Short reads cannot separate strains
that share almost all of their sequence, so both StrainScan and Strain2bScan resolve to *clusters*
of near-identical strains rather than individual genomes, and we evaluate at that resolution. The
clusters-to-genomes ratio at the 0.95 cut is a per-species property: for diverse species such as
*C. acnes* clusters are essentially single strains, whereas for near-clonal species such as
*Mycobacterium tuberculosis* (5 clusters from 40 genomes) the cluster is the appropriate and
honest level of claim. Crucially, resolving to clusters does not cost precision — precision
remained 1.0 even for low-diversity *S. epidermidis*, so the profiler does not over-detect closely
related strains. What varies with diversity is *recall*: where a panel is dominated by
near-identical genomes, fewer distinct clusters can be called. This mirrors StrainScan's own
cluster/CST-level reporting and is, we argue, the correct way to state strain results — to the
resolution the data support.

**Where Strain2bScan wins, and its one honest limit.** Its advantages are per-sample speed and
memory (~8× faster, ~11× lighter than StrainScan) and, above all, scaling in the number of
species: because a sample is digested once and matched against every per-species database, the
marginal cost of an additional species is a hash lookup rather than a re-count. For a community of
*S* species over *N* samples the cost is ≈ *N × (digest + S·ε)* versus ≈ *N × S ×* (k-mer count +
search) for a full k-mer tool run per species — an *S*-fold structural advantage that, measured on
a 55-species community, reached **~132× at 100 samples and ~146× beyond** (minutes versus
projected hours) and grows with community richness. The one genuine limit is recall on near-clonal
species, where the biology itself caps how many strains short reads can distinguish. This limit is
intrinsic and shared by all short-read tools — and where it bites hardest, on *M. tuberculosis*,
StrainScan failed to complete altogether (>3.3 h, >25 GB) while Strain2bScan finished in ~1 s at
precision 1.0, so even at the hard limit Strain2bScan is strictly the more usable of the two.

**Design choices.** Occurrence-based uniqueness — defining a marker as cluster-unique only if it
is absent, at any copy number, from every other cluster's genomes — eliminates false-unique
markers created by single-copy filtering and is what keeps precision at 1.0 in similar-strain
species. MinHash-sketch clustering makes database construction scale to large panels while
producing partitions identical to exact Jaccard, and — as an external check — identical to the
clustering of StrainScan's own pre-built *P. copri* database (112 genomes → 51 clusters).
Assembly-quality filtering counters the completeness bias of Jaccard clustering using tag-count
and contig-count proxies computed from data already in hand. The enzyme count is exposed as an
explicit resolution-versus-cost control: **~4 enzymes** captured essentially all the recall at a
fraction of the markers and compute of the full 16, and single-enzyme BcgI operation is what lets
native BcgI 2bRAD-M libraries be profiled directly.

**Limitations and future work.** Our benchmark reads are simulated, error-free and closed-world —
the truth strain is always present and there is no sequencing error. Closing this gap is the
priority, and the paired shotgun + BcgI-2bRAD saliva chapter is designed to do so on real data:
it tests individual discrimination (strain- versus species-level PERMANOVA R²), open-world
false-positive control, and — through the paired design — validates both real-read performance and
the in-silico-digestion ↔ native-2bRAD equivalence without external ground truth. On the method
side, raising recall on near-clonal species is the clearest algorithmic target: layering a
within-cluster overlap-matrix/regression step, as StrainScan does, on top of occurrence-based
uniqueness could recover additional strains where clusters are near-degenerate. Finally, the
low-input advantage of *wet-lab* 2bRAD libraries — where sequencing is concentrated on the tags,
yielding high per-tag depth at a fraction of the total sequencing of shotgun — is a distinct use
case that Strain2bScan already supports and that warrants dedicated experimental evaluation.

**Conclusion.** Strain2bScan shows that reduced-representation 2bRAD markers, combined with a
StrainScan-style resolution framework and a fast Rust implementation, make accurate strain-level
profiling across many species and many samples practical at a fraction of the compute and memory
of full k-mer methods, while uniquely enabling strain-level analysis of native 2bRAD experimental
data.
