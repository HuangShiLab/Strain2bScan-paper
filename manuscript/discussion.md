# Strain2bScan — Discussion (draft)

Strain2bScan reframes strain-level profiling around a sparse, reproducible marker set. By
porting StrainScan's clustering-plus-unique-marker framework onto 2bRAD tags in Rust, it turns
a ~50–100× smaller marker database into large, structural gains in speed, memory and
scalability, while retaining accurate strain resolution wherever the underlying biology and
sequencing permit. The contribution is therefore not higher accuracy than full k-mer methods,
but a different point on the accuracy–efficiency frontier — one that makes community-wide,
population-scale strain surveillance tractable, and that uniquely accommodates BcgI 2bRAD
experimental data.

**Cluster resolution is the honest unit of strain analysis.** Short reads cannot separate
strains that share almost all of their sequence, so both StrainScan and Strain2bScan resolve
to *clusters* of near-identical strains rather than individual genomes, and we evaluate at that
resolution. We make this explicit and quantitative: the clusters-to-genomes ratio at the 0.95
cut is a per-species *resolvability* score, and Strain2bScan reports when a species is
detectable but not strain-resolvable with the given enzyme set. For diverse species such as
*C. acnes*, clusters are essentially single strains and precision is high; for low-diversity
species such as *S. epidermidis*, the cluster is the appropriate and honest level of claim, and
we report cluster-level detections and abundances there. This mirrors StrainScan's own
cluster/CST-level reporting and, we argue, is the correct way to state strain results: to the
resolution the data support, with a per-species measure of what that resolution is.

**Where Strain2bScan wins, and where it does not.** Its advantages are per-sample speed and
memory and, above all, scaling in the number of species: because a sample is digested once and
matched against every per-species database, the marginal cost of an additional species is a
hash lookup rather than a re-count. For a community of *S* species over *N* samples the cost is
≈ *N × (digest + S·ε)* versus ≈ *N × S ×* (k-mer count + search) for a full k-mer tool run
per species — an *S*-fold structural advantage that, on our 40-species panel, projects to
~200× at 100 samples (minutes vs hours) and scales further with community richness. The costs
of this design are sensitivity at very low per-strain depth (<~1×, where the sparse marker set
lacks statistical power) and precision for near-identical strains in low-diversity species
(where cluster over-detection persists at ~0.5–0.6 precision even after occurrence-based
uniqueness). Both are quantified rather than hidden: the depth sweep defines the usable
coverage range, and the resolvability score flags the hard species.

**Design choices.** Occurrence-based uniqueness — defining a marker as cluster-unique only if
it is absent, at any copy number, from every other cluster's genomes — was the single most
effective accuracy lever for similar-strain species, because it eliminates false-unique markers
created by single-copy filtering. MinHash-sketch clustering makes database construction scale
to large panels with partitions identical to exact Jaccard. Assembly-quality filtering counters
the completeness bias of Jaccard clustering, using tag-count and contig-count proxies computed
from data already in hand. The enzyme count is exposed as an explicit efficiency-versus-
resolution control, with ~8 enzymes near-optimal for the species tested.

**Limitations and future work.** Our simulated reads are error-free; an error model (and
held-out truth strains) would harden the benchmarks. The same-panel head-to-head against
StrainScan required Linux-only build dependencies and is provided as a script rather than run
here. The residual over-detection on near-identical strains is the classical hard case that
StrainScan addresses with a within-cluster overlap-matrix/lasso Layer-2; layering that on top
of occurrence-based uniqueness is the clearest path to higher precision in low-diversity
species. Finally, the low-input advantage of *wet-lab* 2bRAD libraries — where sequencing is
concentrated on tags, inverting the low-depth disadvantage seen for in-silico digestion of
shotgun data — is a distinct use case that Strain2bScan already supports and that warrants a
dedicated experimental evaluation.

**Conclusion.** Strain2bScan shows that reduced-representation 2bRAD markers, combined with a
StrainScan-style resolution framework and a fast Rust implementation, make strain-level
profiling across many species and many samples practical at a fraction of the compute and
memory of full k-mer methods, with an operating envelope that is measured and reported rather
than assumed.
