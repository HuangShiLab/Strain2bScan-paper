# Discussion

Strain2bScan couples a sparse 2bRAD marker space to a StrainScan-style [1] clustering and unique-marker
framework in a fast Rust implementation. Across the tested settings, 2bRAD tags preserved genome-wide strain
ordering that 16S did not (median Spearman 0.94 versus 0.36), precision was 1.0 in simulated single-species
benchmarks, and detection matched StrainScan at 0.5× coverage. In a run-level 15-species StrainScan rerun,
species-median precision was 1.0 for both tools and the confidence intervals for paired recall/F1
differences included zero, although Strain2bScan had higher species-median point estimates
(0.733/0.844 versus 0.667/0.800). Its consistent advantages were database-construction cost, single-pass
community profiling, and retention under high host DNA (Figs 2, 3, S5, 10, 11).

**Native 2bRAD for low-biomass, high-host microbiomes.** Native reduction occurs before host DNA dominates a
library. In the primary MSA-1002 comparison at 99% human DNA and a 10⁻⁴ abundance threshold, Strain2bScan
retained F1 = 1.0 with Bray–Curtis dissimilarity 0.302; StrainScan [1] had thresholded recall 0.05 and abundance
dissimilarity 0.974, whereas inStrain [8] had F1 = 0.333 (Fig 11). In eight saliva subjects, leave-one-timepoint-out host
identification reached 100% for strain profiles versus 78.1% for species profiles, although PERMANOVA R² was
similar (0.757 versus 0.755). That classifier result was within-subject and involved only eight subjects; it
should not be interpreted as independent strain-identity accuracy. In three usable paired samples, 65/65
shotgun strain-cluster calls were present
in native 2bRAD-M, which produced 128–163 additional candidate strain-cluster calls per sample (Fig 8). These
are concordance and sensitivity observations, not independent validation. Tumour, FFPE and skin applications
remain untested but are natural targets because host content and degradation can limit shotgun sequencing.

**Conventional metagenomes at community scale.** In shotgun mode, Strain2bScan was approximately 8× faster
and 11× lighter in single-panel testing. Its main structural advantage is that one digest serves all species,
so the marginal cost of another species is a hash lookup. In a 55-species analysis, measured Strain2bScan
runtime was 121–146× lower than the corresponding projected per-species StrainScan cost (Fig 9C). The primary
high-host mock supported detection and abundance retention (Fig 11), and saliva calls were contained in the
native-2bRAD call set (Fig 8), but the clean-mock precision trade-off means these results support a shared
framework rather than uniform superiority.

**Cluster resolution is the honest unit of strain analysis.** Short reads cannot separate strains that
share almost all of their sequence, so both StrainScan [1] and Strain2bScan resolve to *clusters* of
near-identical strains and we evaluate at that resolution. The clusters-to-genomes ratio at the 0.95 cut
is a per-species property; for diverse species (*C. acnes*) clusters are essentially single strains,
whereas for near-clonal *M. tuberculosis* (5 clusters from 40 genomes) the cluster is the honest level of
claim. Resolving to clusters does not cost precision (it stayed 1.0 even for low-diversity
*S. epidermidis*); what varies with diversity is recall.

**Design choices.** Occurrence-based uniqueness; a marker is cluster-unique only if absent, at any copy
number, from every other cluster; eliminates false-unique markers from single-copy filtering and keeps
precision 1.0 in similar-strain species. MinHash-sketch [9,10] clustering scales database construction while
producing partitions identical to exact Jaccard and to StrainScan's own pre-built *P. copri* clustering
(112 → 51 clusters; `results/panelsize_prevotella.tsv`). Reference incompleteness is the one factor that
genuinely degrades strain identification under Jaccard; an incomplete genome's markers are a subset of a
complete relative's, so the two fall below the 0.95 similarity cut and split (Fig 4). We address this two
ways: the built-in assembly-quality filter drops low-quality genomes before clustering
(`--min-tag-fraction`/`--max-contigs`), and the optional **`--containment` clustering mode**
(max-containment) keeps subset genomes with their complete relatives, restoring precision/recall down to
~80 % completeness and removing the near-clonal *M. tuberculosis* cluster-fragmentation artifact (Fig 4).
The same completeness concern motivated restricting the Fig 2 motivation panel to complete/near-complete
genomes, where the 2bRAD-over-16S advantage is unchanged (median 0.90→0.94), confirming it is not a
draft-assembly artifact. The enzyme count is an explicit resolution/cost control
(~4 enzymes captures most recall; Fig 5), and single-enzyme BcgI operation is what enables native
2bRAD-M libraries.

**Limitations and future work.** (i) **Near-clonality** caps recall on species where short reads cannot
distinguish strains (*M. tuberculosis*); this limit is intrinsic and shared by all short-read tools;
and where it bites hardest StrainScan failed to complete (>3.3 h, >25 GB) while Strain2bScan finished in
~1 s at precision 1.0 (Fig 10). Layering a within-cluster overlap/regression step on top of
occurrence-based uniqueness is the clearest algorithmic target for raising near-clonal recall.
(ii) **Reference panels must be niche-appropriate and genome-rich** for real communities: a generic
pathogen panel with few genomes per species gave no saliva signal because real strains map uniformly
across arbitrary clusters, whereas a genome-rich oral panel recovered the full individual-discrimination
signal; panel design is a real determinant of strain-level performance on open-world data.
(iii) **Host-limited shotgun** cannot reach the low-abundance strain tail on high-host samples; this is a
property of the input, not the tool, and is precisely the gap the native-2bRAD mode fills.
(iv) **Reference incompleteness is improved but not fully solved.** The `--containment` clustering mode
recovers accuracy when an incomplete genome has a complete relative in the panel (Fig 4), but it cannot
recover markers that are simply *absent* (a strain whose only reference is a partial assembly) nor undo
contamination that injects foreign tags, so it converges with Jaccard below ~70 % completeness; and
because it merges more aggressively it trades a little resolution on complete panels (hence opt-in, not
the default). Making strain identification *more resistant* to incomplete references is a clear direction:
a **completeness-aware detection gate** (scale the unique-marker floor by each genome's estimated
completeness, or gate on a *fraction* of a cluster's available markers rather than an absolute count;
analogous to the Layer-1 breadth term) so incomplete strains are not gated out; a
**best-quality-representative** marker set per cluster (define the cluster's markers from its most-complete
member); **upstream completeness/contamination estimation and decontamination** (CheckM2 [11] / GUNC [12]) feeding
the quality filter; and **pangenome-based imputation** of missing markers from complete conspecifics. The
irreducible case (a strain represented only by a low-completeness, contaminated genome) is a data limit
no clustering can overcome. (v) The public-cohort experiment in Table 4 adds real WGS profiling of longitudinal, multi-site and
low-biomass samples, and Table 5 adds an isolate-derived panel compatibility check. These are application-oriented
subsets rather than closed-world truth benchmarks: read sets used to build the isolate panel were also
used for self-recovery, and the real metagenomes lack exhaustive strain truth. The comparator scope in
Table 8 separates the executed StrainScan [1] and inStrain [8] comparisons from the supplementary StrainGST [2]
rerun; StrainGR variant calling and native-2bRAD input for StrainGE remain untested, and a full StrainGR
benchmark is required before claiming StrainGE-equivalent genomic resolution. Deeper external
validation, oral-cancer case/control analysis, FFPE and degraded material, broader multi-tool
comparison (sylph), and completing the Fast2bRAD-M species layer remain future work.

(vi) The PRJNA1517970 body-site panel still produced no calls. This finding shows that species choice alone does not
solve low-marker input: the negative result reflects sparse strain-marker evidence and incomplete body-site
coverage, with panel completeness and sequencing depth remaining limiting. (vii) Leave-one-isolate-out
*E. coli* tests showed that reads from an absent isolate were assigned to a merged cluster formed by the
two closest available relatives. This behaviour avoids inventing a false isolate-specific call but also
means that absence of the true strain cannot be inferred from a nearest-relative assignment. Future work
therefore requires broader pangenome panels and explicit absent-strain models.

(viii) **Strain-level profiling did not improve case/control prediction over species-level on two
native BcgI 2bRAD-M oral datasets** (ECC caries vs healthy, and Lim_ORPI clean vs unclean denture;
`results/real_data_strain_benchmark.md`). Species-level classifiers outperformed strain-level
(AUROC 0.81 vs 0.59 for ECC, 0.72 vs 0.52 for Lim_ORPI), largely because BcgI single-enzyme marker
density produced too few strain calls (9–18 prevalent clusters) for stable ML, and the phenotypes
appeared driven more by species presence than by strain identity. This is a marker-density limit of
the native BcgI protocol, not a failure of the resolution framework: in-silico multi-enzyme or k-mer
marker sources would be needed to test whether strain-level resolution can outperform species-level
classification on these phenotypes.

**Conclusion.** Reduced-representation 2bRAD markers, combined with a StrainScan-style [1] resolution
framework and a fast Rust implementation, make accurate strain-level profiling practical at a fraction of
the compute and memory of full-k-mer methods. The same tool spans two regimes: native 2bRAD-M for
strain-level analysis of low-biomass, high-host microbiomes, and in-silico-digested shotgun for strain
profiling across communities of many species and many samples.
