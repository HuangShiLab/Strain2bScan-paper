**Positioning within the strain-typing landscape.** Strain2bScan is best understood as a
*representation* change to the reference-k-mer family rather than a new inference principle: it keeps
StrainScan's two-layer logic — cluster near-identical strains, then score samples on markers unique to a
cluster — but swaps the full k-mer set for the 2bRAD tag set, which is 50–100× sparser yet still
genome-wide and single-copy. That single change is what buys the two-to-three orders of magnitude in build
and profiling cost measured here (Fig 9–11), and it is orthogonal to the inference logic, so the accuracy
of the framework is preserved (equal precision, matched-or-higher recall). Relative to the other families,
the trade-offs are explicit. *Marker-gene SNP methods* (StrainPhlAn) reconstruct within-clade haplotypes
and can in principle exceed cluster resolution, but require sufficient depth over a fixed marker set and do
not natively serve the low-biomass or the native-2bRAD regime. *Sketch/containment estimators* (sylph)
are extremely fast at the species/genome-containment level but are not designed to resolve co-present
strains within a species at the cluster level. *Reference-k-mer methods* (StrainScan, StrainGE) sit
closest to Strain2bScan in what they claim, and the head-to-head against StrainScan on a common benchmark
is therefore the sharpest test: same genomes, same reads, each tool in its own cluster space. Strain2bScan
matched StrainScan's precision, exceeded its recall, and completed a species (*K. pneumoniae*) StrainScan
could not build — while being dramatically cheaper. A fuller comparison against StrainGE, sylph and
StrainPhlAn on identical inputs is the natural next benchmark; the framework and evaluation harness built
here are designed to accommodate it.

**Interpreting the head-to-head, and its caveats.** Three points bound the interpretation of Fig 11 and
Table 1–3. First, StrainScan is Linux-x86-only and was run under emulation on the Apple-silicon test
machine; its *wall-clock* is therefore an upper bound, which is why the profiling comparison is reported in
a same-environment ratio (an emulated `linux/amd64` build of Strain2bScan) rather than native-vs-emulated,
and why the build-time comparison is framed as an order-of-magnitude structural difference rather than a
precise multiplier. Even after discounting emulation by a generous ~10×, the build-cost gap (minutes-to-
tens-of-minutes and 8–28 GB, versus seconds and ≤0.4 GB) and the memory gap remain large and are structural
— StrainScan's peak memory is set by the full-k-mer matrix, not by emulation. Second, each tool is scored
in *its own* cluster space; because the two tools cluster the same genomes slightly differently (e.g. 27 vs
30 clusters for *E. coli*), the precision/recall values are not paired at the level of individual clusters
but at the level of "was each truth strain's cluster detected" — the only fair comparison when two tools
define clusters independently, and the honest resolution unit for short reads. Third, the recall advantage
is concentrated at shallow depth and on minor co-present strains: at 10× both tools saturate, so the
practical benefit is the ability to recover the low-coverage tail at a given sequencing budget, consistent
with the low-biomass results of Part I. The one species StrainScan could not build under emulation
(*K. pneumoniae*, 47 genomes × 5.5 Mb) is an extreme illustration of the same build-cost problem rather
than a separate failure mode.

**Why cluster resolution is the right unit — and where it can be pushed.** Both tools resolve to clusters
of near-identical strains because short reads cannot separate genomes that share almost all of their
sequence; reporting a specific strain when the data support only a cluster would be over-claiming. The
cluster-to-genome ratio is a per-species property of genuine diversity: for diverse species (*C. acnes*)
clusters are essentially single strains, whereas for near-clonal *M. tuberculosis* a handful of clusters
is the honest ceiling, and this chapter confirms that resolving to clusters costs no precision even in
low-diversity species. The clearest algorithmic route to *sub-cluster* resolution — the main limitation
shared with all short-read callers — is to layer a within-cluster step on top of the occurrence-based
detection used here: for a called cluster, a per-locus overlap or non-negative regression over the SNP-
bearing tags of its members could apportion signal among within-cluster genomes when depth allows, without
disturbing the between-cluster uniqueness logic that guarantees precision.

**Reference incompleteness: solved in part, and the remaining frontier.** The chapter shows that reference
incompleteness — not the sparsity of the tag set — is the one factor that genuinely degrades strain
identification, because an incomplete genome's markers are a subset of a complete relative's and Jaccard
splits them. The `--containment` mode addresses the *clustering* half of this cleanly (Fig 4): it keeps
subset genomes with their complete relatives and restores precision and recall to near-baseline down to
~80–90 % completeness, and it removes the near-clonal *M. tuberculosis* fragmentation artifact. What
containment cannot do is recover markers that are simply *absent* — a strain represented only by a partial
assembly — nor undo contamination that injects foreign tags; below ~70 % completeness it therefore
converges with Jaccard, and because it merges more aggressively it trades a little resolution on complete
panels (hence opt-in). Making strain identification *more resistant* to incomplete references is a concrete
programme: (i) a **completeness-aware detection gate** that scales the unique-marker floor by each genome's
estimated completeness (or gates on a *fraction* of a cluster's available markers rather than an absolute
count), so genuinely incomplete strains are not gated out; (ii) a **best-quality-representative** marker
set per cluster, defining the cluster's markers from its most-complete member; (iii) **upstream
completeness/contamination estimation and decontamination** (CheckM2, GUNC) feeding the quality filter; and
(iv) **pangenome-based imputation** of missing markers from complete conspecifics. The irreducible case — a
strain whose only reference is a low-completeness, contaminated genome — is a data limit no clustering can
overcome.

**Panel design as a first-class determinant.** A recurring, practically important finding is that
strain-level performance on real, open-world communities depends as much on the reference panel as on the
algorithm. A generic pathogen panel with few genomes per species gave a null result on real saliva, because
real strains map uniformly across arbitrary clusters when the panel does not represent the niche's genuine
diversity; a genome-rich, niche-appropriate oral panel recovered the full individual-discrimination signal.
For deployment this means panel construction (niche-appropriate species, many genomes per species, quality-
filtered) is not a preprocessing detail but part of the method, and it is the main determinant of whether
strain resolution is achievable at all on a given sample type.

**The two modes as one method, and the Fast2bRAD-M tie-in.** A conceptual contribution of the chapter is
that the fast shotgun mode and the sensitive native-2bRAD mode are not two tools but two entry points to
one tag space: the shotgun mode is validated by the mock (20/20 at 0 % host) and by the saliva concordance
(its calls are a confirmed subset of the native-2bRAD calls), and the native mode extends the same method
into the regime where shotgun fails. Because the tags are identical to Fast2bRAD-M's, completing the
species layer so that species and strain calls come from a *single* 2bRAD digest — the species gate taken
from a broad Fast2bRAD-M database and strain resolution from the per-species panels — is a natural and
high-value extension, particularly for the native-2bRAD clinical regime where a single library would then
yield both community composition and strain-level detail.