# Strain2bScan — Discussion (draft)

Strain2bScan reframes strain-level profiling around a sparse, reproducible marker set, and does so as a
single engine serving **two input modes**: native 2bRAD-M libraries and in-silico digestion of
conventional shotgun. By porting StrainScan's clustering-plus-unique-marker framework onto 2bRAD tags in
Rust, it turns a ~50–100× smaller marker database into large, structural gains **without sacrificing
accuracy** — across every species tested it held precision 1.0, matched StrainScan's detection onset at
0.5× coverage (Fig 3), and matched or exceeded StrainScan's recall both on StrainScan's own databases
(Fig 10) and on a common 15-species simulated benchmark where the two tools built databases from the same
genomes and profiled the same reads (median recall 0.80 vs 0.67 at equal precision, while building
databases 249–614× faster and profiling 4–105× faster; Fig 11). The foundation for a strain-level *2bRAD* method is that its tags, unlike the 16S gene, carry
genome-wide strain signal: across 15 species 2bRAD between-strain distances tracked the whole genome
(median Spearman 0.94) while 16S did not (0.36), several 16S intervals overlapping zero (Fig 2A,B).
16S resolves species; 2bRAD tags resolve strains.

**Native 2bRAD for low-biomass, high-host microbiomes.** The distinctive advantage of the native-2bRAD
mode is that the reduction happens at the bench, before host DNA can swamp the library. On the ATCC
MSA-1002 mock at 99 % human DNA, native 2bRAD held precision 1.0 and full 20/20 strain recall where
in-silico-digested shotgun of the same material dropped to 12/20, because 2bRAD preserves ~10× more
usable markers under host load (Fig 6). On real saliva this translated into biology: strain-level
profiles discriminated individuals better than species-level profiles (PERMANOVA R² 0.83 vs 0.82;
leave-one-timepoint-out host-ID 100 % vs 94 %), were temporally stable within subject, and — validated
against paired shotgun of the same samples — recovered 128–163 low-abundance strains per sample that
host-limited shotgun could not reach, while calling nothing the shotgun mode contradicted (Fig 7, Fig 8).
This is the setting where strain resolution matters most (oral, tumour/FFPE, skin) and where shotgun is
weakest; a wet-lab reduction that concentrates sequencing on informative tags is the natural fit, and
Strain2bScan is, to our knowledge, the first tool to make native 2bRAD strain-resolved.

**Conventional metagenomes at community scale.** In the shotgun mode, the decisive gains are per-sample
speed and memory (~8× faster, ~11× lighter than StrainScan; Fig 9A) and, above all, scaling in the
number of species: because a sample is digested once and matched against every per-species database, the
marginal cost of an additional species is a hash lookup rather than a re-count. For *S* species over *N*
samples the cost is ≈ *N × (digest + S·ε)* versus ≈ *N × S ×* (k-mer count + search) for a full k-mer
tool run per species — an *S*-fold structural advantage that, measured on a 55-species community, reached
**~132× at 100 samples and ~146× beyond** (minutes versus projected hours) and grows with community
richness (Fig 9C). The two modes are not separate tools: the shotgun mode is validated by the mock
(20/20 at 0 % host, Fig 6) and by the saliva concordance (its calls are a confirmed subset of the
native-2bRAD calls, Fig 8), so the fast shotgun mode and the sensitive 2bRAD mode are two faces of one
method.

**Cluster resolution is the honest unit of strain analysis.** Short reads cannot separate strains that
share almost all of their sequence, so both StrainScan and Strain2bScan resolve to *clusters* of
near-identical strains and we evaluate at that resolution. The clusters-to-genomes ratio at the 0.95 cut
is a per-species property — for diverse species (*C. acnes*) clusters are essentially single strains,
whereas for near-clonal *M. tuberculosis* (5 clusters from 40 genomes) the cluster is the honest level of
claim. Resolving to clusters does not cost precision (it stayed 1.0 even for low-diversity
*S. epidermidis*); what varies with diversity is recall.

**Design choices.** Occurrence-based uniqueness — a marker is cluster-unique only if absent, at any copy
number, from every other cluster — eliminates false-unique markers from single-copy filtering and keeps
precision 1.0 in similar-strain species. MinHash-sketch clustering scales database construction while
producing partitions identical to exact Jaccard and to StrainScan's own pre-built *P. copri* clustering
(112 → 51 clusters; `results/panelsize_prevotella.tsv`). Reference incompleteness is the one factor that
genuinely degrades strain identification under Jaccard — an incomplete genome's markers are a subset of a
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
distinguish strains (*M. tuberculosis*); this limit is intrinsic and shared by all short-read tools —
and where it bites hardest StrainScan failed to complete (>3.3 h, >25 GB) while Strain2bScan finished in
~1 s at precision 1.0 (Fig 10). Layering a within-cluster overlap/regression step on top of
occurrence-based uniqueness is the clearest algorithmic target for raising near-clonal recall.
(ii) **Reference panels must be niche-appropriate and genome-rich** for real communities: a generic
pathogen panel with few genomes per species gave no saliva signal because real strains map uniformly
across arbitrary clusters, whereas a genome-rich oral panel recovered the full individual-discrimination
signal — panel design is a real determinant of strain-level performance on open-world data.
(iii) **Host-limited shotgun** cannot reach the low-abundance strain tail on high-host samples; this is a
property of the input, not the tool, and is precisely the gap the native-2bRAD mode fills.
(iv) **Reference incompleteness is improved but not fully solved.** The `--containment` clustering mode
recovers accuracy when an incomplete genome has a complete relative in the panel (Fig 4), but it cannot
recover markers that are simply *absent* — a strain whose only reference is a partial assembly — nor undo
contamination that injects foreign tags, so it converges with Jaccard below ~70 % completeness; and
because it merges more aggressively it trades a little resolution on complete panels (hence opt-in, not
the default). Making strain identification *more resistant* to incomplete references is a clear direction:
a **completeness-aware detection gate** (scale the unique-marker floor by each genome's estimated
completeness, or gate on a *fraction* of a cluster's available markers rather than an absolute count —
analogous to the Layer-1 breadth term) so incomplete strains are not gated out; a
**best-quality-representative** marker set per cluster (define the cluster's markers from its most-complete
member); **upstream completeness/contamination estimation and decontamination** (CheckM2 / GUNC) feeding
the quality filter; and **pangenome-based imputation** of missing markers from complete conspecifics. The
irreducible case — a strain represented only by a low-completeness, contaminated genome — is a data limit
no clustering can overcome. (v) Aside
from the mock and saliva chapters, the accuracy benchmarks are simulated, error-free and closed-world.
Extensions include oral-cancer case/control analysis (needs the study's sample labels), FFPE and
degraded material, deeper multi-tool comparison (sylph, StrainGE), and completing the Fast2bRAD-M species
layer so species and strain calls come from one 2bRAD digest.

**Conclusion.** Reduced-representation 2bRAD markers, combined with a StrainScan-style resolution
framework and a fast Rust implementation, make accurate strain-level profiling practical at a fraction of
the compute and memory of full-k-mer methods. Uniquely, the same tool spans two regimes: native 2bRAD-M
for strain-level analysis of low-biomass, high-host microbiomes, and in-silico-digested shotgun for
strain profiling across communities of many species and many samples.
