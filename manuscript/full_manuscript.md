# Strain2bScan: strain-level metagenomic profiling on 2bRAD-reduced markers for low-biomass microbiomes and community-scale cohorts

*Assembled manuscript. Authors, affiliations, funding and data-availability accessions to be completed.*

---

## Abstract

Strain-level variation underlies clinically important microbial phenotypes, but 16S rRNA surveys usually
resolve species rather than strains. Shotgun k-mer profilers recover strain signal, yet repeated
per-species queries scale poorly across communities and host-dominated low-biomass specimens often contain
too little microbial DNA. 2bRAD sequencing samples a reproducible genomic fraction before host DNA dominates
a library, but has mainly supported species-level profiling. We present Strain2bScan, a Rust profiler that
clusters conspecific reference genomes and scores unique 25–33 bp type-IIB restriction markers. It accepts
native BcgI 2bRAD-M libraries and in-silico-digested shotgun reads in one marker framework. Across 15
species, 2bRAD-marker distances tracked whole-genome strain distances (median Spearman ρ = 0.94), whereas
16S distances were weaker (median ρ = 0.36). With complete references, precision was 1.0 across simulated
single-species benchmarks, detection reached 0.5× coverage, and accuracy depended strongly on reference
completeness. In the primary native BcgI 2bRAD-M run at 99% human DNA, Strain2bScan retained F1 = 1.0 at an
abundance threshold of 10⁻⁴ and had the lowest Bray–Curtis dissimilarity among tested tools. In 8 saliva
subjects, leave-one-timepoint-out host identification reached 100% with strain features versus 78.1% with
species features, although PERMANOVA R² was similar (0.757 versus 0.755). In three usable paired samples,
all 65 shotgun strain-cluster calls were also present in native BcgI 2bRAD-M, which yielded 128–163
additional candidate strain-cluster calls per sample. In single-panel testing, Strain2bScan was
approximately 8× faster and 11× lighter than StrainScan; a projected 55-species comparison was 121–146×
lower in cost than per-species querying. On a matched 15-species benchmark it matched StrainScan precision
(1.0), improved median recall (0.80 versus 0.67), built databases 249–614× faster and 43–138× lighter, and
completed *Klebsiella pneumoniae*, which StrainScan could not build. In the primary host-contamination
comparison, it was the only tested tool that preserved both detection and abundance at 99% human DNA. In
exploratory public WGS subsets, 22 libraries were profiled in 100.6 s using 343 MiB peak RSS; a
cohort-specific isolate panel showed the expected 6/6 panel/read compatibility result and detected a
persistent *Bifidobacterium bifidum* strain-resolved cluster across three weekly infant stool samples.
These were feasibility and concordance analyses, not independent strain-accuracy validations. Accuracy
therefore depends on reference completeness and niche-matched panel design.

**Availability.** Rust source: https://github.com/HuangShiLab/Strain2bScan. Public reads: PRJNA1131785,
PRJNA288562, PRJNA1517970, PRJNA1191223 and PRJNA1191225. Derived outputs:
https://github.com/HuangShiLab/Strain2bScan-paper.

## Introduction

Microbial strains of the same species can differ sharply in phenotype; virulence, antibiotic
resistance, metabolic capacity and host interaction; so resolving *which* strains are present, and at
what abundance, is central to microbiome and clinical metagenomics. The 16S rRNA gene, the workhorse of
community surveys, cannot see this variation: its between-strain distances are essentially uncorrelated
with genome-wide divergence, so 16S resolves species but not strains (below, and Fig 2). Shotgun
metagenomes do carry strain information, and a family of tools now recovers it; reference k-mer methods
such as StrainScan and StrainGE, marker-gene methods such as StrainPhlAn, and sketch-based abundance
estimators such as sylph; making strain-level profiling routine for individual samples and species.

Two obstacles keep shotgun-based strain profiling from two increasingly important settings.

**Scale.** Modern studies span hundreds of samples and aim to resolve strains across the
dozens-to-hundreds of species in a community. Full k-mer methods index and query the entire k-mer
content of every reference genome; per-sample cost is dominated by counting the sample's k-mers, and is
paid again for every species database queried. For *S* species across *N* samples this scales as
*N × S × (k-mer count + search)*; hours-to-days of compute and hundreds of megabytes to gigabytes of
memory; a practical ceiling for population-scale, multi-species strain surveillance.

**Low biomass and host contamination.** Many of the most interesting clinical niches; saliva and other
oral sites, tumour and FFPE tissue, skin; yield little microbial DNA against a large human background
(often ≥90–99 % host). Shotgun sequencing spends most of its reads on the host, so the informative
bacterial fraction, and with it the rarer strains, is lost; strain resolution collapses exactly where it
would be most valuable.

Reduced-representation sequencing addresses both. Type-IIB restriction ("2bRAD") digestion releases a
sparse, reproducible set of fixed-length tags (roughly 1–2 % of the genome) and, applied as a wet-lab
protocol (2bRAD-M / Fast2bRAD-M), it enriches these informative markers *before* host DNA can swamp the
library, giving accurate profiles from picogram inputs, heavily host-contaminated samples and degraded
material. To date, however, 2bRAD has been used only for **species**-level profiling. Yet strain-resolution
methods do not need the whole genome; they score samples on *unique* markers that distinguish a strain or
a cluster of near-identical strains; and a 2bRAD tag set is exactly this kind of low-redundancy,
taxonomically informative marker set, 50–100× smaller than the full k-mer set.

We present **Strain2bScan**, which ports the StrainScan resolution framework; within-species clustering
into a search structure, followed by unique-marker scoring and abundance estimation; onto 2bRAD tags, in
Rust. Its tag lengths and recognition patterns match Fast2bRAD-M / `2bRADExtraction.pl`
exactly, so the tags are interoperable with the Fast2bRAD-M species layer, and Strain2bScan accepts **two input modes** that map onto the two obstacles above:

1. **Native BcgI 2bRAD experimental libraries**, enabling strain-level analysis of
   low-biomass and high-host microbiomes. Because the reduction happens at the bench, Strain2bScan holds
   precision 1.0 and full strain recall at 99 % host DNA where in-silico-digested shotgun loses most of
   its strains, and on real saliva it resolves individual-specific, temporally stable strain signatures;
   recovering low-abundance strains that host-limited shotgun cannot reach.

2. **In-silico digestion of conventional shotgun metagenomes**, enabling community-scale strain
   profiling. The sample is digested **once** and matched against every per-species database, so
   per-sample cost is independent of the number of species and linear in the number of samples: ~8×
   faster and ~11× lighter per sample than StrainScan, and projected to be 121–146× lower in cost on a 55-species community,
   while matching StrainScan's precision, its 0.5× detection onset and its recall on its own databases;
   and completing near-clonal *Mycobacterium tuberculosis* in ~1 s where StrainScan does not finish.

The two modes are directionally concordant: in three usable paired samples, 65/65 shotgun strain-cluster calls were contained in native calls. Thus a single
tool, on one 2bRAD tag space, spans both the low-biomass clinical regime and the cohort-scale regime.

## Results

#### Overview of Strain2bScan and its two data modes (Fig 1)

Strain2bScan reduces reference genomes to canonical, single-copy 2bRAD tags, clusters genomes within each
species at 0.95 Jaccard similarity, and builds a compact cluster-by-marker database. Profiling digests a
sample once, gates species on panel-specific markers, and scores strain-resolved clusters on unique
markers. The same tag space accepts two entry points: in-silico digestion of conventional shotgun reads and
native BcgI 2bRAD-M libraries whose reads are already tags. The tag definitions match Fast2bRAD-M, so its
species layer and the Strain2bScan strain layer share one marker space.

#### 2bRAD tags carry strain-level signal that 16S cannot (Fig 2)

We compared pairwise strain distances in whole-genome, 2bRAD-tag and 16S spaces across 15 species with
complete or near-complete genomes. Per species, 2bRAD distances tracked whole-genome distances in every
case (median Spearman ρ = 0.94; range 0.59–0.99), whereas 16S distances were much weaker (median ρ = 0.36),
with several confidence intervals overlapping zero (Fig 2). The rank-rank matrices show the mechanism:
2bRAD preserved the ordering of genome-wide strain pairs, whereas 16S collapsed many unrelated pairs to a
few conserved-gene distances. Thus, 16S resolved species but the 2bRAD marker set retained genome-wide
strain signal.

#### Accurate and robust strain profiling (Fig 3–4)

On simulated mixtures from curated reference panels, Strain2bScan achieved precision 1.0 for *Cutibacterium
acnes*, *Staphylococcus aureus* and *Staphylococcus epidermidis*, with high recall and Bray–Curtis
dissimilarity of 0.24, 0.33 and 0.02, respectively (Fig 3A). A single *C. acnes* strain was detected at
0.5× per-strain coverage, matching the StrainScan detection onset (Fig 3B). Precision therefore remained
high even when many genomes collapsed to few 0.95-similarity clusters.

Reference completeness was the dominant accuracy constraint. Degrading truth-strain references from 100% to
70% completeness while holding sample reads fixed reduced median Jaccard-mode precision from 1.0 to 0.71 and
recall from 0.96 to 0.74 across 14 resolvable species, because incomplete genomes fragmented from complete
relatives (Fig 4). The optional `--containment` mode, which links genomes by maximum marker-set containment,
restored precision to 0.92–0.98 and recall to 0.92–0.95 at 90–95% completeness. It also kept the near-clonal
*M. tuberculosis* cluster intact to 90% completeness, whereas Jaccard recall collapsed to approximately
0.05. Because containment can coarsen complete panels, it remains opt-in for uneven-completeness references.

### Part I; Native 2bRAD-M for low-biomass, high-host microbiomes

#### A tunable enzyme set enables native BcgI 2bRAD operation (Fig 5)

Across 14 resolvable simulation-pool species, increasing the enzyme set from 1 to 14 left precision at 1.0
while increasing median recall monotonically from 0.50 to 1.00 and the median strain-specific marker yield
from 976 to 12,647 (Fig 5). Cluster count was invariant, so additional enzymes increased marker density and
recall rather than resolution. BcgI alone under-resolved near-clonal species, whereas approximately four
enzymes recovered most recall at much lower marker cost. Single-enzyme BcgI operation nevertheless enables
direct profiling of native BcgI 2bRAD-M libraries.

#### Strain-level identification and quantification across four DNA mocks (Fig 6)

We next profiled four ATCC whole-cell mocks against a unified 28-species, 164-genome tree in which each mock
species was represented by its ATCC genome and up to five conspecific decoys (ANI > 95%). Native BcgI reads
were clustered with `--containment` and profiled with `--min-abundance 0 --min-coverage 0.2`; the primary
detection threshold was abundance ≥ 10⁻⁴. At 90, 95 and 99% human DNA, F1 was 0.976–1.0 and Bray–Curtis
similarity was 0.745–0.756 (Fig 6). In the DNA-input ladder, F1 was 0.952 at 0.1 ng (precision 0.909, recall
1.0), 1.0 at 1 ng and 0.625 at 0.01 ng; AUPR fell to 0.50 at 0.01 ng, and 0.001 ng was below the practical
detection floor. On staggered MSA-1003, F1 was 0.703–0.824 and AUPR was 0.65–0.70. MSA-1005 retained F1 =
1.0; MSA-1007 had F1 = 0.667–1.0 because a small number of abundance-negligible calls crossed the presence
threshold, although AUPR was 1.0. Together, the mocks reveal a ~1× marker-depth noise floor: true strains
were deeply covered, whereas spurious near-sibling calls occurred near 1× and approximately 0.02% relative
abundance. This floor was larger on the unified 28-species database than on a per-mock 20-species tree
(Fig S3).

#### Real saliva: individual-specific, temporally stable strain signatures (Fig 7)

We profiled native BcgI 2bRAD saliva from 8 subjects sampled at four times of day (32 libraries) against a
19-species oral reference panel. Subject PERMANOVA R² was similar for strain and species profiles (0.757
versus 0.755; both p = 2 × 10⁻⁴), but leave-one-timepoint-out nearest-neighbour classification identified
the host with 100% accuracy from strain features versus 78.1% from species features (Fig 7). The strongest
single-species association was *Neisseria subflava* (subject R² = 0.808); all 13 testable species were
significant at nominal p values. Within-subject Bray–Curtis distance was lower than between-subject distance
(0.327 versus 0.696; p = 8.05 × 10⁻¹⁹), indicating a stable, person-specific strain signature.

#### Native 2bRAD detects candidate low-abundance strains missed by truncated shotgun (Fig 8)

Paired shotgun allowed a directional comparison of input modalities, not an independent accuracy validation.
Only three subjects had usable paired bacterial recovery, and the shotgun comparison used R1 prefix
subsamples. In those samples, all 65 shotgun strain-cluster calls were also present in native BcgI 2bRAD-M
(65/65). Native BcgI 2bRAD-M yielded 128–163 additional candidate strain-cluster calls per sample. Calls
unique to 2bRAD had lower median community relative abundance than shared calls (0.0029 versus 0.0097;
Mann–Whitney p = 1.2 × 10⁻²³, with calls treated as observations rather than independent subjects; Fig 8).
The pattern is consistent with greater sensitivity to low-abundance strains on native high-host material, but
the additional calls require independent validation. Strain2bScan also profiled four exploratory oral
libraries, producing 15–17 species and 115–158 strain-cluster calls per sample in approximately 2 s
(Table S3).

### Part II; Conventional metagenomes at community scale

#### Fast, light, and scalable to whole communities (Fig 9)

For shotgun input, Strain2bScan profiled the *C. acnes* panel in 0.86 s and 78 MB per sample, versus 7.06 s
and 828 MB for StrainScan, an approximately 8× runtime and 11× memory advantage (Fig 9A). Database build and
profiling parallelised to 4.6× and 5.8× speed at 16 threads (Fig 9B). Community profiling was approximately
flat in species count: on a 55-species panel, Strain2bScan profiled samples in 2.6–3.1 s. Because StrainScan
has no multi-species mode, its projected cost was the measured single-species runtime multiplied by species
and sample count; the projected advantage was 121–146× (Fig 9C).

#### Matches or exceeds StrainScan on its own databases (Fig 10)

On StrainScan-curated reference sets, both tools reached precision 1.0 for *A. muciniphila* and *P. copri*
(Fig 10). Strain2bScan matched or exceeded recall (0.93 versus 0.24 and 0.94 versus 0.90), was 17–23×
faster and 15–24× lighter, and completed near-clonal *M. tuberculosis* in 0.89 s, whereas StrainScan did not
complete. Its low *M. tuberculosis* recall reflects the resolution limit of a panel that collapses to five
0.95-similarity clusters.

The shotgun mode was then stress-tested in ATCC mocks and compared with saliva. It preserved detection and
abundance under host contamination in the primary MSA-1002 comparison (Fig 12), and its saliva calls were
contained in the native-2bRAD call set (Fig 8). These are concordance and stress-test results, not
independent validations of strain identity.

#### Systematic head-to-head on a 15-species simulated benchmark (Fig 11, Table 1)

On a common 15-species simulation pool, both tools built databases from the same genomes and profiled the
same simulated reads. Each tool was scored in its own 0.95-similarity cluster space. Across 204 paired
single-species samples from 14 resolvable species, median precision was 1.0 for both, but Strain2bScan
reached full recall by 3× coverage versus 10× for StrainScan; median recall and F1 were 0.80 and 0.89 versus
0.67 and 0.75 (Fig 11, Table 1). Strain2bScan held precision 1.0 in every species, whereas StrainScan fell
to 0.80–0.83 in four species. Multi-species accuracy was comparable, with Strain2bScan stronger at low depth
and StrainScan slightly stronger at high depth (Table 3).

The cost differences were largest during database construction: Strain2bScan used 0.7–5.1 s and 0.1–0.4 GB
per species, versus 5–43 min and 8–28 GB for StrainScan, a 249–614× runtime and 43–138× memory advantage
(Table 1). In the same emulated container, per-sample profiling was 4–33× faster and 5–39× lighter
(Table 2). Because StrainScan lacks a multi-species mode, its community cost was the sum of 14 runs
(100–398 s), versus 1–9 s for one digest-once pass (46–105× faster). StrainScan also failed to build
*Klebsiella pneumoniae* within the resource cap, whereas Strain2bScan built it in 5.1 s. StrainScan build
times were obtained under linux/amd64 emulation and are upper bounds.

#### Strain-level profiling on shotgun, and the advantage under host contamination (Fig 12)

We compared Strain2bScan, StrainScan v1.0.14 and inStrain 1.10.0 on shotgun reads from the same four ATCC
mocks. Strain2bScan used the all-enzyme 164-genome tree; StrainScan used per-species databases; inStrain
used a dereplicated 98%-ANI reference, as its documentation requires. Each tool was scored in its own
0.95-similarity cluster space.

At the primary 10⁻⁴ threshold, each tool scored well on the MSA-1002 0%-host sample. The tools then
separated under host contamination. Strain2bScan retained F1 = 1.0 at 90, 95 and 99% human DNA; Bray–Curtis
similarity was 0.777, 0.772 and 0.698, and L2 similarity was at least 0.839. At 99% host, StrainScan had F1
= 0.095 because thresholded recall collapsed to 0.05 and its abundance diverged (Bray–Curtis similarity
0.026); inStrain had F1 = 0.333 because recall fell to 0.20. In this primary comparison, Strain2bScan was
the only tested tool that retained both strain detection and abundance at 99% host. Clean-mock behaviour was
more mixed. On MSA-1003, Strain2bScan had F1 at 1e-4 = 0.909, 1.0, 1.0 and Bray–Curtis similarity
0.773–0.780. On MSA-1005/1007 it retained full recall and AUPR = 1.0, but F1 was 0.60–0.75 because many
abundance-negligible clusters crossed the threshold; StrainScan had higher thresholded precision on those
two mocks. The host-contamination result therefore should not be read as uniform clean-sample superiority.

#### Exploratory public WGS cohorts: longitudinal observations, low-biomass specificity, and isolate-derived panels (Tables 4–6)

To test feasibility on public human metagenomes, we profiled representative WGS subsets from PRJNA288562,
PRJNA1517970, PRJNA1191223 and PRJNA1191225. These were application-oriented subsets, not complete
epidemiological cohorts. Twenty-two libraries profiled against a legacy 20-species MSA panel completed in
100.6 s with 343 MiB peak RSS (4.57 s per sample).

In PRJNA288562, one pregnancy subject had saliva, vaginal and distal-gut WGS at two gestational ages. Five
of six libraries produced 34 strain-cluster calls across six panel species. Saliva was cluster-rich (13 and
11 calls at the two timepoints), whereas vaginal signal was sparse and gut calls changed over time. Four
*Neisseria* and four *Schaalia* clusters were shared between early and late saliva; a gut *Bifidobacterium
adolescentis* cluster persisted, while other *Bifidobacterium* and *E. coli* clusters were detected only
early (Table 4).

PRJNA1517970 tested vaginal and meconium libraries plus one extraction blank. Neither the generic panel nor
a newly built 13-species body-site panel produced a strain call. The blank yielded 340 distinct markers,
whereas selected vaginal and meconium libraries yielded 7,394–86,980; a one-marker diagnostic gate found only
three *C. acnes* markers in the largest meconium library. These results indicate sparse marker input and
incomplete niche coverage rather than runtime failure, and the negative blank supports specificity at low
marker input.

For direct gut tracking, we assembled six public isolates from PRJNA1191225 (three *E. coli*, one
*B. longum*, one *B. breve* and one *B. bifidum*) and built a cohort-specific panel. All six isolate read
sets produced their expected panel unit (6/6 panel/read compatibility result) with breadth 0.9995–1.0, depth
9.8–25.6× and within-species abundance 0.80–0.96; profiling took 9.53 s and 55.7 MiB (Table 5). Because the
same reads generated the assemblies and the panel, this was a compatibility check rather than an independent
classification benchmark.

Reprofiling preterm-infant P08 across weeks 1–3 with the cohort panel detected a persistent
*Bifidobacterium bifidum* LHCA82 strain-resolved cluster (breadth 0.764/0.773/0.792; depth 237×/214×/238×;
within-panel abundance 64.6%/58.3%/61.0%). *E. coli* was detected only in week 2 at 0.029× depth and sample
fraction 1.40 × 10⁻⁴; because coverage did not distinguish the three study isolates, the tool conservatively
reported the merged unit `C0|C1|C2`. The generic MSA panel detected the week-2 *E. coli* signal but lacked
the cohort-specific *B. bifidum* reference. Because the panel and metagenomes came from the same cohort and
no independent strain truth was available, these are within-panel longitudinal observations. Leave-one-
isolate-out *E. coli* tests reinforced this boundary: reads from an absent isolate covered the remaining
relatives broadly (0.799–0.911) and were reported as merged `C0|C1` rather than creating a false third
cluster. Thus, closed panels can misassign conspecific signal, and absence of a true strain cannot be
inferred from a nearest-relative call (Table 6).

## Discussion

Strain2bScan couples a sparse 2bRAD marker space to a StrainScan-style clustering and unique-marker
framework in a fast Rust implementation. Across the tested settings, 2bRAD tags preserved genome-wide strain
ordering that 16S did not (median Spearman 0.94 versus 0.36), precision was 1.0 in simulated single-species
benchmarks, and detection matched StrainScan at 0.5× coverage. On a common 15-species benchmark, Strain2bScan
improved median recall (0.80 versus 0.67) at equal precision while reducing database build by 249–614× and
profiling by 4–105× (Figs 2, 3, 10 and 11).

**Native 2bRAD for low-biomass, high-host microbiomes.** Native reduction occurs before host DNA dominates a
library. In the primary MSA-1002 comparison at 99% human DNA and a 10⁻⁴ abundance threshold, Strain2bScan
retained F1 = 1.0 with Bray–Curtis dissimilarity 0.302; StrainScan had thresholded recall 0.05 and abundance
dissimilarity 0.974, whereas inStrain had F1 = 0.333 (Fig 12). In saliva, leave-one-timepoint-out host
identification reached 100% for strain profiles versus 78.1% for species profiles, although PERMANOVA R² was
similar (0.757 versus 0.755). In three usable paired samples, 65/65 shotgun strain-cluster calls were present
in native 2bRAD-M, which produced 128–163 additional candidate strain-cluster calls per sample (Fig 8). These
are concordance and sensitivity observations, not independent validation. Tumour, FFPE and skin applications
remain untested but are natural targets because host content and degradation can limit shotgun sequencing.

**Conventional metagenomes at community scale.** In shotgun mode, Strain2bScan was approximately 8× faster
and 11× lighter in single-panel testing. Its main structural advantage is that one digest serves all species,
so the marginal cost of another species is a hash lookup. In a 55-species analysis, measured Strain2bScan
runtime was 121–146× lower than the corresponding projected per-species StrainScan cost (Fig 9C). The primary
high-host mock supported detection and abundance retention (Fig 12), and saliva calls were contained in the
native-2bRAD call set (Fig 8), but the clean-mock precision trade-off means these results support a shared
framework rather than uniform superiority.

**Cluster resolution is the honest unit of strain analysis.** Short reads cannot separate strains that
share almost all of their sequence, so both StrainScan and Strain2bScan resolve to *clusters* of
near-identical strains and we evaluate at that resolution. The clusters-to-genomes ratio at the 0.95 cut
is a per-species property; for diverse species (*C. acnes*) clusters are essentially single strains,
whereas for near-clonal *M. tuberculosis* (5 clusters from 40 genomes) the cluster is the honest level of
claim. Resolving to clusters does not cost precision (it stayed 1.0 even for low-diversity
*S. epidermidis*); what varies with diversity is recall.

**Design choices.** Occurrence-based uniqueness; a marker is cluster-unique only if absent, at any copy
number, from every other cluster; eliminates false-unique markers from single-copy filtering and keeps
precision 1.0 in similar-strain species. MinHash-sketch clustering scales database construction while
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
member); **upstream completeness/contamination estimation and decontamination** (CheckM2 / GUNC) feeding
the quality filter; and **pangenome-based imputation** of missing markers from complete conspecifics. The
irreducible case (a strain represented only by a low-completeness, contaminated genome) is a data limit
no clustering can overcome. (v) The public-cohort experiment in Table 4 adds real WGS profiling of longitudinal, multi-site and
low-biomass samples, and Table 5 adds an isolate-derived panel compatibility check. These are application-oriented
subsets rather than closed-world truth benchmarks: read sets used to build the isolate panel were also
used for self-recovery, and the real metagenomes lack exhaustive strain truth. Deeper external
validation, oral-cancer case/control analysis, FFPE and degraded material, broader multi-tool
comparison (sylph, StrainGE), and completing the Fast2bRAD-M species layer remain future work.

(vi) The PRJNA1517970 body-site panel still produced no calls, showing that species choice alone does not
solve low-marker input; panel completeness and sequencing depth remain limiting. (vii) Leave-one-isolate-out
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

**Conclusion.** Reduced-representation 2bRAD markers, combined with a StrainScan-style resolution
framework and a fast Rust implementation, make accurate strain-level profiling practical at a fraction of
the compute and memory of full-k-mer methods. Uniquely, the same tool spans two regimes: native 2bRAD-M
for strain-level analysis of low-biomass, high-host microbiomes, and in-silico-digested shotgun for
strain profiling across communities of many species and many samples.

## Methods

### Overview

Strain2bScan reimplements the two-layer StrainScan strategy; cluster near-identical strains,
then score samples on markers unique to a strain or cluster; but replaces the full k-mer set
with **2bRAD tags**, and implemented in Rust for speed and parallelism. The primary benchmark used the dependency-free release `benchmark-v0.1.0-f26f234`; later development versions are not part of that frozen configuration. The
pipeline is: (i) digest reference genomes and sample reads into 2bRAD-tag markers; (ii) build,
per species, a within-species cluster database annotated with unique markers; (iii) profile a
sample by detecting present clusters from their unique markers and estimating their abundance;
(iv) at the community level, digest each sample once and match it against every per-species
database, gated by a species-level Layer-1 check.

### 2bRAD tag extraction

Type-IIB restriction enzymes cut on both sides of their recognition site, releasing a
fixed-length fragment (the 2bRAD tag, 25–33 bp depending on enzyme). Each of the 16 enzymes in
the Fast2bRAD-M table is modelled as a set of anchored sequence patterns; literal motifs at
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

### An alternative marker source: FracMinHash sketching

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

The sketch was verified to behave as FracMinHash requires. Marker count tracks 1 / *S* to within
sampling noise across three decades of scale (ratios 0.995–1.008 relative to the unsketched set,
on *E. coli* K-12), and the FNV-1a hash produces no collisions on the 989,962 distinct canonical
31-mers of a 1 Mb window, against 0.027 expected for a uniform 64-bit hash. The apparent
shortfall against genome length (4,523,995 markers from 4,641,622 windows) is entirely the
single-copy filter removing the 30,274 multi-copy 31-mers, matching an independent count exactly.

### Reference database construction

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
the other; the containment estimator used by Mash-screen and sourmash for uneven-completeness genomes.
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
complete twin's, so their Jaccard falls below 1 and they fail to cluster. Because CheckM is
not run in-line, two dependency-free proxies computed from data already at hand are used;
contig count (`--max-contigs`), and single-copy tag count relative to the conspecific median
(`--min-tag-fraction`, a completeness proxy). Genomes far below the median are always flagged;
they are removed only when a threshold is set.

### Layer-2: detection and abundance

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

### Multi-species profiling and species selection

For community samples, the reads are digested **once** into a shared set of tag counts, matched
against every per-species cluster database in parallel; the per-species marginal cost is a
hash-set lookup rather than a re-count, so the total cost is independent of the number of species.

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

### Cross-species restriction on quantification

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

### Ported StrainScan layers, and why neither is the default

Both stages of StrainScan's resolution framework are implemented and selectable, so the
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

### Implementation

The benchmark release `benchmark-v0.1.0-f26f234` is written in Rust with no third-party dependencies, making that frozen binary self-contained. Current development versions may use additional dependencies and are not part of the frozen benchmark. Data-parallelism (genome digestion, sketch
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

**Measured.** On 4 M reads (289 MB) against a two-species panel with all 16 enzymes, wall time is
0.51 s and peak resident memory 11 MB (1.71 s / 282 MB before these optimisations); gzipped input
costs no additional wall-clock. Tree construction scales as ~O(*n*<sup>1.6</sup>) after replacing a
per-merge recomputation of max-linkage with an incrementally updated cluster-level similarity
matrix, and a per-node set-subtraction with a single carrier-set index; on 543 genomes the
pairwise similarity scan that dominates a large build is parallel, taking the build from 71.4 s to
58.4 s. Every optimisation was verified to leave output unchanged; for the tree, node for node
against the serial build at *n* = 543.

Two intuitive optimisations were **rejected on measurement** and are recorded here because both
are commonly assumed to help: replacing the hash-set Jaccard with a sorted-vector two-pointer
intersection is 0.82×; slower, because a hash lookup on `u64` is cheap while the two-pointer must
traverse both arrays; although it would halve the memory; and fusing the multi-enzyme scan into a
single pass is 1.43× on a 2.5 Mb contig but 0.93× on 150 bp reads, so it would have to be
dispatched on sequence length rather than applied globally.

### Benchmarking

**Datasets.** (i) A real *C. acnes* benchmark: a 64-genome reference panel (14 ground-truth
strains + 50 background, NCBI accessions pinned) and five paired-end mock samples from
MockMetagenomes4Benchmark (~100k read pairs each, ~12× total). (ii) A simulated multi-species
benchmark: 55 real species × ~4 strains (218 NCBI genomes) and 30 samples, each mixing strains
from twelve species at log-normal depth ≥1× (plus a low-depth variant, median 0.62×, used for
gate calibration). (iii) Cross-species mocks for *Staphylococcus
aureus* and *S. epidermidis* (60-genome panels each; 2–5 strains/sample, log-normal ≥1×,
matching the *C. acnes* design). (iv) A reference-degradation gradient in which the truth
strains' database genomes are degraded to completeness 100→50 % (with co-varying contamination
0→10 % and fragmentation), samples held fixed. Simulated reads were generated with ART rather than treated as error-free. The systematic 15-species
benchmark used `art_illumina -p -l 250 -m 600 -s 150`, producing 250 bp paired-end reads with the ART
Illumina quality-error model. Earlier diagnostic datasets may have used different read configurations;
the systematic comparison and its figures are defined by this ART configuration.

**Real-data and motivation datasets.** (v) *2bRAD-vs-16S motivation* (Fig 2): 15
pathogenic/commensal species, ~50 genomes each from NCBI accession lists (ENA FASTA), **restricted to
complete/near-complete assemblies** (CheckM completeness ≥ 97 %, contamination ≤ 5 %, assembly level
Complete Genome/Chromosome; `data/genome_qc_16s_panel.tsv`). Between-strain distance was computed in
three spaces; whole-genome (bottom-3000 canonical 21-mer MinHash), 2bRAD (Strain2bScan `build` BcgI
tags) and 16S (longest gene per genome via barrnap 0.9 + HMMER, 21-mer Jaccard); all with the Mash
transform D(J) = −ln(2J/(1+J)); per species the 2bRAD and 16S pairwise vectors were correlated (Spearman)
against the whole-genome vector, with 95 % CIs from 500 genome subsamples. (vi) *ATCC DNA mocks,
strain-level (Fig 6, Fig 12, Fig S3, Fig S4)*: four whole-cell mocks; MSA-1002 (20 strains,
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
design), 5 replicates per cell; 2 025 samples. *Multi-species* samples mixed ~18 co-present species
(one to a few strains each) across three community depth gradients; 60 samples. Reads were simulated
with ART (`art_illumina -p -l 250 -m 600 -s 150`, error-modelled 250-bp paired-end reads) from the truth genomes; truth tables record each strain's
species, genome accession and 0.95-cluster assignment.

Both tools built their databases from the **same genome pool** and profiled the **same reads**.
Strain2bScan databases were built with `cluster --enzyme all --similarity 0.95` and profiled with
`profile` / `multi-profile --enzyme all` (reads decompressed, R1+R2 concatenated). StrainScan (v1.0.14,
bioconda) is Linux-x86-only; it ships `dashing_s128` and `jellyfish-linux` ELF binaries, a Python-3.7
`.so`, and an R reclustering step; so it was run inside a Docker `linux/amd64` container (QEMU emulation
on Apple Silicon; `strainscan_build`, then `strainscan -i R1 -j R2 -d DB`). Because the two tools cluster
genomes independently, **each tool was scored in its own cluster space**: predicted clusters were compared
against the truth strains mapped into that tool's clusters; for Strain2bScan via the truth `cluster`
column, for StrainScan via its `Cluster_Result/hclsMap_95_recls.txt` (report `Cluster_ID` = `C`+cluster
id); and precision/recall/F1 computed over the cluster sets per sample. Strain2bScan profiled all 2 025 +
60 samples; StrainScan profiled a matched subset (different-cluster mixtures, k = 2/3/5, one replicate, all
depths; near-clonal *M. tuberculosis* via its same-cluster samples); 204 depth-matched paired
single-species samples across 14 species, plus 4 multi-species samples per depth. StrainScan has no
multi-species mode, so each community sample was profiled once per species database and the per-sample
cost taken as the sum of wall-clock over species (peak RSS as the maximum).

*Software and timing provenance.* The primary ATCC mock benchmark was frozen at Strain2bScan commit
`f26f234b817ba7772a3f1df59ce720751e9b45b9` (release-build SHA-256
`a4cf7a4043a06c55a99d6abf1e9fd312f6243aa5a5ca021bd41915636fb56b18`, Rust 1.97.0); later software commits
were not used for those mock outputs. The primary configuration, database hashes, metric definitions, and
per-run outputs are in `results/benchmark_configuration.json` and `results/mock_benchmark_f26f234/`.
Strain2bScan build and profile times for the 15-species benchmark were native arm64. To reduce the emulation
confound in profiling speed, a `linux/amd64` Strain2bScan binary was run in the same container on the same
subset, giving the same-environment ratio in Fig 11E/Table 2. StrainScan build times were obtained under
`linux/amd64` QEMU emulation and are therefore upper bounds. DB build for *K. pneumoniae* did not complete
under StrainScan and that species is omitted from the paired accuracy set. Complete raw StrainScan per-sample
outputs and the original `scratchpad/eval` drivers were not retained. The available Strain2bScan per-sample
tables, aggregate tables, and frozen ATCC mock outputs are retained, but the 15-species StrainScan comparison
is therefore reproducible only at the aggregate-table level and is declared a provenance limitation.

**Comparison to StrainScan (curated-DB and per-sample benchmarks).** In addition to the common benchmark
above, StrainScan v1.0.14 was run on its **own** reference databases (Fig 10) and on the same *C. acnes*
per-sample profiling comparison (Fig 9A), using its low-depth modes for the depth series.

**Primary ATCC mock configuration.** The primary Strain2bScan variant was the default flat path on the
164-genome containment tree, with no trace-gap filter and no Layer-1/Layer-2 override. Native BcgI libraries
used the BcgI database and shotgun libraries used the all-enzyme database. Optional trace-gap and port-layer
runs were retained as sensitivity artifacts but are not primary evidence and are not shown in Figures 6 or
12. The primary detection threshold was abundance ≥ 10⁻⁴.

**Metrics.** For the primary ATCC mock comparison, detection precision, recall and F1 use an abundance
threshold of 10⁻⁴; AUPR is the threshold-free abundance-ranking summary. Abundance error is Bray–Curtis
dissimilarity or L2 distance against sequence-abundance truth, and corresponding similarities are reported
as one minus dissimilarity where noted. All mock metrics are evaluated at 0.95-similarity cluster
resolution. Simulation detection metrics are computed over each tool's predicted and truth cluster sets.
Wall-clock time and peak resident set size were measured with `/usr/bin/time` on a 16-core Apple-silicon
machine. Primary ATCC predictions, commands, database hashes, and checksums are archived in
`results/mock_benchmark_f26f234/` and `results/benchmark_configuration.json`. Available figure recipes are
in the repository; the complete 15-species run-driver provenance limitation is stated above.
#### Public real-metagenome cohorts and isolate-derived panels

To supplement the controlled analyses, we profiled representative paired-end WGS subsets from four
public BioProjects: PRJNA288562 (pregnancy saliva/vagina/distal-gut time series), PRJNA1517970 (preterm
vaginal and meconium metagenomes), PRJNA1191223 (preterm-infant stool time series) and PRJNA1191225
(preterm-infant isolate WGS). These were exploratory application subsets selected to test multi-site
stability, low-biomass specificity and isolate-derived strain recovery; they were not complete
epidemiological cohorts. For the generic-panel screen we used the existing 20-species MSA database in
`all`-enzyme mode with `--min-species-markers 20 --min-species-detect 2 --min-support 2
--min-coverage 0.01 --min-abundance 0`. PRJNA288562 analyses used subject T23 at gestational days 84 and
273 for all three body sites. PRJNA1517970 analyses used three vaginal, three meconium and one negative
extraction-control library.

For cohort-specific compatibility testing, six PRJNA1191225 isolate read sets were assembled with SPAdes 4.3.0
(`--isolate`, 8 threads). Three *E. coli* assemblies were clustered at 0.95 similarity into three
resolvable units; one assembly each from *B. longum*, *B. breve* and *B. bifidum* was built as a
single-genome species database. Panels used `--enzyme recommended`. All six isolate read sets were then
profiled with the cohort panel using `--min-species-markers 50 --min-species-detect 3 --min-support 2
--min-coverage 0.01 --min-abundance 0`. Because these reads also generated the assemblies, the 6/6
self-recovery is a panel/read compatibility check, not an independent classification benchmark. P08
W1–W3 metagenomes were reprofiled with the same cohort panel using `--min-species-markers 20
--min-species-detect 2` and otherwise the gates above. Timings were measured on a local Apple Silicon
workstation; full command lines and per-sample outputs are provided in
`results/realworld_cohort_benchmark/`.
To avoid reference mismatch in PRJNA1517970, we built a body-site-oriented panel from 13 species commonly
representing vaginal, neonatal-gut and meconium communities: *Bifidobacterium adolescentis*,
*Cutibacterium acnes*, *Enterococcus faecalis*, *Escherichia coli*, *Gardnerella vaginalis*,
*Lactobacillus gasseri*, *Lactobacillus jensenii*, *Prevotella bivia*, *Staphylococcus aureus*,
*Staphylococcus epidermidis*, *Streptococcus agalactiae*, *Streptococcus mitis* and *Streptococcus
mutans*. We used up to six RefSeq or ATCC genomes per species (146 genomes; two *S. mitis* genomes were
available). Databases were built with `--enzyme recommended` and profiled with the same thresholds as the
generic-panel screen. To test behaviour when a true strain is absent from the panel, we also built three
leave-one-isolate-out *E. coli* databases, each containing only the other two *E. coli* assemblies plus
the three Bifidobacterium controls, and profiled the held-out isolate reads.

## Figure legends

### Main figures

**Figure 1. Strain2bScan overview and the two input modes.**
Pipeline schematic. (A) Reference construction: type-IIB (2bRAD) digestion of reference genomes into
single-copy 25–33 bp tags → within-species clustering at 0.95 Jaccard (MinHash-accelerated) → a
cluster × marker database of species-core, cluster-specific and strain-specific markers. (B) Profiling:
a sample is digested once into canonical markers, gated on species-specific markers (Layer-1), and
strains are detected and quantified within each present species from unique markers (Layer-2). The two
input modes; in-silico digestion of conventional shotgun, or native BcgI 2bRAD-M libraries whose reads
are the tags; enter the same tag space, shared with the Fast2bRAD-M species layer.

**Figure 2. 2bRAD tags track genome-wide strain distance; 16S does not.** **(A)** Per-species Spearman correlation between
whole-genome (bottom-3000 21-mer MinHash) between-strain distance and the corresponding 2bRAD-tag (blue)
or 16S rRNA (red) distance, across 15 species (complete/near-complete genomes only, CheckM ≥97 %/≤5 %;
14–50 genomes per species, shown at right); error bars are 95 % CIs over genome subsamples; genus and
species on separate label lines; 2bRAD median 0.94 vs 16S median 0.36. **(B)** 3×5 matrix of per-species
rank–rank scatters; rank of each strain pair's whole-genome distance (x) vs its 2bRAD (blue) or 16S (red)
distance (y): 2bRAD hugs the diagonal in every species, 16S forms flat rank-bands (most extreme in
*P. dorei*, *M. tuberculosis*).

**Figure 3. Accurate strain profiling and depth sensitivity.** (A) ;  precision/recall/abundance error (Bray–Curtis) on real reference panels with simulated mixtures for
*C. acnes*, *S. aureus*, *S. epidermidis* (precision 1.0 throughout). (B) ;  detection of a single *C. acnes* strain across a coverage ladder; onset at 0.5×, matching StrainScan.

**Figure 4. Reference-genome completeness controls strain-identification accuracy (all 15 species).** Sample reads held fixed while the truth strains' reference genomes are
degraded (completeness 100→50 %, with co-varying contamination + fragmentation), DB rebuilt and re-profiled.
Problem → mechanism → solution. **(A)** schematic (): an incomplete
genome's tags are a subset of a complete relative's, so Jaccard (6/10 = 0.60) drops below the 0.95 cut and
splits them; the shared tags land in two clusters and are demoted to non-discriminating *SharedPartial*;
max-containment (6/6 = 1.0) merges them into one cluster whose marker set is the **union** of members, so
the discriminating *ClusterSpecific* tags (and the complete member's markers) are preserved. **(B)** median
precision, **(C)** median recall vs reference completeness (14 resolvable
species; faint lines = per-species containment), comparing **default Jaccard** clustering (grey dashed)
with the **`--containment`** mode (solid, shaded gap). Degrading references collapses Jaccard accuracy
(precision 1.0→0.84→0.71, recall 0.96→0.80→0.74 at 100→90→70 %) because incomplete genomes split from
complete relatives; `--containment` (max-containment) keeps them together, restoring precision to
0.98/0.92 and recall to 0.95/0.92 at 95/90 %, converging with Jaccard only at ≤70 % (genuinely
low-quality). Inset (C): near-clonal ***M. tuberculosis*** recall; Jaccard collapses to ≈0.05 on any
degradation (single cluster shatters), containment holds 1.0 to 90 % (artifact fixed).

**Figure 5. The 2bRAD enzyme set is a resolution/cost knob (14 species).** Enzyme ladder 1/2/4/8/14 across all 14 resolvable sim-pool species. **(A)** median precision (1.0 throughout) and recall (0.50→1.00) vs enzyme number, faint lines per species. **(B)** strain-specific marker yield vs enzyme number (median 976→12 647). Cluster count is invariant. Single-enzyme BcgI enables native BcgI 2bRAD-M.

**Figure 6. Native BcgI 2bRAD strain-level identification and abundance across four ATCC DNA mocks.** Native BcgI reads were profiled against the primary 28-species,
164-genome containment tree (`--min-abundance 0 --min-coverage 0.2`), with one Strain2bScan row per sample.
Right-hand columns show AUPR, precision/recall/F1 at abundance >= 1e-4, Bray-Curtis similarity and L2
similarity. At 1e-4, F1 was 0.952 at 0.1 ng and 0.625 at 0.01 ng; it remained 1.0 at 99% host. MSA-1003
showed the 1x marker-depth noise floor, and MSA-1005/1007 retained recall and AUPR = 1.0.

**Figure 7. Real saliva: strain profiles discriminate individuals and are temporally stable.**
**(A–C)** : PCoA (species vs strain) coloured by subject,
and per-species strain-level subject R² (8 subjects × 4 timepoints; strain R² 0.757 versus species 0.755,
p 2e-4; *Neisseria subflava* R² 0.808; 13/13 testable species significant). **(D–E)**
: within- vs between-subject strain distance (0.327 vs 0.696,
p 8.05e-19) and leave-one-timepoint-out host-ID accuracy (100 % strain vs 78.1 % species).

**Figure 8. Native 2bRAD contains all callable shotgun strains and adds candidate low-abundance calls.**
; paired shotgun↔2bRAD on the same saliva samples. (A) Per sample,
strain-resolved clusters detected: 65/65 shotgun calls were also present in native BcgI 2bRAD-M (three paired samples; prefix-subsampled shotgun) plus 2bRAD-only
(128–163/sample). These are directional concordance and candidate-call results, not independent validation. (B) Community relative abundance of shared vs 2bRAD-only strains (2bRAD-only
significantly lower, p 1.2e-23).

**Figure 9. Fast, light, and scalable to whole communities.** **(A)** :
per-sample time/memory vs StrainScan (~8×/~11×). **(B)** : thread scaling of
build and profile (4.6×/5.8× to 16 threads). **(C)** : 55-species
community; measured Strain2bScan cost versus projected per-species querying, 121–146× faster.

**Figure 10. Matches or exceeds StrainScan on its own databases.** :
head-to-head on StrainScan's reference sets: *A. muciniphila*, *P. copri* (precision 1.0 both; recall
0.93 vs 0.24, 0.94 vs 0.90; ~17–23× faster, ~15–24× lighter) and near-clonal *M. tuberculosis*, where
StrainScan did not complete (>3.3 h, >25 GB) and Strain2bScan finished in 0.89 s.

**Figure 11. Systematic head-to-head on the 15-species simulated benchmark.** Both tools build databases from the same
genome pool and profile the same simulated reads, each scored in its own cluster space. **(A–C)** median
single-species precision/recall/F1 vs per-strain sequencing depth (14 species, 204 depth-matched paired
samples): both hold precision 1.0, Strain2bScan reaches full recall by 3× vs StrainScan's 10×.
**(D)** per-species database build time (log): Strain2bScan 0.7–5.1 s vs StrainScan 5–43 min (249–614×).
**(E)** per-sample profile time vs depth in the same emulated container (4–33× faster).
**(F)** multi-species profiling time per community sample; Strain2bScan's one digest-once pass vs
StrainScan's sum over 14 per-species runs (46–105×).

**Figure 12. Primary WMS mock comparison: Strain2bScan, StrainScan and inStrain.** Whole-metagenome shotgun reads from four ATCC mocks were
profiled by the three named tools. Strain2bScan used the primary all-enzyme 164-genome tree; StrainScan used
its per-species databases; inStrain used a dereplicated 98%-ANI reference. Each tool was scored in its own
0.95-similarity cluster space, and right-hand columns use abundance >= 1e-4 for P/R/F1. At 99% host,
Strain2bScan retained F1 = 1.0 and Bray-Curtis dissimilarity 0.302, whereas StrainScan had F1 = 0.095 and
inStrain F1 = 0.333. MSA-1005/1007 show that Strain2bScan retained recall and AUPR = 1.0 but had F1 =
0.60-0.75 because trace-abundance clusters crossed the threshold. Optional trace-gap and port-layer runs are
archived as sensitivity outputs but are not shown.

### Supplementary figures and tables

**Figure S3. Cost of unified-database expansion (20- vs 28-species combined tree).** Same row-per-sample
layout as Fig 6/12, WMS modality: Strain2bScan on the 20-species (`120`) vs 28-species (`164`) combined
tree, for the mocks present in both (MSA-1002, MSA-1003). MSA-1002 (even) is identical on both trees;
DB expansion costs nothing when strains are well above the marker-depth floor, even under host. MSA-1003
(staggered) drops single-threshold precision to 0.74–0.83 on the 28-species tree vs 1.00 on the
20-species tree, while AUPR stays ≈ 0.96 either way; the cost is confined to the ~1× noise floor.

**Figure S4. inStrain requires a dereplicated reference (MSA-1002 shotgun).**
 (generated 2026-10-02 from  and the non-dereplicated reference control). inStrain on the
non-dereplicated 164-genome reference vs a 98 %-ANI dereplicated reference (), across
the MSA-1002 host ladder. The non-dereplicated reference multi-maps reads across the > 95 %-ANI decoys
and inflates false positives (precision 0.29 at 0 % host); dereplication restores precision to 1.0. This
is why Fig 12 reports inStrain on its documented dereplicated reference.

(The former Fig S1 (per-species rank–rank scatter) is now **Fig 2 panel B**.)

**Figure S1. ATCC MSA-1002 DNA-input titration (native BcgI 2bRAD-M).** :
precision 1.0 and full recall to 0.1 ng input.

**Figure S2. Layer-1 gate calibration on the 55-species panel.** ;
species precision/recall vs the marker floor at normal and low depth.

**Table S3. Exploratory clinical oral cohort profiling.** : 4 clinical
oral 2bRAD-M samples, 15–17 species and 115–158 strain calls each, ~2 s/sample (no case/control labels
public). Doc: `docs/clinical_oral_exploratory.md`.

**Table S4. Genome quality-control for the 16S/2bRAD motivation panel.** ;
assembly level, CheckM completeness/contamination, contig count and length per genome; high-quality flag.

### Numbered directory manifest

 records the source figure or combined panels, PNG/PDF names and
SHA-256 checksums for every draft-backed numbered figure. The current set is Figures 1–12 and S1–S4.
Fig07 was regenerated from corrected `sample_fraction` saliva metrics. FigS4 is the inStrain
dereplicated-reference control; the obsolete `FigS4_sim_tracegap` and `FigS5_sim_fig12_style` names
were removed from the draft-backed numbered set.

## Tables

**Table 1. Strain2bScan vs StrainScan on the 15-species simulated benchmark, per species.** Single-species accuracy is the median over depth-matched paired samples (2/3/5-strain mixtures across the 0.5–10× ladder), each tool scored in its own cluster space. Database build cost is per species (Strain2bScan native, arm64; StrainScan `linux/amd64` under emulation). n = genomes in the pool.

| Species | n | S2B P | S2B R | S2B F1 | SS P | SS R | SS F1 | S2B build | SS build |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| *A. muciniphila* | 50 | 1.000 | 0.800 | 0.889 | 0.833 | 0.667 | 0.667 | 2.8 s / 0.16 GB | 17 min / 15 GB |
| *C. difficile* | 47 | 1.000 | 0.800 | 0.889 | 0.800 | 0.667 | 0.667 | 3.2 s / 0.28 GB | 30 min / 14 GB |
| *C. acnes* | 43 | 1.000 | 0.667 | 0.800 | 1.000 | 0.667 | 0.800 | 1.5 s / 0.16 GB | 11 min / 8 GB |
| *E. coli* | 50 | 1.000 | 1.000 | 1.000 | 1.000 | 0.600 | 0.667 | 4.3 s / 0.33 GB | 43 min / 28 GB |
| *F. nucleatum* | 25 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 0.7 s / 0.10 GB | 5 min / 8 GB |
| *L. plantarum* | 50 | 1.000 | 1.000 | 1.000 | 1.000 | 0.600 | 0.750 | 5.0 s / 0.23 GB | 24 min / 17 GB |
| *M. tuberculosis* | 29 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.6 s / 0.25 GB | 16 min / 15 GB |
| *P. dorei* | 15 | 1.000 | 0.667 | 0.800 | 0.833 | 0.583 | 0.667 | 1.0 s / 0.25 GB | 5 min / 13 GB |
| *P. gingivalis* | 43 | 1.000 | 0.667 | 0.800 | 1.000 | 0.667 | 0.800 | 2.3 s / 0.15 GB | 19 min / 10 GB |
| *P. copri* | 19 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.8 s / 0.18 GB | 8 min / 25 GB |
| *S. enterica* | 45 | 1.000 | 1.000 | 1.000 | 0.833 | 0.800 | 0.800 | 4.2 s / 0.36 GB | 34 min / 16 GB |
| *S. aureus* | 48 | 1.000 | 0.667 | 0.800 | 1.000 | 0.500 | 0.667 | 1.7 s / 0.16 GB | 17 min / 12 GB |
| *S. epidermidis* | 50 | 1.000 | 0.667 | 0.800 | 1.000 | 0.600 | 0.750 | 1.6 s / 0.13 GB | 15 min / 8 GB |
| *S. pneumoniae* | 50 | 1.000 | 0.667 | 0.800 | 1.000 | 0.600 | 0.750 | 1.5 s / 0.15 GB | 11 min / 11 GB |

**Median (14 resolvable species):** Strain2bScan P 1.00 / R 0.80 / F1 0.89; StrainScan P 1.00 / R 0.67 / F1 0.75. Build speed-up 249–614×; build memory 43–138× lighter.

**Table 2. Accuracy and per-sample cost vs sequencing depth (single-species, 14 species, 204 paired samples; medians).** Profile time/memory in the same emulated container (Strain2bScan `linux/amd64` vs StrainScan).

| Depth (×) | n | S2B P/R/F1 | SS P/R/F1 | S2B time/mem | SS time/mem |
|--:|--:|:--:|:--:|:--:|:--:|
| 0.5 | 40 | 1.000/0.500/0.667 | 1.000/0.667/0.800 | 0.16 s / 21 MB | 5.3 s / 831 MB |
| 1 | 41 | 1.000/0.667/0.800 | 1.000/0.333/0.500 | 0.24 s / 30 MB | 3.5 s / 831 MB |
| 3 | 41 | 1.000/1.000/1.000 | 1.000/0.600/0.750 | 0.56 s / 60 MB | 4.2 s / 831 MB |
| 5 | 41 | 1.000/1.000/1.000 | 1.000/0.667/0.800 | 0.89 s / 91 MB | 4.8 s / 831 MB |
| 10 | 41 | 1.000/1.000/1.000 | 1.000/1.000/0.889 | 1.66 s / 161 MB | 6.9 s / 832 MB |

**Table 3. Multi-species community profiling (4 samples/depth, matched to the 14 species with StrainScan databases; medians).** Strain2bScan profiles each community in one digest-once pass; StrainScan (no multi-species mode) profiles once per species, so its cost is the sum over 14 databases.

| Community depth | S2B P/R/F1 | SS P/R/F1 | S2B time/mem | SS time/mem |
|---|:--:|:--:|:--:|:--:|
| low | 0.926/0.678/0.782 | 0.898/0.767/0.827 | 1.0 s / 311 MB | 100 s / 1112 MB |
| med | 0.863/0.853/0.869 | 0.911/0.856/0.895 | 4.3 s / 670 MB | 228 s / 1696 MB |
| high | 0.773/0.872/0.819 | 0.814/0.972/0.895 | 8.7 s / 1119 MB | 398 s / 2028 MB |

**Table 4. Public real-metagenome application subsets profiled with Strain2bScan.** These were
informative application subsets rather than complete cohort analyses. The generic panel was the legacy
20-species MSA database; calls used all-enzyme mode with `--min-species-markers 20
--min-species-detect 2 --min-support 2 --min-coverage 0.01 --min-abundance 0`.

| Project | Design and selected libraries | Calls | Observation |
|---|---|--:|---|
| PRJNA288562 | Pregnancy subject T23; saliva, vaginal swab and distal gut at GD84 and GD273 (6 WGS libraries) | 34 calls in 5/6 libraries | Saliva was cluster-rich; four *Neisseria* and four *Schaalia* clusters were shared across timepoints. Gut *B. adolescentis* C4 persisted, whereas three other *Bifidobacterium* clusters and three *E. coli* clusters were GD84-only. |
| PRJNA1517970 | Preterm-birth vaginal/meconium subset plus extraction blank (7 WGS libraries) | 0 calls with the MSA panel; 0 calls with a 13-species body-site panel | The blank had 340 distinct markers, and selected libraries had 7,394–86,980. At a one-marker diagnostic gate, one meconium library showed only three *C. acnes* markers, supporting specificity while also indicating incomplete niche coverage. |
| PRJNA1191223 | Preterm infant P08 stool at W1, W2 and W3 (3 WGS libraries) | 9 generic-panel calls | Generic panel showed weekly turnover: one *S. aureus* cluster (W1), two *E. coli* clusters (W2), and six *E. faecalis* clusters (W3). |
| PRJNA1191225 | Six preterm-infant isolate WGS read sets used for cohort-panel validation | 6/6 expected self-calls | The cohort-specific panel recovered all expected *E. coli* clusters or Bifidobacterium species units. |

**Table 5. Isolate assemblies and cohort-specific panel validation.** Isolates were assembled with
SPAdes `--isolate`. *E. coli* assemblies were clustered at 0.95 similarity; each Bifidobacterium
species was built as a single-genome database. All isolate read sets recovered the expected panel unit.

| Isolate | Species | Contigs | Assembly (Mb) | N50 (kb) | Panel unit | Panel markers | Self-call |
|---|---|--:|--:|--:|---|--:|:--:|
| LHCA45 | *Escherichia coli* | 185 | 4.96 | 240.0 | C0 | 7,786 unique | ✓ |
| LHCA56 | *Escherichia coli* | 448 | 5.15 | 222.6 | C1 | 5,466 unique | ✓ |
| LHCA72 | *Escherichia coli* | 457 | 5.04 | 159.2 | C2 | 5,665 unique | ✓ |
| LHCA43 | *Bifidobacterium longum* | 164 | 2.36 | 55.6 | LHCA43 | 18,426 | ✓ |
| LHCA81 | *Bifidobacterium breve* | 90 | 2.39 | 228.3 | LHCA81 | 19,214 | ✓ |
| LHCA82 | *Bifidobacterium bifidum* | 180 | 2.30 | 75.6 | LHCA82 | 17,985 | ✓ |

**Table 6. P08 weekly stool profiling with the cohort-specific isolate panel.** Breadth is database
coverage; abundance is within-species abundance in Strain2bScan output.

| Week | Species | Panel unit | Breadth | Depth (×) | Abundance | Fraction of sample |
|--:|---|---|--:|--:|--:|--:|
| W1 | *Bifidobacterium bifidum* | LHCA82 | 0.764 | 236.7 | 1.000 | 0.6456 |
| W2 | *Bifidobacterium bifidum* | LHCA82 | 0.773 | 213.9 | 1.000 | 0.5834 |
| W2 | *Escherichia coli* | C0\|C1\|C2 | 0.026 | 0.029 | 1.000 | 0.000140 |
| W3 | *Bifidobacterium bifidum* | LHCA82 | 0.792 | 238.2 | 1.000 | 0.6099 |

The generic MSA panel detected the W2 *E. coli* signal as two low-abundance generic clusters
(fraction 6.1 × 10⁻⁵ and 5.7 × 10⁻⁵) but did not contain the persistent *B. bifidum* LHCA82 unit. This
contrast illustrates why body-site- or cohort-specific panels are needed for biological interpretation.


**Table 7. Leave-one-isolate-out behaviour for the three *Escherichia coli* panel isolates.** For each test, the
held-out isolate was absent from the database. Reads were profiled against a panel containing the other
two *E. coli* assemblies and the three Bifidobacterium species controls.

| Held-out isolate | *E. coli* genomes in panel | Observed *E. coli* unit | Breadth | Depth (×) | Fraction of sample | Interpretation |
|---|--:|---|--:|--:|--:|---|
| LHCA45 | 2 | C0\|C1 | 0.799 | 7.87 | 0.612 | Reads from the absent isolate were assigned to merged relatives rather than creating a false third cluster. |
| LHCA56 | 2 | C0\|C1 | 0.911 | 10.55 | 0.649 | Same conspecific-assignment behaviour. |
| LHCA72 | 2 | C0\|C1 | 0.891 | 7.25 | 0.729 | Same conspecific-assignment behaviour. |

This test shows that closed panels can misattribute a truly absent conspecific strain to near relatives.
It therefore supports the paper's conservative treatment of the low-coverage W2 *E. coli* signal as an
unresolved `C0|C1|C2` unit rather than assigning it to one isolate.

## References

*Working bibliography compiled from the tools and methods cited in the text. Author-facing note:
verify exact volume/issue/page numbers and publication years against the primary sources before
submission; DOIs/accession numbers should be added.*

1. Liao H, Ji Y, Sun Y. **StrainScan: a highly accurate and sensitive computational tool for strain-level
   detection of the microbiome.** *Microbiome* 2023; 11:186.
2. Sun Z, Huang S, Zhu P, *et al.* **Species-resolved sequencing of low-biomass or degraded microbiomes
   using 2bRAD-M.** *Genome Biology* 2022; 23:36.
3. Wang S, Meyer E, McKay JK, Matz MV. **2b-RAD: a simple and flexible method for genome-wide
   genotyping.** *Nature Methods* 2012; 9(8):808–810.
4. Truong DT, Tett A, Pasolli E, Huttenhower C, Segata N. **Microbial strain-level population structure
   and genetic diversity from metagenomes.** *Genome Research* 2017; 27(4):626–638.
5. van Dijk LR, Walker BJ, Straub TJ, *et al.* **StrainGE: a toolkit to track and characterize
   low-abundance strains in complex microbial communities.** *Genome Biology* 2022; 23:74.
6. Shaw J, Yu YW. **Rapid species-level metagenome profiling and containment estimation with sylph.**
   *Nature Biotechnology* 2024 (advance online).
7. Ondov BD, Treangen TJ, Melsted P, *et al.* **Mash: fast genome and metagenome distance estimation using
   MinHash.** *Genome Biology* 2016; 17:132.
8. Ondov BD, Starrett GJ, Sappington A, *et al.* **Mash Screen: high-throughput sequence containment
   estimation for genome discovery.** *Genome Biology* 2019; 20:232.
9. Brown CT, Irber L. **sourmash: a library for MinHash sketching of DNA.** *Journal of Open Source
   Software* 2016; 1(5):27.
10. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. **CheckM: assessing the quality of
    microbial genomes recovered from isolates, single cells, and metagenomes.** *Genome Research*
    2015; 25(7):1043–1055.
11. Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. **CheckM2: a rapid, scalable and accurate tool for
    assessing microbial genome quality using machine learning.** *Nature Methods* 2023; 20(8):1203–1212.
12. Orakov A, Fullam A, Coelho LP, *et al.* **GUNC: detection of chimerism and contamination in
    prokaryotic genomes.** *Genome Biology* 2021; 22:178.
13. Bowers RM, Kyrpides NC, Stepanauskas R, *et al.* **Minimum information about a single amplified genome
    (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea.** *Nature Biotechnology*
    2017; 35(8):725–731.
14. Huang W, Li L, Myers JR, Marth GT. **ART: a next-generation sequencing read simulator.**
    *Bioinformatics* 2012; 28(4):593–594.
15. Marçais G, Kingsford C. **A fast, lock-free approach for efficient parallel counting of occurrences of
    k-mers.** *Bioinformatics* 2011; 27(6):764–770.
16. Baker DN, Langmead B. **Dashing: fast and accurate genomic distances with HyperLogLog.** *Genome
    Biology* 2019; 20:265.
17. Minkin I, Medvedev P. **Scalable multiple whole-genome alignment and locally collinear block
    construction with SibeliaZ.** *Nature Communications* 2020; 11:6327.
18. Eddy SR. **Accelerated profile HMM searches.** *PLoS Computational Biology* 2011; 7(10):e1002195.
19. Seemann T. **Barrnap: basic rapid ribosomal RNA predictor.** Software, https://github.com/tseemann/barrnap.
20. Anderson MJ. **A new method for non-parametric multivariate analysis of variance.** *Austral Ecology*
    2001; 26(1):32–46. (PERMANOVA)
