# Strain2bScan: strain-level profiling from 2bRAD markers across microbiomes

*Thesis chapter. Extended background, technical methods and discussion relative to the manuscript version (`full_manuscript.md`).*

---

## Summary

Strain-level variation underlies clinically important microbial phenotypes, but 16S rRNA surveys usually
resolve species rather than strains, whereas shotgun profilers can scale poorly and lose signal in
host-dominated, low-biomass samples. We present Strain2bScan, a Rust profiler that clusters conspecific
reference genomes and scores unique 25–33 bp type-IIB restriction markers. It profiles native BcgI
2bRAD-M libraries and in-silico-digested shotgun reads in one marker framework. Across 15 species,
2bRAD-marker distances tracked whole-genome strain distances more closely than 16S distances (median
Spearman ρ = 0.94 versus 0.36). Native BcgI profiling retained F1 = 1.0 at an abundance threshold of 10⁻⁴
in the 99%-human-DNA mock. In the primary shotgun comparison, Strain2bScan was the only tested tool that
preserved detection and abundance at 99% host. In a
run-level StrainScan rerun, accuracy was comparable: species-median precision was 1.0 for both tools,
median recall/F1 were 0.733/0.844 versus 0.667/0.800, and paired-difference intervals included zero.
Archived cross-environment runs reduced database-construction cost by 249–614× in wall time and 43–138×
in memory. A
supplementary StrainGST rerun used a different 0.90-reference space; StrainGR was not run. In eight
saliva subjects, a leave-one-timepoint-out classifier reached 100% with strain features versus 78.1%
with species features; the paired accuracy-difference interval included zero. This was a small
within-subject classification result, not independent strain-accuracy validation. Accuracy depends on
reference completeness and niche-matched panel design.

**Availability.** Rust source: https://github.com/HuangShiLab/Strain2bScan. Public reads: PRJNA1131785,
PRJNA288562, PRJNA1517970, PRJNA1191223 and PRJNA1191225. Derived outputs:
https://github.com/HuangShiLab/Strain2bScan-paper.

## 1. Introduction

Bacteria of a single named species are not interchangeable. Isolates that a 16S survey; or even a
species-level shotgun profile; would report as one taxon can differ by tens of thousands of single-
nucleotide variants, by the presence or absence of whole genomic islands, plasmids and prophages, and,
consequently, in phenotypes that matter directly to health: virulence, antibiotic-resistance carriage,
toxin production, metabolic capacity and the ability to colonise a particular host. *Escherichia coli*
ranges from harmless commensal to enterohaemorrhagic pathogen; *Cutibacterium acnes* phylotypes partition
with healthy versus acne-associated skin; *Klebsiella pneumoniae* lineages differ sharply in
carbapenem resistance and hypervirulence. Resolving *which* strains are present in a sample, and at what
relative abundance, is therefore a central problem for microbiome science and for clinical metagenomics ;
and it is a fundamentally different problem from cataloguing which species are present.

#### The limits of 16S and the promise of shotgun strain typing

The 16S rRNA gene, the historical workhorse of community surveys, is structurally unable to resolve this
variation. The gene is short, highly conserved, and often present in multiple, non-identical copies per
genome; its between-strain distances are essentially uncorrelated with genome-wide divergence. This
chapter quantifies that failure directly (Fig 2): across fifteen species, the rank correlation between
16S distance and whole-genome distance has a median of only 0.36, with several species indistinguishable
from zero. 16S resolves species, not strains.

Whole-genome shotgun metagenomics does carry strain information, and over the past decade a family of
tools has learned to extract it. They fall into a few families. *Reference k-mer methods*; StrainScan,
StrainGE; index the k-mer content of a curated set of reference genomes per species and assign a sample's
k-mers to the best-supported combination of reference strains (or clusters of near-identical strains).
*Marker-gene methods* (StrainPhlAn) call SNPs within a fixed set of clade-specific marker genes and
reconstruct per-sample haplotypes. *Sketch/containment estimators*; sylph, and the Mash/sourmash
lineage they build on; estimate the containment of reference genomes in a sample from subsampled k-mer
sketches, trading per-SNP resolution for speed. These tools have made strain-level profiling routine for
individual samples and individual species, and each embodies a different point on the accuracy/speed/
generality trade-off surface.

#### Two settings where shotgun strain typing struggles

Two obstacles keep shotgun-based strain profiling out of two settings that are becoming central to the
field.

**Scale.** Modern microbiome studies span hundreds to thousands of samples and increasingly aim to
resolve strains across the dozens-to-hundreds of species that co-occur in a community, not one species at
a time. Full k-mer methods index and query the entire k-mer content of every reference genome; the
dominant per-sample cost is counting the k-mers of the sample, and (critically) that cost is *paid
again for every species database queried*, because there is no shared per-sample representation across
species. For *S* species across *N* samples the work scales as *N × S ×* (k-mer count + search): the
species dimension multiplies rather than amortises. In practice this means hours-to-days of compute and
hundreds of megabytes to gigabytes of resident memory for a community-scale, multi-species survey; a
ceiling that this chapter shows is not intrinsic to the strain-typing problem but to the full-k-mer
representation.

**Low biomass and host contamination.** Many of the clinical niches where strain resolution would be most
valuable yield very little microbial DNA against an overwhelming human background. Saliva and other oral
sites, tumour and formalin-fixed paraffin-embedded (FFPE) tissue, and skin routinely present ≥90–99 %
host DNA. Ordinary shotgun sequencing spends its reads in proportion to DNA abundance, so at 99 % host a
sequencing run devotes ~99 % of its reads to the human genome; the informative bacterial fraction, and
with it the rarer strains, is simply not sequenced deeply enough to be resolved. Strain resolution
collapses in exactly the samples where it would matter most, and no amount of downstream computation can
recover reads that were never generated.

#### Reduced-representation 2bRAD sequencing

Reduced-representation sequencing offers a route around both obstacles. Type-IIB restriction enzymes
(the basis of the "2bRAD" method) cut on *both* sides of a short, degenerate recognition site, excising a
fixed-length fragment (the 2bRAD tag, here 25–33 bp) at every occurrence of the site in a genome. The
result is a sparse, reproducible, genome-wide sample of roughly 1–2 % of the genome. Because the tag set
is defined by the recognition sequence rather than by abundance, the *same* loci are recovered from any
genome that contains them, making tags directly comparable across samples and reference genomes.

Two properties make this attractive for the two problem settings above. First, the tag set is 50–100×
smaller than the full k-mer set of a genome while remaining genome-wide and taxonomically structured ;
precisely the low-redundancy, informative marker set that a strain-resolution framework needs, since such
frameworks never use the whole genome but only the *unique* markers that distinguish a strain or a cluster
of near-identical strains. Second, and decisively for the low-biomass problem, 2bRAD can be realised as a
*wet-lab* protocol (2bRAD-M and its faster variant Fast2bRAD-M): the reduction to informative tags happens
at the bench, during library preparation, *before* host DNA can swamp the sample. A native 2bRAD library
of a 99 %-host sample concentrates sequencing on bacterial tags rather than on the human genome, giving
accurate profiles from picogram inputs and heavily contaminated or degraded material.

To date, however, 2bRAD has been used only for **species**-level profiling. The tags carry strain-level
signal (this chapter demonstrates that they track whole-genome divergence where 16S does not) but no
method had yet exploited that signal for strain resolution, and it was not obvious *a priori* that a marker
set 50–100× sparser than a full k-mer index would retain enough discriminating, single-copy loci to
separate near-identical strains at realistic sequencing depths.

#### Contribution of this chapter

This chapter presents **Strain2bScan**, a tool that ports the two-layer StrainScan resolution framework ;
within-species clustering into a search structure, followed by unique-marker detection and abundance
estimation; onto 2bRAD tags, implemented in dependency-free Rust for speed and portability. Its tag
lengths and recognition patterns match Fast2bRAD-M / `2bRADExtraction.pl` exactly, so its markers are
interoperable with the Fast2bRAD-M species layer and the two operate in one shared tag space. Uniquely,
Strain2bScan accepts **two input modes**, one for each of the obstacles above:

1. **Native BcgI 2bRAD experimental libraries**, whose reads *are* the tags; enabling, for the first
   time, strain-level analysis of low-biomass, high-host microbiomes. Because the reduction happens at the
   bench, Strain2bScan holds precision 1.0 and full strain recall at 99 % host DNA, where in-silico
   digestion of ordinary shotgun of the same material loses most of its strains; on real saliva it
   resolves individual-specific, temporally stable strain signatures and recovers low-abundance strains
   host-limited shotgun cannot reach.

2. **In-silico digestion of conventional shotgun metagenomes**; enabling community-scale strain
   profiling. The sample is digested **once** into a shared tag representation and matched against every
   per-species database, so the marginal cost of an additional species is a hash-set lookup rather than a
   re-count, and per-sample cost is essentially independent of the number of species. This turns the
   *N × S ×* (count + search) scaling of full-k-mer methods into *N ×* (digest + *S·ε*), an *S*-fold
   structural saving that this chapter measures directly.

The two modes are shown to agree; in-silico and native digestion of the same material recover the same
strains; so a single tool, operating on one 2bRAD tag space, spans both the low-biomass clinical regime
and the cohort-scale regime. The chapter first establishes the shared foundation (the 2bRAD-versus-16S
motivation, core accuracy, sensitivity, and the effect of reference-genome quality), then develops the two
input-mode pillars in turn, and closes with a systematic head-to-head against StrainScan on a common
simulated benchmark that isolates the accuracy and the cost of the two approaches under identical inputs.

## 2. Results

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

#### Real saliva: individual-specific, within-day strain signatures (Fig 7)

We profiled native BcgI 2bRAD saliva from 8 subjects sampled at four times of day (32 libraries) against a
19-species oral reference panel. Subject PERMANOVA R² was similar for strain and species profiles (0.757
versus 0.755; both p = 2 × 10⁻⁴). Leave-one-timepoint-out nearest-neighbour classification identified the
host with 100% accuracy from strain features versus 78.1% from species features (Fig 7), but this was a
small within-subject classification result across eight subjects and does not constitute independent
strain-identity validation. The strongest
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

#### Matches or exceeds StrainScan on its own databases (Supplementary Fig S5)

On StrainScan-curated reference sets [1], both tools reached precision 1.0 for *A. muciniphila* and
*P. copri* (Supplementary Fig S5). Strain2bScan matched or exceeded recall (0.93 versus 0.24 and 0.94
versus 0.90), was 17–23×
faster and 15–24× lighter, and completed near-clonal *M. tuberculosis* in 0.89 s, whereas StrainScan did not
complete. Its low *M. tuberculosis* recall reflects the resolution limit of a panel that collapses to five
0.95-similarity clusters.

The shotgun mode was then stress-tested in ATCC mocks and compared with saliva. It preserved detection and
abundance under host contamination in the primary MSA-1002 comparison (Fig 11), and its saliva calls were
contained in the native-2bRAD call set (Fig 8). These are concordance and stress-test results, not
independent validations of strain identity.

#### Systematic head-to-head on a 15-species simulated benchmark (Fig 10, Tables 1–3)

On a common 15-species simulation pool, both tools built databases from the same genomes and profiled the
same simulated reads. Each tool was scored in its own 0.95-similarity cluster space. The reproducible
comparator rerun comprised 225 matched single-species entries; 204 completed and were paired to the
corresponding Strain2bScan runs, whereas 21 ended without calls, including all *Salmonella enterica*
entries. Across the 204 completed runs from 14 species, species-median precision was 1.0 for both tools;
species-median recall and F1 were 0.733 and 0.844 for Strain2bScan versus 0.667 and 0.800 for StrainScan
(Fig 10, Table 1). In a species-cluster bootstrap, the paired mean differences in recall (−0.014; 95% CI
−0.065 to 0.033) and F1 (−0.023; 95% CI −0.062 to 0.013) included zero. Strain2bScan reached full median
recall at 3× and remained there, whereas StrainScan's depth profile was non-monotonic (Table 2).
Multi-species F1 was lower for Strain2bScan at medium and high depth but higher at low depth (Table 3).

The cost differences were largest during database construction: Strain2bScan used 0.7–5.1 s and 0.1–0.4 GB
per species, versus 5–43 min and 8–28 GB for StrainScan, a 249–614× runtime and 43–138× memory advantage
(Table 1). In the same emulated container, per-sample profiling was 4–33× faster and 5–39× lighter
(Table 2). Because StrainScan lacks a multi-species mode, its community cost was the sum of 14 runs
(100–398 s), versus 1–9 s for one digest-once pass (46–105× faster). StrainScan also failed to build
*Klebsiella pneumoniae* in the archived timing run; in the reproducible accuracy rerun, however, StrainScan
built and profiled this species. This discrepancy underscores resource sensitivity. StrainScan build
times were obtained under linux/amd64 emulation and are upper bounds.

#### Strain-level profiling on shotgun, and the advantage under host contamination (Fig 11)

We compared Strain2bScan, StrainScan v1.0.14 [1] and inStrain 1.10.0 [8] on shotgun reads from the same
four ATCC mocks. Strain2bScan used the all-enzyme 164-genome tree; StrainScan used per-species databases;
inStrain used a dereplicated 98%-ANI reference, as its documentation requires. Each tool was scored in its own
0.95-similarity cluster space. StrainScan and inStrain were therefore the executed strain-resolved
comparators. A supplementary StrainGST rerun [2] was added separately because it uses StrainGE-specific
0.90-reference databases and was not part of the frozen primary configuration; StrainGR was not run
(Tables 8–9).

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
three *C. acnes* markers in the largest meconium library. This negative result therefore localises the
failure to sparse strain-marker input and incomplete niche coverage rather than runtime failure: even the
largest library provided only three markers for a candidate species, far below any strain-resolution gate,
and the blank supports specificity at low marker input. These low-biomass cohorts require body-site-complete
pangenome panels and deeper marker coverage before absence calls can be interpreted biologically.

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
inferred from a nearest-relative call (Table 7).

#### Comparator scope

The primary strain-resolved comparisons used StrainScan [1] and inStrain [8] because both operate directly
on the frozen shotgun benchmark reads and could be scored under the frozen cluster-space or read-level
workflow. StrainScan is the closest methodological comparator because it uses reference-guided
strain-cluster resolution. inStrain represents widely used read-level shotgun profiling. A supplementary
StrainGST rerun covers the related reference-guided StrainGE workflow [2]; because StrainGE's documented
input is WMS rather than native BcgI 2bRAD, and StrainGR variant calling was not run, this rerun is
reported separately rather than pooled with the primary same-assay comparisons (Table 8).

#### StrainGST rerun on simulations and ATCC WMS mocks (Table 9)

StrainGST completed all 237 simulation runs and all four ATCC WMS mocks. On the 225 matched
single-species simulations, median precision, recall and F1 were 1.000 in StrainGE's 0.90-reference
space; median cost was 18.66 s and 1.30 GB per sample. On the 12 multi-species communities, median
precision was 0.901, recall was 1.000 and F1 was 0.941; median cost was 236.98 s and 7.90 GB. This
runtime comparison is specific to the recommended species-database StrainGST workflow and includes one
sample k-merization plus searches against all relevant species databases. On the same matched simulations,
Strain2bScan medians were 0.35 s and 0.051 GB for single-species samples and 4.33 s and 0.662 GB for
communities. StrainGST therefore used 53.31× more wall time and 25.49× more peak RSS for single-species
samples, and 54.73× more wall time and 12.07× more peak RSS for multi-species communities. Accuracy was
not pooled because the two tools use different cluster spaces.

ATCC behaviour was mock dependent. At the 1e-4 abundance threshold, StrainGST F1 was 0.700 on MSA-1002
at 99% host, 0.542 on MSA-1003, 0.476 on MSA-1005 and 0.600 on MSA-1007. Cost was 1007.22-2158.75 s and
14.36-21.52 GB peak RSS per sample, including one sample k-merization and 28 species searches. These
results are not directly pooled with the primary comparator curves because each tool is scored in its own
reference space and StrainGR was not run; they show that the reference-guided StrainGST workflow is
feasible on these shotgun mocks but does not remove the input-scope distinction from native 2bRAD
analysis.

## 3. Discussion

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

**Positioning within the strain-typing landscape.** Strain2bScan is best understood as a
*representation* change to the reference-k-mer family rather than a new inference principle: it keeps
StrainScan's two-layer logic; cluster near-identical strains, then score samples on markers unique to a
cluster; but swaps the full k-mer set for the 2bRAD tag set, which is 50–100× sparser yet still
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
could not build; while being dramatically cheaper. A fuller comparison against StrainGE, sylph and
StrainPhlAn on identical inputs is the natural next benchmark; the framework and evaluation harness built
here are designed to accommodate it.

**Interpreting the head-to-head, and its caveats.** Three points bound the interpretation of Fig 11 and
Table 1–3. First, StrainScan is Linux-x86-only and was run under emulation on the Apple-silicon test
machine; its *wall-clock* is therefore an upper bound, which is why the profiling comparison is reported in
a same-environment ratio (an emulated `linux/amd64` build of Strain2bScan) rather than native-vs-emulated,
and why the build-time comparison is framed as an order-of-magnitude structural difference rather than a
precise multiplier. Even after discounting emulation by a generous ~10×, the build-cost gap (minutes-to-
tens-of-minutes and 8–28 GB, versus seconds and ≤0.4 GB) and the memory gap remain large and are structural
;  StrainScan's peak memory is set by the full-k-mer matrix, not by emulation. Second, each tool is scored
in *its own* cluster space; because the two tools cluster the same genomes slightly differently (e.g. 27 vs
30 clusters for *E. coli*), the precision/recall values are not paired at the level of individual clusters
but at the level of "was each truth strain's cluster detected"; the only fair comparison when two tools
define clusters independently, and the honest resolution unit for short reads. Third, the recall advantage
is concentrated at shallow depth and on minor co-present strains: at 10× both tools saturate, so the
practical benefit is the ability to recover the low-coverage tail at a given sequencing budget, consistent
with the low-biomass results of Part I. The one species StrainScan could not build under emulation
(*K. pneumoniae*, 47 genomes × 5.5 Mb) is an extreme illustration of the same build-cost problem rather
than a separate failure mode.

**Why cluster resolution is the right unit; and where it can be pushed.** Both tools resolve to clusters
of near-identical strains because short reads cannot separate genomes that share almost all of their
sequence; reporting a specific strain when the data support only a cluster would be over-claiming. The
cluster-to-genome ratio is a per-species property of genuine diversity: for diverse species (*C. acnes*)
clusters are essentially single strains, whereas for near-clonal *M. tuberculosis* a handful of clusters
is the honest ceiling, and this chapter confirms that resolving to clusters costs no precision even in
low-diversity species. The clearest algorithmic route to *sub-cluster* resolution; the main limitation
shared with all short-read callers; is to layer a within-cluster step on top of the occurrence-based
detection used here: for a called cluster, a per-locus overlap or non-negative regression over the SNP-
bearing tags of its members could apportion signal among within-cluster genomes when depth allows, without
disturbing the between-cluster uniqueness logic that guarantees precision.

**Reference incompleteness: solved in part, and the remaining frontier.** The chapter shows that reference
incompleteness (not the sparsity of the tag set) is the one factor that genuinely degrades strain
identification, because an incomplete genome's markers are a subset of a complete relative's and Jaccard
splits them. The `--containment` mode addresses the *clustering* half of this cleanly (Fig 4): it keeps
subset genomes with their complete relatives and restores precision and recall to near-baseline down to
~80–90 % completeness, and it removes the near-clonal *M. tuberculosis* fragmentation artifact. What
containment cannot do is recover markers that are simply *absent*; a strain represented only by a partial
assembly; nor undo contamination that injects foreign tags; below ~70 % completeness it therefore
converges with Jaccard, and because it merges more aggressively it trades a little resolution on complete
panels (hence opt-in). Making strain identification *more resistant* to incomplete references is a concrete
programme: (i) a **completeness-aware detection gate** that scales the unique-marker floor by each genome's
estimated completeness (or gates on a *fraction* of a cluster's available markers rather than an absolute
count), so genuinely incomplete strains are not gated out; (ii) a **best-quality-representative** marker
set per cluster, defining the cluster's markers from its most-complete member; (iii) **upstream
completeness/contamination estimation and decontamination** (CheckM2, GUNC) feeding the quality filter; and
(iv) **pangenome-based imputation** of missing markers from complete conspecifics. The irreducible case; a
strain whose only reference is a low-completeness, contaminated genome; is a data limit no clustering can
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
species layer so that species and strain calls come from a *single* 2bRAD digest; the species gate taken
from a broad Fast2bRAD-M database and strain resolution from the per-species panels; is a natural and
high-value extension, particularly for the native-2bRAD clinical regime where a single library would then
yield both community composition and strain-level detail.

**Conclusion.** Reduced-representation 2bRAD markers, combined with a StrainScan-style [1] resolution
framework and a fast Rust implementation, make accurate strain-level profiling practical at a fraction of
the compute and memory of full-k-mer methods. The same tool spans two regimes: native 2bRAD-M for
strain-level analysis of low-biomass, high-host microbiomes, and in-silico-digested shotgun for strain
profiling across communities of many species and many samples.

## 4. Methods

### Overview

Strain2bScan reimplements the two-layer StrainScan strategy [1]; cluster near-identical strains,
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

The sketch was verified to behave as FracMinHash [10] requires. Marker count tracks 1 / *S* to within
sampling noise across three decades of scale (ratios 0.995–1.008 relative to the unsketched set,
on *E. coli* K-12), and the FNV-1a hash produces no collisions on the 989,962 distinct canonical
31-mers of a 1 Mb window, against 0.027 expected for a uniform 64-bit hash. The apparent
shortfall against genome length (4,523,995 markers from 4,641,622 windows) is entirely the
single-copy filter removing the 30,274 multi-copy 31-mers, matching an independent count exactly.

### Reference database construction

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
MockMetagenomes4Benchmark [16,17] (~100k read pairs each, ~12× total). (ii) A simulated multi-species
benchmark: 55 real species × ~4 strains (218 NCBI genomes) and 30 samples, each mixing strains
from twelve species at log-normal depth ≥1× (plus a low-depth variant, median 0.62×, used for
gate calibration). (iii) Cross-species mocks for *Staphylococcus
aureus* and *S. epidermidis* (60-genome panels each; 2–5 strains/sample, log-normal ≥1×,
matching the *C. acnes* design). (iv) A reference-degradation gradient in which the truth
strains' database genomes are degraded to completeness 100→50 % (with co-varying contamination
0→10 % and fragmentation), samples held fixed. Simulated reads were generated with ART [18] rather than treated as error-free. The systematic 15-species
benchmark used `art_illumina -p -l 250 -m 600 -s 150`, producing 250 bp paired-end reads with the ART
Illumina quality-error model. Earlier diagnostic datasets may have used different read configurations;
the systematic comparison and its figures are defined by this ART configuration.

**Real-data and motivation datasets.** (v) *2bRAD-vs-16S motivation* (Fig 2): 15
pathogenic/commensal species, ~50 genomes each from NCBI accession lists (ENA FASTA), **restricted to
complete/near-complete assemblies** (CheckM [15] completeness ≥ 97 %, contamination ≤ 5 %, assembly level
Complete Genome/Chromosome, consistent with high-quality MISAG/MIMAG criteria [19]; `data/genome_qc_16s_panel.tsv`). Between-strain distance was computed in
three spaces; whole-genome (bottom-3000 canonical 21-mer MinHash), 2bRAD (Strain2bScan `build` BcgI
tags) and 16S (longest gene per genome via barrnap [20] 0.9 + HMMER [21], 21-mer Jaccard); all with the
Mash [9] transform D(J) = −ln(2J/(1+J)); per species the 2bRAD and 16S pairwise vectors were correlated
(Spearman)
against the whole-genome vector, with 95 % CIs from 500 genome subsamples. (vi) *ATCC DNA mocks,
strain-level (Fig 6, Fig 11, Fig S3, Fig S4)*: four whole-cell mocks; MSA-1002 (20 strains,
even; native BcgI 2bRAD and shotgun WMS across a 0/90/95/99/99.9 % human-DNA ladder and a 1→0.001 ng
low-biomass ladder, SRA PRJNA1131785), MSA-1003 (20 strains, staggered), MSA-1005 and MSA-1007 (6 strains
each). A single unified combined tree was built from **28 species × up to 6 genomes = 164 genomes** (each
mock species = its ATCC genome + up to 5 high-quality conspecific decoys, CheckM [15] completeness ≥ 90 %,
contamination ≤ 5 %, within-species ANI 95–99.9 % to the ATCC reference by skani [22]), clustered at 0.95
similarity with `--containment`; native 2bRAD used the BcgI tree and shotgun used the all-enzyme tree.
Strain2bScan was run with `--min-abundance 0 --min-coverage 0.2`. On the shotgun samples it was compared
against **StrainScan** 1.0.14 [1] (per-species databases, `linux/amd64` container) and **inStrain** 1.10.0
[8] (Bowtie2 [23], rather than BWA-MEM [24], → `inStrain profile` against a 98 %-ANI dereplicated reference; the non-dereplicated reference is
shown as a control in Fig S4). Each tool was scored in its own 0.95-similarity cluster space
against the mock ground truth (`Ground_truth/*`, sequence abundance), reporting precision, recall, F1,
AUPR (abundance-threshold sweep, Ye et al. [25]), and Bray–Curtis and L2 similarity to the truth profile;
scorer `scripts/score_all.py`, figures `scripts/plot_figs_h.py`. The primary frozen comparators were
therefore StrainScan and inStrain. A separate StrainGST rerun [2] was added after the frozen run because it
requires StrainGE-specific 0.90-reference databases and a WMS-oriented workflow; it is described under
*StrainGST rerun*. StrainGR was not run. (vii) *Real saliva* (Fig 7, Fig 8): native BcgI
2bRAD (and paired shotgun WMS) saliva from PRJNA1131785, 8 subjects × 4 within-day timepoints, profiled
against a 19-species oral-commensal panel (up to 25 genomes/species). Strain- and species-level relative
abundances → Bray–Curtis → PERMANOVA [26] (adonis, subject/timepoint factors) and leave-one-timepoint-out
1-NN host classification; shotgun R1 (in-silico BcgI) compared to native 2bRAD calls per sample. Full
per-dataset procedures and accessions are in `docs/` (`motivation_16s.md`,
`saliva_individual_discrimination.md`, `saliva_temporal_ml.md`, `saliva_concordance.md`).

**Systematic head-to-head on a 15-species simulated benchmark (Fig 10, Table 1–3).** A common
benchmark was built from a fixed pool of 15 pathogenic/commensal species (15–50 complete/near-complete
NCBI genomes each; `figure_raw_data/sim_pool_manifest.tsv`). *Single-species* samples were generated for
every species as 2/3/5 co-present strains drawn either from the same or from different 0.95 clusters, at
per-strain coverages 0.5/1/3/5/10× with uneven abundance ratios (following StrainScan's simulation
design), 5 replicates per cell; 2 025 samples. *Multi-species* samples mixed ~18 co-present species
(one to a few strains each) across three community depth gradients; 60 samples. Reads were simulated
with ART [18] (`art_illumina -p -l 250 -m 600 -s 150`, error-modelled 250-bp paired-end reads) from the truth genomes; truth tables record each strain's
species, genome accession and 0.95-cluster assignment.

Both tools [1] built their databases from the **same genome pool** and profiled the **same reads**.
Strain2bScan databases were built with `cluster --enzyme all --similarity 0.95` and profiled with
`profile` / `multi-profile --enzyme all` (reads decompressed, R1+R2 concatenated). StrainScan (v1.0.14,
bioconda) is Linux-x86-only; it ships Dashing [27] and jellyfish [28] executables, a Python-3.7
`.so`, and an R reclustering step; so it was run inside a Docker `linux/amd64` container (QEMU emulation
on Apple Silicon; `strainscan_build`, then `strainscan -i R1 -j R2 -d DB`). Because the two tools cluster
genomes independently, **each tool was scored in its own cluster space**: predicted clusters were compared
against the truth strains mapped into that tool's clusters; for Strain2bScan via the truth `cluster`
column, for StrainScan via its `Cluster_Result/hclsMap_95_recls.txt` (report `Cluster_ID` = `C`+cluster
id); and precision/recall/F1 computed over the cluster sets per sample. Strain2bScan profiled all 2 025 +
60 samples; the reproducible StrainScan rerun used a matched subset (different-cluster mixtures, k = 2/3/5,
one replicate, all depths; near-clonal *M. tuberculosis* via its same-cluster samples). Of these 225
intended entries, 204 completed and were paired; 21 ended without calls, including all *S. enterica*
entries. The paired accuracy set therefore comprised 14 species. For communities, 12 samples were rerun and
158 of 180 possible species-by-sample final reports were generated; absent reports were treated as no
detection. StrainScan has no multi-species mode, so archived cost was the sum of wall-clock over species
databases (peak RSS as the maximum).

*StrainGST rerun.* We additionally ran StrainGST 1.3.9, the reference-search component of StrainGE [2], on
the 225 matched single-species simulations, the 12 multi-species communities and the four primary ATCC
WMS mocks. For each species panel, all genomes were k-merized with StrainGE's default k = 23, near-subset
references were removed, remaining references were clustered at Jaccard 0.90, and one StrainGST
species-level pan-genome database was created with the recommended workflow. Each sample was k-merized
once and searched against every relevant species database (maximum five iterations for single-species
samples, eight for communities and 32 for mocks). Calls used StrainGST score >= 0.02. Accuracy was scored
by mapping truth genomes and reported references through StrainGE's 0.90-reference clusters; DNA-mock
detection additionally used abundance >= 1e-4, as in the primary comparison. Wall time is the sum of one
sample k-merization and all species searches, and peak RSS is the maximum across those stages. StrainGR
was not run, so these are StrainGST identification and abundance results rather than full StrainGE
variant calling.

*Software and timing provenance.* The primary ATCC mock benchmark was frozen at Strain2bScan commit
`f26f234b817ba7772a3f1df59ce720751e9b45b9` (release-build SHA-256
`a4cf7a4043a06c55a99d6abf1e9fd312f6243aa5a5ca021bd41915636fb56b18`, Rust 1.97.0); later software commits
were not used for those mock outputs. The primary configuration, database hashes, metric definitions, and
per-run outputs are in `results/benchmark_configuration.json` and `results/mock_benchmark_f26f234/`.
Strain2bScan build and profile times for the 15-species benchmark were native arm64. To reduce the emulation
confound in profiling speed, a `linux/amd64` Strain2bScan binary was run in the same container on the same
subset, giving the same-environment ratio in Fig 10E/Table 2. StrainScan build times were obtained under
`linux/amd64` QEMU emulation and are therefore upper bounds. DB build for *K. pneumoniae* did not complete
in the archived timing run but did complete in the reproducible accuracy rerun; these results are therefore
reported separately rather than pooled. Original `scratchpad/eval` drivers were not retained, but the new
comparator rerun has a manifest, container-derived database mappings, logs, run outputs, run-level scores,
and aggregate tables in `work/strainscan_headtohead_rerun/` and `results/strainscan_rerun/`. Uncertainty was
estimated by resampling species with replacement for 20 000 replicates and reporting 2.5th and 97.5th
percentiles of paired mean Strain2bScan-minus-StrainScan differences.

**Comparison to StrainScan (curated-DB and per-sample benchmarks).** In addition to the common benchmark
above, StrainScan v1.0.14 [1] was run on its **own** reference databases (Supplementary Fig S5) and on the same *C. acnes*
per-sample profiling comparison (Fig 9A), using its low-depth modes for the depth series.

**Primary ATCC mock configuration.** The primary Strain2bScan variant was the default flat path on the
164-genome containment tree, with no trace-gap filter and no Layer-1/Layer-2 override. Native BcgI libraries
used the BcgI database and shotgun libraries used the all-enzyme database. Optional trace-gap and port-layer
runs were retained as sensitivity artifacts but are not primary evidence and are not shown in Figures 6 or
11. The primary detection threshold was abundance ≥ 10⁻⁴.

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

For cohort-specific compatibility testing, six PRJNA1191225 isolate read sets were assembled with SPAdes 4.3.0 [29]
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

### Extended algorithmic detail, data structures and complexity

The subsections above summarise the pipeline; the following give the algorithmic detail, data
structures and complexity in full, as implemented in the Rust source (`src/`).

#### Type-IIB tag model and enzyme patterns

Each type-IIB restriction enzyme recognises a short, partially degenerate motif and cleaves a fixed
distance to either side, releasing a tag of constant length. An enzyme is modelled as a triple
(*upstream gap*, *anchored pattern set*, *downstream gap*): for BcgI the excised fragment is 32 bp with a
central `CGA…TGC`-type recognition anchor and *N*-runs on either flank; the sixteen enzymes of the
Fast2bRAD-M table differ in tag length (25–33 bp) and anchor. Digestion scans every offset of a sequence
and tests the anchor set; because type-IIB sites are palindromically constrained, providing both the
forward and reverse anchor patterns lets a single left-to-right pass over one strand recover the tags that
would be produced from both strands, which halves the work and (importantly) yields exactly one canonical
marker per site rather than the two-fold-inflated set an explicit both-strand scan produces (the earlier
both-strand workaround was retired for this reason, ~3.5× fewer, correct-length markers). For conventional
shotgun input the union of a chosen enzyme set is applied (`--enzyme all` uses all sixteen), enriching the
marker yield ~*n*-fold for *n* enzymes; for native BcgI 2bRAD input the single BcgI pattern is used, because
the reads already are BcgI tags.

#### Canonicalisation and hashing

Every tag is reduced to a single 64-bit integer *marker* by (i) taking the lexicographically smaller of
the tag and its reverse complement (the *canonical* form, so a tag and its complement collide), and (ii)
hashing that canonical string with FNV-1a. Genome tags and sample-read tags pass through the identical
canonicalisation and hash, so marker values are internally consistent across the reference database and any
sample; comparison and counting are then exact integer-set operations. Using a 64-bit hash rather than the
raw 2·*L*-bit packed sequence keeps the marker width constant across enzymes (tag lengths differ) and lets
every downstream structure be a hash set or hash map keyed by `u64`; the collision probability at the
marker counts used here (10⁴–10⁶ per species) is negligible.

#### Single-copy filtering

For reference genomes only *single-copy* tags (those occurring exactly once in the genome) are retained,
following the practice of StrainScan and Fast2bRAD-M. Multi-copy tags carry copy-number rather than
presence/absence information and would bias both clustering (a repeat expansion inflates set overlap) and
abundance estimation (a repeat contributes disproportionate counts); dropping them makes every retained
marker a clean presence/absence locus. Sample reads are *not* single-copy-filtered; a marker's observed
count in a sample is the quantity of interest; but detection thresholds are stated in tag (marker) units
throughout.

#### Within-species clustering: single-linkage, union-find, and MinHash acceleration

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

#### Containment clustering for uneven-completeness panels

Jaccard penalises incompleteness. If genome B is an incomplete assembly of the same strain as complete
genome A, its marker set is approximately a subset of A's, so |A∩B| ≈ |B| but |A∪B| ≈ |A|, giving
*J* ≈ |B|/|A|; which falls below τ as soon as B is materially smaller than A, spuriously splitting the two
into different clusters. The consequence is not merely coarser clustering: the shared markers, now present
in *two* clusters, are demoted from *cluster-specific* (discriminating) to *shared-partial*
(non-discriminating), so reads from the complete strain match *both* fragments, injecting false positives
and missing calls (Fig 4A).

The optional `--containment` mode replaces Jaccard with **max-containment**,

&nbsp;&nbsp;&nbsp;&nbsp; *C*(A,B) = |A∩B| / min(|A|,|B|),

which stays ≈ 1 when one marker set is contained in the other and so keeps the incomplete genome clustered
with its complete relative; the merged cluster's marker set is the *union* of its members, so the complete
member supplies the markers the incomplete one lacks and the discriminating tags are preserved (Fig 4A).
This is the containment estimator used by Mash-screen and sourmash for genomes of uneven completeness. It
is exact for small panels; for large panels the intersection is recovered from the sketch-estimated Jaccard
*Ĵ* and the exact set sizes via |A∩B| = *Ĵ*·(|A|+|B|)/(1+*Ĵ*), then divided by min(|A|,|B|). Because
*C* ≥ *J* always, containment merges at least as aggressively as Jaccard and can coarsen clusters on
already-complete panels; it is therefore opt-in (recommended for reference sets of mixed completeness),
with the default remaining Jaccard complemented by the assembly-quality filter.

#### Marker classification and the database

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
to streaming a sample's markers through the inverted index and incrementing per-cluster counters; the
per-marker work is a single hash lookup.

#### Layer-1: which species to strain-profile

Strain markers are unique only *within* a species, so a species that is absent from a sample can be hit
spuriously by shared tags of a present relative. Species selection is therefore made on **absolute
species-specific marker evidence**, never on relative abundance (which conflates community composition with
depth). Let *total* be the number of species-specific markers a species carries (tags unique to that
species across the panel; the same tag space as the Fast2bRAD-M species layer) and *present* the subset
observed in the sample at count ≥ 2. The gate is

&nbsp;&nbsp;&nbsp;&nbsp; *resolve_gate* = max(*G*, ⌈*f* · *total*⌉),&nbsp;&nbsp;&nbsp; *detect_gate* = min(*d*, *resolve_gate*),

with an absolute floor *G* (default 200), a breadth fraction *f* (default 0) that scales the bar to each
species' panel size, and a low detection floor *d* (default 10). This produces three outcomes per species:
**strain-resolved** (*present* ≥ *resolve_gate*; Layer-2 runs), **detected but not strain-resolvable**
(*detect_gate* ≤ *present* < *resolve_gate*; reported at species level with its observed marker breadth,
no strain claim), or **absent**. The middle tier is the honest treatment of a low-abundance species (present but too faint for a strain claim) rather than a binary drop or an over-call. The breadth term *f*
is scale insurance: on the panels used here *f* = 0 (the shipped default) already gives species precision
1.0, but as a panel grows large enough for a fixed floor to be outrun by cross-species leakage, a small
*f* (≈ 0.02) restores precision by raising the bar in proportion to panel size, where large-panel leakage
concentrates.

#### Layer-2: strain detection and abundance

Within each strain-resolved species, a cluster is **called present** iff at least *N* of its unique markers
are observed at count ≥ 2 (default *N* = 10, in tag units; the full-k-mer StrainScan floor of ~1240 k-mers
is inappropriate for the ~50–100× sparser tag set). Using *only* unique markers makes detection immune to
the shared-marker cross-talk that would otherwise let a greedy set-cover over a large conspecific panel
strip shared markers and starve true strains. Each present cluster's **relative abundance** is estimated
from the *median* sample count over its detected unique markers; robust to repeat and contamination
outliers, and less prone than a joint regression over shared markers to mis-attributing signal between very
similar co-present strains; a non-negative Elastic-Net solver over the marker×cluster incidence matrix is
also provided for users who prefer a regression estimate. Calls are then filtered by a minimum coverage
fraction of their unique markers (`--min-coverage`, default 0.1; suppresses spurious detection of large,
similar clusters whose absolute unique-marker count clears the floor at a tiny coverage fraction) and a
minimum relative abundance (0.02), and renormalised. When no cluster passes, the species is reported as
detectable but not strain-resolvable at the given enzyme set.

#### Complexity and the community-scale argument

Let *n* be genomes per species, *m̄* the mean single-copy marker count, *S* species, *N* samples, and *R*
the reads per sample. **Building** a species database costs O(*n*·(genome length)) to digest, O(*n*·*m̄*)
to sketch, and O(*n*²·*k*) to cluster; dominated in practice by digestion, and independent across species
(embarrassingly parallel). **Profiling** one sample against one species costs O(*R*·*L*) to digest the
reads into markers once, plus O(#markers) hash lookups against the inverted index. The decisive point is
the *community* cost. A full-k-mer tool has no shared per-sample representation across species, so it pays
the sample-side count once *per species*: total ≈ *N × S ×* (k-mer count + search). Strain2bScan digests
each sample **once** into a shared marker multiset and matches it against every species' inverted index at a
marginal cost *ε* of a hash-set intersection: total ≈ *N × (digest + S·ε)*. Because *ε* ≪ (k-mer count),
the ratio grows with *S*; an *S*-fold structural advantage confirmed empirically at ~132× on 100 samples
of a 55-species community (Fig 9C) and reproduced in the head-to-head, where StrainScan, lacking a
multi-species mode, must run each community sample once per species and so pays 100–398 s per sample versus
1–9 s for a single Strain2bScan pass (Fig 11F).

#### Implementation and determinism

Strain2bScan is written in Rust with **no third-party dependencies**. data-parallelism (genome digestion,
sketch construction, pairwise clustering, read digestion) uses scoped `std` threads
(`STRAIN2BSCAN_THREADS`, default = all cores). All hashing is deterministic (fixed FNV-1a seed) and
clustering is order-independent (union-find over a symmetric edge set), so a given database and reads
produce identical output across runs, thread counts and platforms; verified here by cross-checking the
native arm64 binary against a `linux/amd64` build of the same source, which produced identical cluster
calls on the benchmark samples.

## 5. Figure legends

### Main figures **Figure 1. Strain2bScan overview and the two input modes.** Pipeline schematic. (A) Reference construction: type-IIB (2bRAD) digestion of reference genomes into
single-copy 25–33 bp tags → within-species clustering at 0.95 Jaccard (MinHash-accelerated) → a
cluster × marker database of species-core, cluster-specific and strain-specific markers. (B) Profiling:
a sample is digested once into canonical markers, gated on species-specific markers (Layer-1), and
strains are detected and quantified within each present species from unique markers (Layer-2). The two
input modes; in-silico digestion of conventional shotgun, or native BcgI 2bRAD-M libraries whose reads
are the tags; enter the same tag space, shared with the Fast2bRAD-M species layer. **Figure 2. 2bRAD tags track genome-wide strain distance; 16S does not.** **(A)** Per-species Spearman correlation between
whole-genome (bottom-3000 21-mer MinHash) between-strain distance and the corresponding 2bRAD-tag (blue)
or 16S rRNA (red) distance, across 15 species (complete/near-complete genomes only, CheckM ≥97 %/≤5 %;
14–50 genomes per species, shown at right); error bars are 95 % CIs over genome subsamples; genus and
species on separate label lines; 2bRAD median 0.94 vs 16S median 0.36. **(B)** 3×5 matrix of per-species
rank–rank scatters; rank of each strain pair's whole-genome distance (x) vs its 2bRAD (blue) or 16S (red)
distance (y): 2bRAD hugs the diagonal in every species, 16S forms flat rank-bands (most extreme in
*P. dorei*, *M. tuberculosis*). **Figure 3. Accurate strain profiling and depth sensitivity.** (A) precision/recall/abundance error (Bray–Curtis) on real reference panels with simulated mixtures for
*C. acnes*, *S. aureus*, *S. epidermidis* (precision 1.0 throughout). (B) detection of a single *C. acnes* strain across a coverage ladder; onset at 0.5×, matching StrainScan. **Figure 4. Reference-genome completeness controls strain-identification accuracy (all 15 species).** Sample reads held fixed while the truth strains' reference genomes are
degraded (completeness 100→50 %, with co-varying contamination + fragmentation), DB rebuilt and re-profiled.
Problem → mechanism → solution. **(A)** (containment-mechanism schematic): an incomplete
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
degradation (single cluster shatters), containment holds 1.0 to 90 % (artifact fixed). **Figure 5. The 2bRAD enzyme set is a resolution/cost knob (14 species).** Enzyme ladder 1/2/4/8/14 across all 14 resolvable sim-pool species. **(A)** median precision (1.0 throughout) and recall (0.50→1.00) vs enzyme number, faint lines per species. **(B)** strain-specific marker yield vs enzyme number (median 976→12 647). Cluster count is invariant. Single-enzyme BcgI enables native BcgI 2bRAD-M. **Figure 6. Native BcgI 2bRAD strain-level identification and abundance across four ATCC DNA mocks.** Native BcgI reads were profiled against the primary 28-species,
164-genome containment tree (`--min-abundance 0 --min-coverage 0.2`), with one Strain2bScan row per sample.
Right-hand columns show AUPR, precision/recall/F1 at abundance >= 1e-4, Bray-Curtis similarity and L2
similarity. At 1e-4, F1 was 0.952 at 0.1 ng and 0.625 at 0.01 ng; it remained 1.0 at 99% host. MSA-1003
showed the 1x marker-depth noise floor, and MSA-1005/1007 retained recall and AUPR = 1.0. **Figure 7. Real saliva: individual-specific, within-day strain signatures.**
**(A–C)** PCoA (species vs strain) coloured by subject,
and per-species strain-level subject R² (8 subjects × 4 timepoints; strain R² 0.757 versus species 0.755,
p 2e-4; *Neisseria subflava* R² 0.808; 13/13 testable species significant). **(D–E)** within- vs between-subject strain distance (0.327 vs 0.696,
p 8.05e-19) and leave-one-timepoint-out host-ID accuracy (100 % strain vs 78.1 % species) in this small
within-subject classifier. **Figure 8. Native 2bRAD contains callable shotgun strain calls and adds candidate low-abundance calls.** paired shotgun↔2bRAD on the same saliva samples. (A) Per sample,
strain-resolved clusters detected: 65/65 shotgun calls were also present in native BcgI 2bRAD-M (three paired samples; prefix-subsampled shotgun) plus 2bRAD-only
(128–163/sample). These are directional concordance and candidate-call results, not independent validation. (B) Community relative abundance of shared vs 2bRAD-only strains (2bRAD-only
significantly lower, p 1.2e-23). **Figure 9. Fast, light, and scalable to whole communities.** **(A)** per-sample time/memory vs StrainScan (~8×/~11×). **(B)** thread scaling of
build and profile (4.6×/5.8× to 16 threads). **(C)** 55-species
community; measured Strain2bScan cost versus projected per-species querying, 121–146× faster. **Figure 10. Systematic head-to-head on the 15-species simulated benchmark.** Both tools build databases from the same
genome pool and profile the same simulated reads, each scored in its own cluster space. **(A–C)** median
single-species precision/recall/F1 vs per-strain sequencing depth (204 completed paired samples across 14
species): both have median precision 1.0; Strain2bScan reaches full median recall at 3× and remains there,
whereas StrainScan is non-monotonic.
**(D)** per-species database build time (log): Strain2bScan 0.7–5.1 s vs StrainScan 5–43 min (249–614×).
**(E)** per-sample profile time vs depth in the same emulated container (4–33× faster).
**(F)** multi-species profiling time per community sample; Strain2bScan's one digest-once pass vs
StrainScan's sum over 14 per-species runs (46–105×). **Figure 11. Primary WMS mock comparison: Strain2bScan, StrainScan and inStrain.** Whole-metagenome shotgun reads from four ATCC mocks were
profiled by the three named tools. Strain2bScan used the primary all-enzyme 164-genome tree; StrainScan used
its per-species databases; inStrain used a dereplicated 98%-ANI reference. Each tool was scored in its own
0.95-similarity cluster space, and right-hand columns use abundance >= 1e-4 for P/R/F1. At 99% host,
Strain2bScan retained F1 = 1.0 and Bray-Curtis dissimilarity 0.302, whereas StrainScan had F1 = 0.095 and
inStrain F1 = 0.333. MSA-1005/1007 show that Strain2bScan retained recall and AUPR = 1.0 but had F1 =
0.60-0.75 because trace-abundance clusters crossed the threshold. Optional trace-gap and port-layer runs are
archived as sensitivity outputs but are not shown. ## Supplementary figures and tables **Figure S1. ATCC MSA-1002 DNA-input titration (native BcgI 2bRAD-M).** precision 1.0 and full recall to 0.1 ng input. **Figure S2. Layer-1 gate calibration on the 55-species panel.** species precision/recall vs the marker floor at normal and low depth. **Figure S3. Cost of unified-database expansion (20- vs 28-species combined tree).** Same row-per-sample
layout as Fig 6/11, WMS modality: Strain2bScan on the 20-species (`120`) vs 28-species (`164`) combined
tree, for the mocks present in both (MSA-1002, MSA-1003). MSA-1002 (even) is identical on both trees;
DB expansion costs nothing when strains are well above the marker-depth floor, even under host. MSA-1003
(staggered) has one lower-threshold replicate on the 28-species tree (precision 0.83), whereas all three
20-species-tree replicates remain at 1.00; AUPR stays approximately 0.95–1.00 either way; the cost is confined to the ~1× noise floor. **Figure S4. inStrain requires a dereplicated reference (MSA-1002 shotgun).** (generated 2026-10-02 from the archived metrics table and the non-dereplicated reference control). inStrain on the
non-dereplicated 164-genome reference vs a 98 %-ANI dereplicated reference (), in the MSA-1002 0%-host shotgun control. The non-dereplicated reference multi-maps reads across the > 95 %-ANI decoys
and inflates false positives (precision 0.29; recall 1.00; F1 0.46); dereplication restores precision to 0.97,
recall to 0.97 and F1 to 0.99. This
is why Fig 11 reports inStrain on its documented dereplicated reference. **Figure S5. Matches or exceeds StrainScan on its own databases.** (demoted from the former main Figure 10):
head-to-head on StrainScan's reference sets: *A. muciniphila*, *P. copri* (precision 1.0 both; recall
0.93 vs 0.24, 0.94 vs 0.90; ~17–23× faster, ~15–24× lighter) and near-clonal *M. tuberculosis*, where
StrainScan did not complete (>3.3 h, >25 GB) and Strain2bScan finished in 0.89 s.

## 6. Tables

**Table 1. Reproducible StrainScan-rerun accuracy on the 15-species simulated benchmark.** Accuracy is the median over completed run-level paired samples (different-cluster mixtures for 14 species; same-cluster mixtures for near-clonal *M. tuberculosis*; k = 2/3/5; depths 0.5–10×); each tool was scored in its own 0.95-similarity cluster space. Build cost is the archived timing benchmark (Strain2bScan native arm64; StrainScan `linux/amd64` under emulation). *S. enterica* produced no completed StrainScan calls in the rerun; *K. pneumoniae*, which did not finish in the archived timing run, was successfully built and profiled in the reproducible rerun.

| Species | Paired n | Pool genomes | S2B P | S2B R | S2B F1 | SS P | SS R | SS F1 | S2B build | SS build |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| *A. muciniphila* | 14 | 50 | 1.000 | 0.800 | 0.889 | 1.000 | 1.000 | 1.000 | 2.8 s / 0.16 GB | 17 min / 16 GB |
| *C. difficile* | 15 | 47 | 1.000 | 0.800 | 0.889 | 1.000 | 0.667 | 0.800 | 3.2 s / 0.28 GB | 30 min / 14 GB |
| *C. acnes* | 15 | 43 | 1.000 | 0.667 | 0.800 | 1.000 | 0.800 | 0.889 | 1.5 s / 0.17 GB | 11 min / 8 GB |
| *E. coli* | 15 | 50 | 1.000 | 1.000 | 1.000 | 1.000 | 0.800 | 0.889 | 4.3 s / 0.34 GB | 43 min / 28 GB |
| *F. nucleatum* | 15 | 25 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 0.7 s / 0.10 GB | 5 min / 8 GB |
| *K. pneumoniae* | 15 | 47 | 1.000 | 0.600 | 0.750 | 1.000 | 0.800 | 0.889 | 5.1 s / 0.40 GB | DNF (timing run) |
| *L. plantarum* | 15 | 50 | 1.000 | 1.000 | 1.000 | 1.000 | 0.800 | 0.889 | 5.0 s / 0.24 GB | 24 min / 18 GB |
| *M. tuberculosis* | 15 | 29 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.6 s / 0.25 GB | 16 min / 15 GB |
| *P. dorei* | 10 | 15 | 1.000 | 0.667 | 0.800 | 1.000 | 0.667 | 0.800 | 1.0 s / 0.25 GB | 5 min / 13 GB |
| *P. gingivalis* | 15 | 43 | 1.000 | 0.667 | 0.800 | 1.000 | 0.667 | 0.800 | 2.3 s / 0.15 GB | 19 min / 11 GB |
| *P. copri* | 15 | 19 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.8 s / 0.19 GB | 8 min / 26 GB |
| *S. enterica* | 0 | 45 | -- | -- | -- | -- | -- | -- | 4.2 s / 0.37 GB | 34 min / 16 GB |
| *S. aureus* | 15 | 48 | 1.000 | 0.667 | 0.800 | 1.000 | 0.667 | 0.800 | 1.7 s / 0.17 GB | 17 min / 12 GB |
| *S. epidermidis* | 15 | 50 | 1.000 | 0.667 | 0.800 | 1.000 | 0.667 | 0.800 | 1.6 s / 0.14 GB | 15 min / 8 GB |
| *S. pneumoniae* | 15 | 50 | 1.000 | 0.667 | 0.800 | 1.000 | 0.667 | 0.800 | 1.5 s / 0.16 GB | 11 min / 11 GB |

**Median across 14 species with completed paired runs:** Strain2bScan P 1.000 / R 0.733 / F1 0.844; StrainScan P 1.000 / R 0.667 / F1 0.800. The species-cluster bootstrap paired mean differences and 95% intervals are in `results/uncertainty_summary.tsv`.

**Table 2. Accuracy and archived per-sample cost vs sequencing depth (single-species rerun; medians).** Accuracy is restricted to the 204 completed run-level StrainScan pairs. Cost comes from the archived same-environment emulation subset (Strain2bScan `linux/amd64` vs StrainScan).

| Depth (×) | n | S2B P/R/F1 | SS P/R/F1 | S2B time/mem | SS time/mem |
|--:|--:|:--:|:--:|:--:|:--:|
| 0.5 | 40 | 1.000/0.500/0.667 | 1.000/0.667/0.800 | 0.16 s / 21 MB | 5.30 s / 831 MB |
| 1 | 41 | 1.000/0.667/0.800 | 1.000/1.000/1.000 | 0.24 s / 30 MB | 3.48 s / 831 MB |
| 3 | 41 | 1.000/1.000/1.000 | 1.000/0.600/0.750 | 0.56 s / 60 MB | 4.17 s / 831 MB |
| 5 | 41 | 1.000/1.000/1.000 | 1.000/0.667/0.800 | 0.89 s / 91 MB | 4.79 s / 831 MB |
| 10 | 41 | 1.000/1.000/1.000 | 1.000/1.000/1.000 | 1.66 s / 161 MB | 6.92 s / 832 MB |

Strain2bScan reached median recall 1.0 at 3× and remained there; StrainScan was non-monotonic (full median recall at 1× and 10×, but not at intermediate depths). The depth table aggregates heterogeneous strain mixtures and is not a causal depth-onset estimate.

**Table 3. Multi-species community accuracy (reproducible rerun; medians) and archived profiling cost.** Four samples per depth. Strain2bScan profiles each community in one digest-once pass; StrainScan has no multi-species mode, so archived cost is the sum over per-species runs.

| Community depth | S2B P/R/F1 | SS P/R/F1 | S2B time/mem | SS time/mem |
|---|:--:|:--:|:--:|:--:|
| low | 0.872/0.650/0.741 | 0.938/0.702/0.798 | 1.0 s / 311 MB | 100 s / 1112 MB |
| med | 0.857/0.868/0.860 | 0.954/0.903/0.927 | 4.3 s / 670 MB | 228 s / 1696 MB |
| high | 0.778/0.884/0.827 | 0.972/0.958/0.965 | 8.7 s / 1119 MB | 398 s / 2028 MB |

Of 180 species-by-community opportunities (15 databases × 12 samples), 158 StrainScan final reports were generated; absent reports were treated as no detection. The per-sample coverage audit is in `results/strainscan_rerun/multi_persample.tsv`.

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

**Table 8. Comparator and assay scope.** This table distinguishes tools executed in the benchmark from
related tools discussed for context. It is a scope statement, not a performance comparison. Native 2bRAD
denotes BcgI-derived experimental libraries; shotgun denotes conventional WMS reads, which Strain2bScan
can also process after in-silico digestion.

| Tool | Executed here | Input scope | Typical resolution | Scope and role |
|---|:--:|---|---|---|
| Strain2bScan | Yes | Native BcgI 2bRAD reads and in-silico-digested shotgun reads | Strain clusters | Method under test; one marker framework supports both reduced-representation and shotgun input. |
| StrainScan v1.0.14 | Yes | Shotgun reads against per-species reference databases | Strain clusters | Direct reference-guided comparator on shotgun reads (Fig 10–11, Supplementary Fig S5; Tables 1–3). |
| inStrain v1.10.0 | Yes | Shotgun reads aligned to a dereplicated reference | Strain and SNV | Direct read-level shotgun comparator (Fig 11 and Fig S4). |
| StrainGE / StrainGST / StrainGR | StrainGST only | Shotgun WMS reads against StrainGE 0.90-reference databases; native BcgI 2bRAD is not documented input | Strain level, reference guided | Related shotgun comparator rerun in Table 9; StrainGR was not run. |
| 2bRAD-M / Fast2bRAD-M | No | Native BcgI 2bRAD reads | Mainly species | Same restriction-assay family; provides the species-level upstream context for Strain2bScan. |
| sylph | No | Shotgun sketches | Mainly species and ANI | Context for rapid species-level profiling; not a strain-resolved benchmark here. |
| StrainPhlAn | No | Shotgun marker genes | Strain types | Context for marker-gene strain profiling; not rerun here. |

**Table 9. StrainGST (StrainGE 1.3.9) rerun accuracy and efficiency.** Simulations are the matched
`_rep1_` subset and all 12 multi-species communities. Mocks are the primary WMS replicate used in
Figure 11; MSA-1002 is the 99%-host sample. Accuracy is in StrainGE's 0.90-reference space; mock calls
use abundance >= 1e-4. Wall time is one sample k-merization plus all per-species StrainGST searches; RSS
is the maximum across those stages. StrainGR was not run.

| Dataset | n | Precision | Recall | F1 | Median wall time/sample (s) | Median peak RSS/sample (GB) |
|---|--:|--:|--:|--:|--:|--:|
| Simulated single-species | 225 | 1.000 | 1.000 | 1.000 | 18.66 | 1.30 |
| Simulated multi-species | 12 | 0.901 | 1.000 | 0.941 | 236.98 | 7.90 |
| MSA-1002, 99% host | 1 | 0.700 | 0.700 | 0.700 | 1007.22 | 14.36 |
| MSA-1003 | 1 | 0.464 | 0.650 | 0.542 | 2158.75 | 21.52 |
| MSA-1005 | 1 | 0.333 | 0.833 | 0.476 | 1784.50 | 20.10 |
| MSA-1007 | 1 | 0.429 | 1.000 | 0.600 | 1555.76 | 17.69 |

**Table S3. Exploratory native BcgI 2bRAD-M oral-library profiling.** Four libraries from the
PRJNA1131785 exploratory oral cohort were profiled against the same 19-species panel used for the saliva
analysis. No case/control labels were available in public metadata. Values are from
`results/clinical_exploratory.tsv`.

| Sample | Distinct markers | Species resolved | Strain calls | Runtime (s) | Top species |
|---|--:|--:|--:|--:|---|
| S_0313206 | 271,322 | 15 | 115 | 2.3 | *Neisseria subflava*; *Haemophilus parainfluenzae*; *Actinomyces odontolyticus* |
| S_1413231 | 288,019 | 17 | 129 | 1.8 | *Neisseria subflava*; *Streptococcus mitis*; *Actinomyces odontolyticus* |
| S_8123213 | 370,379 | 16 | 149 | 2.0 | *Actinomyces odontolyticus*; *Neisseria subflava*; *Haemophilus parainfluenzae* |
| S_8221232 | 488,180 | 17 | 158 | 1.8 | *Actinomyces odontolyticus*; *Neisseria subflava*; *Rothia mucilaginosa* |

## 7. References

1. Liao H, Ji Y, Sun Y. **High-resolution strain-level microbiome composition analysis from short reads.**
   *Microbiome* 2023; 11:183.
   doi:10.1186/s40168-023-01615-w.
2. van Dijk LR, Walker BJ, Straub TJ, *et al.* **StrainGE: a toolkit to track and characterize
   low-abundance strains in complex microbial communities.** *Genome Biology* 2022; 23:74.
   doi:10.1186/s13059-022-02630-0.
3. Truong DT, Tett A, Pasolli E, Huttenhower C, Segata N. **Microbial strain-level population structure
   and genetic diversity from metagenomes.** *Genome Research* 2017; 27(4):626–638.
4. Shaw J, Yu YW. **Rapid species-level metagenome profiling and containment estimation with sylph.**
   *Nature Biotechnology* 2024; 43(8):1348–1359. doi:10.1038/s41587-024-02412-y.
5. Wang S, Meyer E, McKay JK, Matz MV. **2b-RAD: a simple and flexible method for genome-wide
   genotyping.** *Nature Methods* 2012; 9(8):808–810.
6. Sun Z, Huang S, Zhu P, *et al.* **Species-resolved sequencing of low-biomass or degraded microbiomes
   using 2bRAD-M.** *Genome Biology* 2022; 23:36.
7. Huang S. **Fast2bRAD-M: high-performance Rust reimplementation of the 2bRAD-M microbiome profiling
    pipeline.** Software, https://github.com/HuangShiLab/Fast2bRAD-M.
8. Olm MR, Crits-Christoph A, Bouma-Gregson K, Firek BA, Morowitz MJ, Banfield JF. **inStrain profiles
    population microdiversity from metagenomic data and sensitively detects shared microbial strains.**
    *Nature Biotechnology* 2021; 39(6):727–736. doi:10.1038/s41587-020-00797-0.
9. Ondov BD, Treangen TJ, Melsted P, *et al.* **Mash: fast genome and metagenome distance estimation using
   MinHash.** *Genome Biology* 2016; 17:132.
10. Brown CT, Irber L. **sourmash: a library for MinHash sketching of DNA.** *Journal of Open Source
   Software* 2016; 1(5):27.
11. Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. **CheckM2: a rapid, scalable and accurate tool for
    assessing microbial genome quality using machine learning.** *Nature Methods* 2023; 20(8):1203–1212.
    doi:10.1038/s41592-023-01940-w.
12. Orakov A, Fullam A, Coelho LP, *et al.* **GUNC: detection of chimerism and contamination in
    prokaryotic genomes.** *Genome Biology* 2021; 22:178.
    doi:10.1186/s13059-021-02393-0.
13. Minkin I, Medvedev P. **Scalable multiple whole-genome alignment and locally collinear block
    construction with SibeliaZ.** *Nature Communications* 2020; 11:6327.
    doi:10.1038/s41467-020-19777-8.
14. Ondov BD, Starrett GJ, Sappington A, *et al.* **Mash Screen: high-throughput sequence containment
   estimation for genome discovery.** *Genome Biology* 2019; 20:232.
   doi:10.1186/s13059-019-1841-x.
15. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. **CheckM: assessing the quality of
    microbial genomes recovered from isolates, single cells, and metagenomes.** *Genome Research*
    2015; 25(7):1043–1055. doi:10.1101/gr.186072.114.
16. Sun Z, Huang S, Zhang M, *et al.* **Challenges in benchmarking metagenomic profilers.** *Nature
    Methods* 2021; 18(6):618–626. doi:10.1038/s41592-021-01141-3.
17. Meyer F, Fritz A, Deng ZL, *et al.* **Critical Assessment of Metagenome Interpretation: the second
    round of challenges.** *Nature Methods* 2022; 19(4):429–440. doi:10.1038/s41592-022-01431-4.
18. Huang W, Li L, Myers JR, Marth GT. **ART: a next-generation sequencing read simulator.**
    *Bioinformatics* 2012; 28(4):593–594.
19. Bowers RM, Kyrpides NC, Stepanauskas R, *et al.* **Minimum information about a single amplified genome
    (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea.** *Nature Biotechnology*
    2017; 35(8):725–731. doi:10.1038/nbt.3893.
20. Seemann T. **Barrnap: basic rapid ribosomal RNA predictor.** Software, https://github.com/tseemann/barrnap.
21. Eddy SR. **Accelerated profile HMM searches.** *PLoS Computational Biology* 2011; 7(10):e1002195.
    doi:10.1371/journal.pcbi.1002195.
22. Shaw J, Yu YW. **Fast and robust metagenomic sequence comparison through sparse chaining with skani.**
    *Nature Methods* 2023; 20(11):1661–1665. doi:10.1038/s41592-023-02018-3.
23. Langmead B, Salzberg SL. **Fast gapped-read alignment with Bowtie 2.** *Nature Methods* 2012;
    9(4):357–359. doi:10.1038/nmeth.1923.
24. Li H. **Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM.** arXiv:1303.3997
    (2013).
25. Ye SH, Siddle KJ, Park DJ, Sabeti PC. **Benchmarking metagenomics tools for taxonomic
    classification.** *Cell* 2019; 178(4):779–794. doi:10.1016/j.cell.2019.07.010.
26. Anderson MJ. **A new method for non-parametric multivariate analysis of variance.** *Austral Ecology*
    2001; 26(1):32–46. doi:10.1046/j.1442-9993.2001.01070.x.
27. Baker DN, Langmead B. **Dashing: fast and accurate genomic distances with HyperLogLog.** *Genome
    Biology* 2019; 20:265. doi:10.1186/s13059-019-1875-0.
28. Marçais G, Kingsford C. **A fast, lock-free approach for efficient parallel counting of occurrences of
    k-mers.** *Bioinformatics* 2011; 27(6):764–770. doi:10.1093/bioinformatics/btr011.
29. Bankevich A, Nurk S, Antipov D, *et al.* **SPAdes: a new genome assembly algorithm and its
    applications to single-cell sequencing.** *Journal of Computational Biology* 2012; 19(5):455–477.
    doi:10.1089/cmb.2012.0021.
