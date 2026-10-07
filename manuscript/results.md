## Results

### Overview of Strain2bScan and its two data modes (Fig 1)

Strain2bScan reduces reference genomes to canonical, single-copy 2bRAD tags, clusters genomes within each
species at 0.95 Jaccard similarity, and builds a compact cluster-by-marker database. Profiling digests a
sample once, gates species on panel-specific markers, and scores strain-resolved clusters on unique
markers. The same tag space accepts two entry points: in-silico digestion of conventional shotgun reads and
native BcgI 2bRAD-M libraries whose reads are already tags. The tag definitions match Fast2bRAD-M, so its
species layer and the Strain2bScan strain layer share one marker space.

### 2bRAD tags carry strain-level signal that 16S cannot (Fig 2)

We compared pairwise strain distances in whole-genome, 2bRAD-tag and 16S spaces across 15 species with
complete or near-complete genomes. Per species, 2bRAD distances tracked whole-genome distances in every
case (median Spearman ρ = 0.94; range 0.59–0.99), whereas 16S distances were much weaker (median ρ = 0.36),
with several confidence intervals overlapping zero (Fig 2). The rank-rank matrices show the mechanism:
2bRAD preserved the ordering of genome-wide strain pairs, whereas 16S collapsed many unrelated pairs to a
few conserved-gene distances. Thus, 16S resolved species but the 2bRAD marker set retained genome-wide
strain signal.

### Accurate and robust strain profiling (Fig 3–4)

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

## Part I; Native 2bRAD-M for low-biomass, high-host microbiomes

### A tunable enzyme set enables native BcgI 2bRAD operation (Fig 5)

Across 14 resolvable simulation-pool species, increasing the enzyme set from 1 to 14 left precision at 1.0
while increasing median recall monotonically from 0.50 to 1.00 and the median strain-specific marker yield
from 976 to 12,647 (Fig 5). Cluster count was invariant, so additional enzymes increased marker density and
recall rather than resolution. BcgI alone under-resolved near-clonal species, whereas approximately four
enzymes recovered most recall at much lower marker cost. Single-enzyme BcgI operation nevertheless enables
direct profiling of native BcgI 2bRAD-M libraries.

### Strain-level identification and quantification across four DNA mocks (Fig 6)

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

### Real saliva: individual-specific, temporally stable strain signatures (Fig 7)

We profiled native BcgI 2bRAD saliva from 8 subjects sampled at four times of day (32 libraries) against a
19-species oral reference panel. Subject PERMANOVA R² was similar for strain and species profiles (0.757
versus 0.755; both p = 2 × 10⁻⁴), but leave-one-timepoint-out nearest-neighbour classification identified
the host with 100% accuracy from strain features versus 78.1% from species features (Fig 7). The strongest
single-species association was *Neisseria subflava* (subject R² = 0.808); all 13 testable species were
significant at nominal p values. Within-subject Bray–Curtis distance was lower than between-subject distance
(0.327 versus 0.696; p = 8.05 × 10⁻¹⁹), indicating a stable, person-specific strain signature.

### Native 2bRAD detects candidate low-abundance strains missed by truncated shotgun (Fig 8)

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

## Part II; Conventional metagenomes at community scale

### Fast, light, and scalable to whole communities (Fig 9)

For shotgun input, Strain2bScan profiled the *C. acnes* panel in 0.86 s and 78 MB per sample, versus 7.06 s
and 828 MB for StrainScan, an approximately 8× runtime and 11× memory advantage (Fig 9A). Database build and
profiling parallelised to 4.6× and 5.8× speed at 16 threads (Fig 9B). Community profiling was approximately
flat in species count: on a 55-species panel, Strain2bScan profiled samples in 2.6–3.1 s. Because StrainScan
has no multi-species mode, its projected cost was the measured single-species runtime multiplied by species
and sample count; the projected advantage was 121–146× (Fig 9C).

### Matches or exceeds StrainScan on its own databases (Fig 10)

On StrainScan-curated reference sets, both tools reached precision 1.0 for *A. muciniphila* and *P. copri*
(Fig 10). Strain2bScan matched or exceeded recall (0.93 versus 0.24 and 0.94 versus 0.90), was 17–23×
faster and 15–24× lighter, and completed near-clonal *M. tuberculosis* in 0.89 s, whereas StrainScan did not
complete. Its low *M. tuberculosis* recall reflects the resolution limit of a panel that collapses to five
0.95-similarity clusters.

The shotgun mode was then stress-tested in ATCC mocks and compared with saliva. It preserved detection and
abundance under host contamination in the primary MSA-1002 comparison (Fig 12), and its saliva calls were
contained in the native-2bRAD call set (Fig 8). These are concordance and stress-test results, not
independent validations of strain identity.

### Systematic head-to-head on a 15-species simulated benchmark (Fig 11, Tables 1–3)

On a common 15-species simulation pool, both tools built databases from the same genomes and profiled the
same simulated reads. Each tool was scored in its own 0.95-similarity cluster space. The reproducible
comparator rerun comprised 225 matched single-species entries; 204 completed and were paired to the
corresponding Strain2bScan runs, whereas 21 ended without calls, including all *Salmonella enterica*
entries. Across the 204 completed runs from 14 species, species-median precision was 1.0 for both tools;
species-median recall and F1 were 0.733 and 0.844 for Strain2bScan versus 0.667 and 0.800 for StrainScan
(Fig 11, Table 1). In a species-cluster bootstrap, the paired mean differences in recall (−0.014; 95% CI
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

### Strain-level profiling on shotgun, and the advantage under host contamination (Fig 12)

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

### Exploratory public WGS cohorts: longitudinal observations, low-biomass specificity, and isolate-derived panels (Tables 4–6)

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
