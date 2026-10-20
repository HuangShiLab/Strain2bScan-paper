## Results

### Overview of Strain2bScan and its two data modes (Fig 1)

Strain2bScan reduces reference genomes to canonical, single-copy 2bRAD tags, clusters genomes within each species at 0.95 Jaccard similarity, and builds a compact cluster-by-marker database. Profiling digests a sample once, gates species on panel-specific markers, and scores strain-resolved clusters on unique markers. The same tag space accepts native BcgI 2bRAD-M libraries and in-silico-digested shotgun reads. Tag definitions match Fast2bRAD-M, so species and strain profiling share one marker space.

### 2bRAD tags retain strain-level signal that 16S does not (Fig 2)

Across 15 species with complete or near-complete genomes, 2bRAD-tag distances tracked whole-genome strain distances in every species (median Spearman ρ = 0.94; range 0.59–0.99), whereas 16S distances were much weaker (median ρ = 0.36). Thus, 16S resolved species, while 2bRAD markers retained genome-wide strain relationships.

**Part I. Native 2bRAD-M performance.**

#### In-silico marker-density simulation (Fig 3)

We reran the enzyme ladder on the 225 matched simulations, giving 1,125 sample-level predictions across 15 species. Increasing enzymes from 1 to 14 increased median strain-specific markers from 816 to 12,298 while median cluster count remained 19, showing that the ladder increased marker density rather than cluster resolution. Across species, median precision/recall/F1 improved from 0.667/0.333/0.500 with BcgI alone to 1.000/1.000/1.000 with 14 enzymes; four enzymes already reached median precision 1.000 and F1 0.667 (Fig 3; `results/enzyme_sweep_225/`). BcgI therefore enables native 2bRAD-M operation, although single-enzyme panels are panel-dependent and less precise.

#### DNA mocks (Fig 4)

Four ATCC whole-cell mocks were profiled against a unified 28-species, 164-genome tree. In the primary MSA-1002 comparison, native BcgI 2bRAD-M retained F1 at 1e-4 = 0.909 at 0.1 ng DNA and F1 = 1.0 with AUPR = 1.0 at 90–99% human DNA. Across the four mocks, F1 was 0.976–1.0 at 90–99% host DNA and Bray–Curtis similarity was 0.745–0.756. F1 was 0.952 at 0.1 ng DNA and 1.0 at 1 ng, whereas 0.001 ng was below the practical detection floor. These mocks therefore support robust strain detection under high host content and low input, with a remaining noise floor near approximately 1× marker depth.

#### Real saliva time-course (Fig 5)

Thirty-two native BcgI 2bRAD-M saliva libraries from eight subjects and four within-day timepoints were profiled against a 19-species oral panel. Strain features identified the subject in every leave-one-timepoint-out test (100% versus 78.1% for species features), and within-subject Bray–Curtis distance was lower than between-subject distance (0.327 versus 0.696; p = 8.05 × 10⁻¹⁹). These data show a stable, person-specific strain signature, although the small within-subject classifier is not independent strain-identity validation.

#### Paired saliva: 2bRAD-M versus shotgun (Fig 6)

In three usable paired samples, all 65 shotgun strain-cluster calls were also present in native BcgI 2bRAD-M (65/65). Native libraries yielded 128–163 additional candidate strain-cluster calls per sample; these calls had lower median relative abundance than shared calls (0.0029 versus 0.0097; Mann–Whitney p = 1.2 × 10⁻²³). This directional comparison supports greater recovery of low-abundance strains by native 2bRAD-M, but the additional calls require independent validation.

**Part II. Conventional shotgun performance.**

#### Complete 15-species simulation pool (Fig 7)

We evaluated current-source Strain2bScan on all 2,025 single-species ART error-modelled shotgun metagenomes generated from the same 15-species genome pool. Each sample contained 2, 3 or 5 co-present strains at 0.5–10× per-strain depth. Median precision, recall and F1 were 1.000 across the full pool; mean precision, recall and F1 were 0.890, 0.846 and 0.862. Mean F1 increased from 0.786 at 0.5× to 0.888 at 10×, while median precision remained 1.000. These simulations establish the accuracy envelope of the 2bRAD marker framework.

#### Matched 225-sample multi-tool benchmark (Fig 8, Tables 1–3)

On 225 matched single-species simulations, Strain2bScan and StrainScan both achieved species-median precision 1.0 across the 204 completed StrainScan pairs. Species-median recall and F1 were 1.000 and 1.000 for Strain2bScan versus 0.667 and 0.800 for StrainScan; paired recall and F1 differences included zero (+0.063 recall; 95% CI −0.037 to 0.133; +0.005 F1; 95% CI −0.091 to 0.069). StrainScan’s depth profile was non-monotonic, whereas Strain2bScan reached full median recall from 1× onward. In a failure-as-absence sensitivity that retained all 225 intended entries, Strain2bScan had +0.127 higher recall (95% CI 0.007 to 0.248), because StrainScan lacked calls for 21 runs. Strain2bScan was also faster and lighter: database construction used 249–614× less wall time and 43–138× less memory, and per-sample profiling was 4–33× faster and 5–39× lighter. Multi-species F1 was lower for Strain2bScan at all three depths, so the scale advantage should not be interpreted as uniform accuracy superiority.

StrainGST and inStrain extended this comparison. StrainGST completed all simulations and ATCC shotgun mocks but used 53.31× more wall time and 25.49× more peak memory for single-species samples; accuracy was reported separately in its 0.90-reference space. inStrain completed all 225 matched samples and retained high recall, but median precision declined from 0.750 at 0.5× to 0.185 at 10× after mapping genomes to truth 0.95 clusters. Because each tool uses a different reference or cluster space, accuracies were not pooled.

#### Computational scaling (Fig 9)

For shotgun input, Strain2bScan profiled the *C. acnes* panel in 0.86 s and 78 MB per sample versus 7.06 s and 828 MB for StrainScan, an approximately 8× runtime and 11× memory advantage. Database build and profiling parallelised to 4.6× and 5.8× speed at 16 threads. On a 55-species panel, runtime was approximately flat in species count, giving a projected 121–146× cohort-level advantage because StrainScan lacks native multi-species profiling.

#### Supplementary curated-reference comparison (Supplementary Fig S5)

On StrainScan-curated reference sets, both tools reached precision 1.0 for *A. muciniphila* and *P. copri*. Strain2bScan matched or exceeded recall, was 17–23× faster and 15–24× lighter, and completed near-clonal *M. tuberculosis* in 0.89 s, whereas StrainScan did not complete. Low *M. tuberculosis* recall reflected collapse of the panel to five 0.95-similarity clusters.

#### DNA mocks and real shotgun cohorts (Fig 10, Tables 4–7)

On four ATCC shotgun mocks, Strain2bScan retained F1 = 1.0 in MSA-1002 at 90, 95 and 99% host DNA. At 99% host, StrainScan and inStrain had F1 = 0.095 and 0.333. Clean-mock thresholded precision varied with input composition, so the primary advantage was retention under host contamination rather than uniform clean-mock superiority.

Public WGS cohorts showed three application boundaries. A pregnancy multi-site subset produced 34 strain-cluster calls across five of six libraries, with saliva cluster-rich and longitudinal profiles. Preterm vaginal/meconium libraries produced no strain calls because strain-marker evidence was sparse and panels were incomplete. Preterm-infant stool profiling detected weekly turnover with a generic panel and persistent *Bifidobacterium bifidum* LHCA82 with a cohort-specific panel. Six isolate WGS read sets recovered their expected panel units (6/6), but this was a panel/read compatibility result rather than independent classification.

### Comparator scope

The primary executed strain-resolved comparators were StrainScan and inStrain; StrainGST was run separately because it uses StrainGE 0.90-reference databases. StrainPhlAn was not forced into the P/R/F1 framework because it infers strain phylogeny/types rather than cluster abundance calls. Accuracies were therefore compared only within prespecified, tool-specific resolution spaces.

**Conclusion.** Together, these benchmarks show that 2bRAD markers retain strain-resolved signal, that native 2bRAD-M remains robust in high-host and low-input settings, and that Strain2bScan extends reduced-representation sequencing to scalable shotgun strain profiling with substantially lower database and profiling cost.
