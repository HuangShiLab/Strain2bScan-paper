## Results

### Overview of Strain2bScan and its two data modes (Fig 1)

Strain2bScan reduces reference genomes to canonical, single-copy 2bRAD tags, clusters genomes within each species at 0.95 Jaccard similarity, and builds a compact cluster-by-marker database. Profiling digests a sample once, gates species on panel-specific markers, and scores strain-resolved clusters on unique markers. The same tag space accepts native BcgI 2bRAD-M libraries and in-silico-digested shotgun reads. Tag definitions match Fast2bRAD-M, so species and strain profiling share one marker space.

### 2bRAD tags retain strain-level signal that 16S does not (Fig 2)

Across 15 species with complete or near-complete genomes, 2bRAD-tag distances tracked whole-genome strain distances in every species (median Spearman ρ = 0.94; range 0.59–0.99), whereas 16S distances were much weaker (median ρ = 0.36). Thus, 16S resolved species, while 2bRAD markers retained genome-wide strain relationships.

**Part I. Native 2bRAD-M performance.**

#### In-silico marker-density simulation (Fig 3)

We reran the enzyme ladder on 225 matched shotgun-derived simulations, giving 1,125 sample-level predictions across 15 species. Increasing enzymes from 1 to 14 increased median strain-specific markers from 816 to 12,298 while median cluster count remained 19, showing that the ladder increased marker density rather than cluster resolution. Across species, median precision/recall/F1 improved from 0.667/0.333/0.500 with BcgI alone to 1.000/1.000/1.000 with 14 enzymes; four enzymes already reached median precision 1.000 and F1 0.667 (Fig 3; `results/enzyme_sweep_225/`). This in-silico experiment identifies marker density as a design constraint; native BcgI operation is tested below with real DNA mocks and saliva libraries.

#### DNA mocks (Fig 4)

Four ATCC whole-cell mocks were profiled against a unified 28-species, 164-genome tree. In the primary MSA-1002 comparison, native BcgI 2bRAD-M retained F1 = 0.952 at 0.1 ng DNA and F1 = 1.0 with AUPR = 1.0 at 90–99% human DNA. At 90, 95 and 99% host DNA, MSA-1002 F1 was 0.976–1.0 and Bray–Curtis similarity was 0.745–0.756. F1 was 1.0 at 1 ng, whereas 0.001 ng was below the practical detection floor. These MSA-1002 mocks support robust strain detection under high host content and low input, with a remaining noise floor near approximately 1× marker depth; clean MSA-1003/1005/1007 profiles tested composition rather than the host ladder.

#### Real saliva time-course (Fig 5)

Thirty-two native BcgI 2bRAD-M saliva libraries from eight subjects and four within-day timepoints were profiled against a 19-species oral panel. Strain features identified the subject in every leave-one-timepoint-out test (100%; species features 78.1%, 95% CI 0.563–1.000). The subject-aware strain-minus-species accuracy difference was 0.219 (95% CI 0.000–0.438), and within-subject Bray–Curtis distance was lower than between-subject distance (0.327 versus 0.696; exploratory PERMANOVA p = 8.05 × 10⁻¹⁹). These data show a stable, person-specific strain signature, although the small within-subject classifier does not establish a significant improvement over species features or independent strain identity.

#### Paired saliva: 2bRAD-M versus shotgun (Fig 6)

In an exploratory R1-prefix comparison (three of four downloaded libraries usable; WMS capped at 4 × 10⁸ reads and profiled with stricter species gates), all 65 shotgun strain-cluster calls were also present in native BcgI 2bRAD-M (65/65). Native libraries yielded 128–163 additional candidate strain-cluster calls per sample. This small, non-depth-matched observation is directional only and does not establish a systematic native advantage.

A uniform gate-matched rebuild (Supplementary Table S4) then profiled all 32 paired WMS libraries and all 32 native 2bRAD libraries with the same oral panel, species gate, support and coverage thresholds. Full-depth WMS produced 2,487 strain-cluster calls (median 79 per sample) versus 1,368 native calls (median 35.5); 495 calls were shared. Thus, 19.9% of WMS calls were present in native profiles and 36.2% of native calls were recovered in WMS. Because sequencing depth, representation and host content differ between modalities, these data show substantial modality-specific sensitivity rather than accuracy validation.


**Part II. Conventional shotgun performance.**

#### Complete 15-species simulation pool (Fig 7)

We evaluated current-source Strain2bScan on all 2,025 single-species ART error-modelled shotgun metagenomes generated from the same 15-species genome pool. Each sample contained 2, 3 or 5 co-present strains at 0.5–10× per-strain depth. Median precision, recall and F1 were 1.000 across the full pool; mean precision, recall and F1 were 0.890, 0.846 and 0.862. Mean F1 increased from 0.786 at 0.5× to 0.888 at 10×, while median precision remained 1.000. Performance was species-dependent: mean F1 was lowest for *Klebsiella pneumoniae* (0.519) and *Salmonella enterica* (0.725), whereas 11 species exceeded 0.90. These simulations establish the accuracy envelope of the 2bRAD marker framework without implying uniform species-level accuracy.

#### Matched 225-sample multi-tool benchmark (Fig 8, Tables 1–3)

On 225 matched single-species simulations, Strain2bScan and StrainScan both achieved species-median precision 1.0 across the 204 completed StrainScan pairs. Species-median recall and F1 were 1.000 and 1.000 for Strain2bScan versus 0.667 and 0.800 for StrainScan; paired recall and F1 differences included zero (+0.063 recall; 95% CI −0.037 to 0.133; +0.005 F1; 95% CI −0.091 to 0.069). StrainScan’s depth profile was non-monotonic, whereas Strain2bScan reached full median recall from 1× onward. In a failure-as-absence sensitivity (Supplementary Table S2) that retained all 225 intended entries, Strain2bScan had +0.127 higher recall (95% CI 0.007 to 0.248), because StrainScan lacked calls for 21 runs. In the archived database benchmark, which compared native Strain2bScan with emulated linux/amd64 StrainScan builds, construction used 249–614× less wall time and 43–138× less memory. In the stricter same-container profiling subset, per-sample profiling was 4–33× faster and 5–39× lighter. Multi-species F1 was lower for Strain2bScan than StrainScan at all three depths, so the scale advantage should not be interpreted as uniform accuracy advantage.

StrainGST and inStrain extended this comparison descriptively. StrainGST completed all simulations and ATCC shotgun mocks but used 53.31× more wall time and 25.49× more peak memory for single-species samples; accuracy was reported separately in its 0.90-reference space. inStrain completed all 225 matched samples and retained high recall, but median precision declined from 0.750 at 0.5× to 0.185 at 10× after mapping genome calls to truth 0.95 clusters. Because StrainGST uses 0.90-reference units, inStrain uses dereplicated genome references before truth-cluster mapping, and Strain2bScan/StrainScan use 0.95 clusters, accuracies are shown side by side but were not pooled.

#### Computational scaling (Fig 9)

For shotgun input in the same-container benchmark, Strain2bScan profiled the *C. acnes* panel in 0.86 s and 78 MB per sample versus 7.06 s and 828 MB for StrainScan, an approximately 8× runtime and 11× memory advantage. Database build and profiling parallelised to 4.6× and 5.8× speed at 16 threads. On a 55-species panel, runtime was approximately flat over the tested species range, giving a projected 121–146× cohort-level advantage because StrainScan lacks native multi-species profiling.

#### Supplementary curated-reference comparison (Supplementary Fig S5)

On StrainScan-curated reference sets, both tools reached precision 1.0 for *A. muciniphila* and *P. copri*. Strain2bScan matched or exceeded recall, was 17–23× faster and 15–24× lighter, and completed near-clonal *M. tuberculosis* in 0.89 s, whereas StrainScan did not complete. Low *M. tuberculosis* recall reflected collapse of the panel to five 0.95-similarity clusters.

#### DNA mocks (Fig 10, Table 9)

In the MSA-1002 shotgun host ladder, Strain2bScan retained F1 = 1.0 at 90, 95 and 99% host DNA. At 99% host, StrainScan and inStrain had F1 = 0.095 and 0.333 in their own call spaces. Clean-mock thresholded precision varied with input composition, so this was a descriptive comparator contrast rather than a pooled like-for-like accuracy estimate.

#### Real shotgun cohorts (Tables 4–7)

Public WGS cohorts showed three application boundaries (Supplementary Table S3). A pregnancy multi-site subset produced 34 strain-cluster calls across five of six libraries, with saliva cluster-rich two-time-point descriptive profiles. Preterm vaginal/meconium libraries produced no strain calls because strain-marker evidence was sparse and panels were incomplete. Preterm-infant stool profiling detected weekly turnover with a generic panel and persistent *Bifidobacterium bifidum* LHCA82 with a cohort-specific panel. Six isolate WGS read sets recovered their expected panel units (6/6), but this was a panel/read compatibility result rather than independent classification.

#### Real rhesus macaque GI cohort: strain-resolved spatial ecology in 148 samples (Fig 11, Supplementary Tables S5–S10)

We next asked whether the same framework could resolve strain structure in a complete real spatial cohort rather than an exploratory subset. We profiled 148 rhesus macaque libraries: 128 GI spatial samples (66 luminal and 62 mucosal) plus 20 body controls. A cohort-specific panel built from 315 study MAGs represented 18 gut species and 303 0.95 strain-cluster units; after cross-species marker restriction, 29,162 of 30,008 BcgI markers remained available for quantification. With the prespecified species and cluster gates, 113 of 128 GI samples (88.3%) yielded strain calls: 56/66 luminal and 57/62 mucosal samples. Only two body controls produced calls, with three total strain calls, indicating that the main GI signals were not explained by body-site contamination (Supplementary Table S5).

Strain calls partitioned along the GI axis. Eighty-nine clusters met the prespecified LI-enrichment screen; ten high-confidence examples included *Treponema_D_succinifaciens* C25, *WG_1_sp900539665* C22, *Treponema_D_succinifaciens* C2, *Succinivibrio_sp000431835* C6 and *Cryptobacteroides_sp900546445* C1 (Supplementary Table S7). In 50 matched mucosal/luminal positions, zero-filled within-species compositions identified 20 nominally significant clusters and 13 that remained significant after Bonferroni correction; all 13 were luminal-enriched (Supplementary Table S8). The strongest examples included *T._succinifaciens* C25, *UBA4372_sp900543815* C13/C5, *Cryptobacteroides_sp900546445* C13/C3 and *WG_1* C20. Thus, closely related strains of the same species were not uniformly distributed across mucosal and luminal niches.

The 66 luminal GI libraries also had shotgun and HUMAnN profiles, enabling a modality and function overlay. Among the 51 samples with calls in both modalities, species profiles were more concordant than exact cluster profiles (median Jaccard 0.889 versus 0.605), and shotgun yielded more calls per sample (median 123 versus 84). Only 234 of 580 dominant-cluster comparisons matched exactly (40.3%). These results support native 2bRAD-M as a reduced-depth strain-aware measurement while indicating that shotgun provides denser cluster resolution (Supplementary Table S9).

In matched LI samples, strain composition was coupled to bile-acid composition and selected metabolic pathways after residualizing monkey identity and GI position (Supplementary Table S10). Among 219 clusters and four precursor-mass bile-acid classes, 410 of 876 associations passed BH FDR < 0.05; 409 were negative, showing that the LI strain reservoir mainly tracked a lower conjugated/high-mass bile-acid state. For example, *Cryptobacteroides_sp000433355* C1 was negatively associated with the conjugated/high-mass class (ρ = −0.943, q = 1.0 × 10⁻⁸), whereas *UBA636_sp900546285* C11 was positively associated with the dihydroxy class (ρ = 0.648, q = 2.5 × 10⁻³). Functional coupling was also cluster-specific: *UBA4372_sp900543815* C4 tracked pyruvate- and acetyl-CoA-to-butyrate pathways, *Cryptobacteroides_sp000433355* C5 tracked pyruvate-to-propionate, and *Treponema_D_succinifaciens* C2 and *Succinivibrio_sp000431835* C4 tracked N-acetylneuraminate and N-acetylglucosamine-related degradation pathways. These associations describe ecological coupling in real samples; they do not assign individual pathway genes to a strain cluster or establish causality.

#### Strain sharing across timepoints and GI sites (Fig 12, Supplementary Tables S11–S13)

We next summarized the same 0.95 cluster calls as pairwise sharing rates. The denominator was the number
of species with strain calls in both samples. The primary metric counted a species when its most abundant
cluster was identical; a secondary coexistence-aware metric counted at least one shared cluster label.
In the 32 saliva libraries, dominant-cluster sharing was 0.538 within subjects across timepoints versus
0.217 between subjects (mean difference +0.320; pair-level bootstrap 95% CI 0.264–0.375), and any-cluster
sharing was 0.801 versus 0.568 (difference +0.233; 95% CI 0.191–0.271). Thus, within-day saliva samples
retained host-specific strain signatures.

In the 128 GI libraries, within-monkey between-position sharing exceeded between-monkey same-age
same-region/compartment sharing (dominant-cluster mean 0.379 versus 0.159; difference +0.220; 95% CI
0.176–0.262). Between old and young monkeys at corresponding sites, large-intestinal luminal sharing was
high for any-cluster labels (0.870) but much lower for dominant clusters (0.222), whereas LI mucosal
dominant sharing was 0.071. These patterns support spatial retention of gut strain reservoirs but do not
by themselves establish transmission, because exact label sharing is at 0.95 marker-cluster resolution
rather than SNV/popANI identity.

### Comparator scope

The primary executed strain-resolved comparators were StrainScan and inStrain; StrainGST was run separately because it uses StrainGE 0.90-reference databases. StrainPhlAn was not included in this P/R/F1 benchmark because it infers strain phylogeny and types rather than cluster abundance calls. Accuracies are therefore descriptive within prespecified, tool-specific resolution spaces and are not pooled.

**Conclusion.** Together, these benchmarks show that 2bRAD markers retain strain-resolved signal, that native 2bRAD-M remained robust in the MSA-1002 high-host and low-input settings, and that Strain2bScan extends reduced-representation sequencing to scalable shotgun strain-cluster profiling with lower same-container profiling cost.
