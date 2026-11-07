## Results

### Overview of Strain2bScan and its two data modes (Fig 1)

Strain2bScan converts single-copy type-IIB 2bRAD tags into a strain-cluster-resolved marker database and supports both native BcgI 2bRAD-M reads and in-silico digested shotgun reads (Fig 1). The same canonical marker counts feed species gating, cluster detection, depth–breadth consistency filtering and abundance estimation. This overview motivates the two evaluation streams below: native 2bRAD-M data first, followed by shotgun data.

### Part I. Native 2bRAD data

#### 2bRAD tags retain strain-level signal that 16S does not (Fig 2)

Across 15 species with complete or near-complete genomes, between-strain 2bRAD-tag distance tracked whole-genome distance (median Spearman ρ = 0.94), whereas 16S distance did not (median ρ = 0.36; Fig 2A). Per-species rank–rank scatter plots showed the same separation (Supplementary Fig S7): 2bRAD preserved within-species strain ordering, whereas 16S collapsed divergent strains into discrete 16S-distance bands. Thus, the reduced 2bRAD marker space retained information that is lost by 16S.

#### DNA mock communities (Fig 3, Supplementary Fig S1)

Native BcgI 2bRAD-M profiling of four ATCC DNA mocks recovered strain-cluster composition in high-host and low-input settings. In MSA-1002, F1 was 0.976–1.0 at 90–99% human DNA and Bray–Curtis similarity was 0.745–0.756 (Fig 3). A 0.1-ng input retained F1 = 0.952, whereas 0.01 ng approached the marker-depth noise floor and 0.001 ng was below the detection floor. At 90–99% host DNA, AUPR = 1.0. Supplementary Fig S1 shows the separate 62-species species-level titration and is not directly comparable with the primary cluster-resolved benchmark. Clean MSA-1003/1005/1007 profiles tested composition rather than the host ladder. The species-to-strain mapping for all mock communities is provided in Supplementary Table S14.

#### Real saliva contains individual-specific strain signatures (Fig 4)

Thirty-two native BcgI 2bRAD-M saliva libraries from eight subjects and four within-day timepoints were profiled against a 19-species oral panel. Species- and strain-level Bray–Curtis profiles gave similar community-level subject associations (PERMANOVA R² = 0.755 and 0.757), but strain features identified the subject in every leave-one-timepoint-out test (100% versus 78.1% for species; Fig 4). Within-subject Bray–Curtis distance was 0.327 versus 0.696 between subjects. At the species level, 15 species had at least ten detectable samples; strain profiles explained 43–85% of subject-level variation, with the largest R² for *Rothia mucilaginosa* (0.852), indicating that much of the individual signal was species- and strain-resolved. Pairwise sharing gave the same conclusion: dominant-cluster sharing was 0.538 within subjects across timepoints versus 0.217 between subjects (difference +0.320; pair-level bootstrap 95% CI 0.264–0.375), and any-cluster sharing was 0.801 versus 0.568 (Supplementary Table S11).

### Part II. Shotgun data

#### Accuracy across the complete 15-species simulation pool (Fig 5, Supplementary Table S15)

Strain2bScan was evaluated on all 2,025 ART error-modelled shotgun simulations from 15 species, followed by in-silico 2bRAD digestion. The pool comprised 2/3/5 co-present 0.95-cluster strains per species, five depth levels (0.5–10× per strain) and five replicates; Supplementary Table S15 enumerates the complete design. Median precision, recall and F1 were 1.000. Mean precision, recall and F1 were 0.890, 0.846 and 0.862, and mean F1 increased from 0.786 at 0.5× to 0.888 at 10×. Mean F1 was lowest for *Klebsiella pneumoniae* (0.519) and *Salmonella enterica* (0.725), whereas 11 species exceeded 0.90 (Fig 5).

#### Enzyme ladder is a marker-density control (Fig 6)

A cumulative enzyme ladder on 225 matched simulations asked whether additional restriction sites improve strain resolution. Increasing enzymes from one to 14 increased median strain-specific markers from 816 to 12,298 while median cluster count remained 19. Median precision/recall/F1 improved from 0.667/0.333/0.500 to 1.000/1.000/1.000; four enzymes reached median precision 1.000 and F1 0.667 (Fig 6). This experiment identifies marker density as a design constraint; the native 2bRAD and primary shotgun comparisons below use the prespecified BcgI and all-enzyme settings, respectively.

#### Matched 225-sample multi-tool comparison (Fig 7, Tables 1–3)

On 225 matched single-species shotgun simulations generated in all-enzyme marker space, Strain2bScan, StrainScan, StrainGST and inStrain were compared in their own reference or call spaces. Strain2bScan and StrainScan both achieved species-median precision 1.0 across the 204 completed StrainScan pairs. Strain2bScan reached full median recall from 1× onward, whereas StrainScan was non-monotonic. Among completed pairs, paired mean differences included zero (+0.063 recall; 95% CI −0.037 to 0.133; +0.005 F1; 95% CI −0.091 to 0.069). StrainGST had precision and recall 1.0 at all depths in its 0.90-reference space, including 0.5×, because its strain-search reference set retained more separate low-divergence references than the 0.95-cluster panels. inStrain showed high recall but reduced precision at ≥1× because mappings to the full genome pool were assigned to 0.95 truth clusters (Fig 7A–C). In a failure-as-absence sensitivity for StrainScan, Strain2bScan had +0.127 higher recall (95% CI 0.007 to 0.248). Multi-species F1 was lower for Strain2bScan than StrainScan at all three depths, so the scale advantage should not be interpreted as a uniform accuracy advantage.

#### Opt-in FracMinHash mode

An opt-in FracMinHash k-mer mode (k = 31; scale 100) was evaluated on all 225 matched simulations. Mean
precision, recall and F1 were 0.824, 0.750 and 0.776; medians were 1.000, 1.000 and 0.889. Mean F1 was
0.626 at 0.5× and reached 0.830–0.835 at ≥3×. This mode retained full-k-mer sensitivity for many
high-diversity panels but was less robust than the 2bRAD primary mode across low-diversity and
near-clonal species.

#### Computational performance across 225 simulations and 15 species (Fig 8)

Computational performance was summarized from the 15-species simulation design rather than a single-species panel (Fig 8). Across the matched per-depth Strain2bScan/StrainScan runs, median wall time and peak RSS were 0.56 s and 0.059 GB for Strain2bScan versus 4.79 s and 0.811 GB for StrainScan. StrainGST required 18.66 s and 1.241 GB per single-species sample because it k-merized the sample and searched species pan-genome databases. inStrain BAM profiling required a median 1.80 s but 18.8 GB peak RSS; this timing excludes Bowtie2 alignment and is therefore not an end-to-end runtime. Database construction showed the same pattern across the 15 species panels, and multi-species profiling remained approximately flat in species count for Strain2bScan, whereas per-species comparator workflows scaled with the number of species searched.

#### Shotgun DNA mocks separate tools under host contamination (Fig 9, Table 9)

Four ATCC shotgun mocks were profiled with Strain2bScan, StrainScan and inStrain in their own call units. In MSA-1002 at 99% host, Strain2bScan retained F1 1.0 and Bray–Curtis dissimilarity 0.302, whereas StrainScan and inStrain had F1 0.095 and 0.333 (Fig 9). These are descriptive tool-specific contrasts, not pooled estimates. Clean-mock thresholded precision varied with input composition.

#### Real shotgun saliva shows dense calls but not denser strain-level host discrimination (Fig 10)

The complete 32-pair PRJNA1131785 WMS/native comparison used the same oral panel and gates. Full-depth WMS produced 2,487 strain-cluster calls (median 79 per sample) versus 1,368 native calls (median 35.5); 495 calls were shared (Fig 10A). Thus, 19.9% of WMS calls were present in native profiles and 36.2% of native calls were recovered in WMS. Unique calls were abundance-poor: WMS-only calls represented 80.1% of WMS calls but only 14.8% of WMS-detected sample fraction, whereas native-only calls represented 63.8% of native calls but 26.6% of native-detected sample fraction (Fig 10B).

Full-depth WMS also reproduced the saliva time-course analysis. WMS species and strain profiles gave subject PERMANOVA R² of 0.628 and 0.572, respectively, and leave-one-timepoint-out subject accuracy of 78.1% and 84.4%. Thus, denser WMS calling did not produce stronger strain-level host discrimination than native 2bRAD. Because sequencing depth, representation and host content differ, this is a modality-scope comparison rather than accuracy validation.

#### Real shotgun cohorts (Tables 4–7)

Public WGS cohorts showed three application boundaries (Supplementary Table S3). A pregnancy multi-site subset produced 34 strain-cluster calls across five of six libraries, with saliva cluster-rich two-time-point descriptive profiles. Preterm vaginal/meconium libraries produced no strain calls because strain-marker evidence was sparse and panels were incomplete. Preterm-infant stool profiling detected weekly turnover with a generic panel and persistent *Bifidobacterium bifidum* LHCA82 with a cohort-specific panel. Six isolate WGS read sets recovered their expected panel units (6/6), but this was a panel/read compatibility result rather than independent classification.

#### Real rhesus macaque GI cohort combines 2bRAD and shotgun evidence (Fig 11, Supplementary Tables S5–S13)

A 148-library rhesus macaque GI cohort included 128 spatial samples (66 luminal and 62 mucosal) and 20 body controls. A cohort-specific panel built from 315 study MAGs represented 18 gut species and 303 0.95 strain-cluster units; after cross-species marker restriction, 29,162 of 30,008 BcgI markers remained available. Strain calls were recovered in 113/128 GI samples (88.3%): 56/66 luminal and 57/62 mucosal. Only two body controls produced calls, with three total calls (Supplementary Table S5).

Strain calls partitioned along the GI axis. Eighty-nine clusters met the LI-enrichment screen; in 50 matched mucosal/luminal positions, 13 clusters remained luminal-enriched after Bonferroni correction (Fig 11C,D; Supplementary Tables S7–S8). Among 51 luminal samples with calls in both modalities, species profiles were more concordant than exact cluster profiles (median Jaccard 0.889 versus 0.605), and shotgun yielded more calls per sample (median 123 versus 84), while only 234/580 dominant-cluster comparisons matched exactly (Supplementary Table S9).

In matched LI samples, strain composition was coupled to bile-acid composition and metabolic pathways after residualizing monkey identity and GI position (Supplementary Table S10). Among 219 clusters and four bile-acid classes, 410/876 associations passed BH FDR < 0.05; 409 were negative, consistent with an LI luminal strain reservoir tracking a lower conjugated/high-mass bile-acid state. Functional coupling was cluster-specific for butyrate, propionate and N-acetyl-sugar pathways. Pairwise sharing showed spatial structure: within-monkey between-position dominant-cluster sharing was 0.379 versus 0.159 between monkeys at same-age same-region/compartment sites (difference +0.220; 95% CI 0.176–0.262). Between old and young monkeys, LI luminal any-cluster sharing was 0.870 but LI mucosal dominant sharing was 0.071 (Fig 11G,H; Supplementary Tables S12–S13). These patterns support spatial retention of strain reservoirs but do not establish SNV/popANI-level transmission.

### Comparator scope

The executed strain-resolved comparators were StrainScan, StrainGST and inStrain. StrainPhlAn was not included in precision/recall/F1 benchmarking because it infers strain phylogeny and types rather than cluster abundance calls. Accuracies are descriptive within prespecified, tool-specific resolution spaces and are not pooled.

**Conclusion.** Together, these benchmarks show that 2bRAD markers retain strain-resolved signal, that native 2bRAD-M remained robust in high-host and low-input DNA mocks and saliva, and that Strain2bScan extends reduced-representation sequencing to scalable shotgun strain-cluster profiling with lower computational cost.
