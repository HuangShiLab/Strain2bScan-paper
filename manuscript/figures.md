## Main figures

**Figure 1. Strain2bScan overview and quality control.** `figures/source_panels/overview.{png,pdf}`
**(A)** Reference construction: type-IIB (2bRAD) digestion into single-copy 25–33 bp tags, within-species clustering at 0.95 Jaccard, marker classification and a cluster × marker database. **(B)** Profiling: one sample digest, species-level marker gating, cross-species marker restriction, cluster detection and abundance estimation for conventional shotgun or native BcgI input. **(C)** Depth–breadth consistency separates a genuine low-depth strain from a partial-marker shadow.

**Figure 2. 2bRAD tags track genome-wide strain distance; 16S does not.** `figures/source_panels/mash_2brad_vs_16s.{png,pdf}`
**(A)** Per-species Spearman correlation between whole-genome strain distance and 2bRAD-tag or 16S distance across 15 complete/near-complete genome panels; 2bRAD median 0.94 versus 16S median 0.36. **(B)** Rank-rank scatter plots show that 2bRAD preserves genome-wide strain ordering whereas 16S collapses many unrelated pairs.

**Figure 3. Native BcgI 2bRAD-M across four ATCC DNA mocks.** `figures/source_panels/fig06_native_2brad_mocks.{png,pdf}`
MSA-1002, MSA-1003, MSA-1005 and MSA-1007 were profiled against a 28-species, 164-genome containment tree. At 90, 95 and 99% human DNA, F1 was 0.976–1.0 and Bray–Curtis similarity 0.745–0.756. A 0.1-ng input retained F1 0.952; 0.01 ng approached the marker-depth noise floor and 0.001 ng was below the detection floor. Species-level colour keys correspond to the ATCC strain identifiers listed in Supplementary Table S14.

**Figure 4. Real saliva contains individual-specific strain signatures.** `figures/source_panels/saliva_fig07.{png,pdf}`
**(A,B)** Bray–Curtis PCoA of species-level and strain-level profiles coloured by subject. **(C)** Per-species strain-level subject R² for the 15 species with at least ten samples; the dashed line marks whole-community strain signal. **(D)** Within- versus between-subject Bray–Curtis distance. **(E)** Leave-one-timepoint-out subject classification. **(F,G)** Dominant-cluster and any-shared-cluster rates for pairwise saliva comparisons; denominators were species with calls in both samples. Thirty-two native BcgI 2bRAD-M libraries from eight subjects sampled four times in one day were profiled against a 19-species oral panel.

**Figure 5. Accuracy across the complete 2,025-sample single-species simulation pool.** `figures/source_panels/sim2025_summary.{png,pdf}`
**(A)** Per-species mean accuracy. **(B)** Accuracy across per-strain depth. **(C)** Accuracy by different-cluster versus same-cluster design. Strain2bScan was evaluated on all 2,025 ART error-modelled shotgun simulations from 15 species, followed by in-silico 2bRAD digestion. Median precision, recall and F1 were 1.000; mean precision, recall and F1 were 0.890, 0.846 and 0.862. Mean F1 increased from 0.786 at 0.5× to 0.888 at 10×.

**Figure 6. The 2bRAD enzyme set is a marker-density and recall control.** `figures/source_panels/enzyme_sweep.{png,pdf}`
**(A)** Species-median precision, recall and F1 across the cumulative enzyme ladder. **(B)** Median strain-specific marker yield and cluster count. On 225 matched simulations per enzyme level (1,125 sample-level predictions across 15 species), increasing enzymes from 1 to 14 raised median strain-specific markers from 816 to 12,298 while median cluster count stayed 19. Median precision/recall/F1 improved from 0.667/0.333/0.500 to 1.000/1.000/1.000; four enzymes reached median precision 1.000 and F1 0.667. This is an in-silico design experiment, not a native wet-lab library benchmark.

**Figure 7. Matched 225-sample simulation comparison.** `figures/source_panels/sim_headtohead.{png,pdf}`
**(A-C)** Single-species precision, recall and F1 across per-strain depth for Strain2bScan, StrainScan, StrainGST and inStrain. Strain2bScan/StrainScan use 0.95 clusters, StrainGST uses 0.90 references and inStrain uses full genome pools mapped to truth 0.95 clusters; accuracy is not pooled across these units. StrainGST remained accurate at 0.5× because its 0.90 reference space retained more separate low-divergence references. **(D)** Database-build time per species. **(E)** Same-container profiling cost. **(F)** Multi-species profiling time.

**Figure 8. Computational performance across 225 simulations and 15 species.** `figures/source_panels/performance.{png,pdf} + scalability.{png,pdf} + community_throughput.{png,pdf}`
The panels summarize median per-sample wall time and peak RSS across depth, database-build cost across the 15 species panels and multi-species throughput. Strain2bScan remained faster and lighter than StrainScan across matched per-depth runs. StrainGST included sample k-merization and species searches. inStrain points show BAM profiling only and exclude Bowtie2 alignment.

**Figure 9. Shotgun DNA mocks separate tools under host contamination.** `figures/source_panels/fig11_wms_toolcompare.{png,pdf}`
Four ATCC shotgun mocks were profiled with Strain2bScan, StrainScan and inStrain in their own call units. In MSA-1002 at 99% host, Strain2bScan retained F1 1.0 and Bray–Curtis dissimilarity 0.302, whereas StrainScan and inStrain had F1 0.095 and 0.333. These are descriptive tool-specific contrasts, not pooled estimates.

**Figure 10. Full-depth paired saliva WMS reanalysis.** `figures/source_panels/saliva_fig08.{png,pdf}`
**(A)** Shared, WMS-only and native-2bRAD-only strain-cluster calls for all 32 gate-matched WMS/native pairs. **(B)** Within-modality partition of detected sample fraction: WMS-only calls represented 80.1% of WMS calls but only 14.8% of WMS-detected sample fraction; native-only calls represented 63.8% of native calls but 26.6% of native-detected sample fraction. **(C,D)** Full-depth WMS time-course reanalysis shows subject-specific structure but lower strain-level community resolution than native 2bRAD.

**Figure 11. Real rhesus macaque GI cohort reveals LI-luminal strain reservoirs.** `figures/numbered/Fig11_monkey_spatial_cohort.{png,pdf}`
**(A)** Distribution of 128 GI libraries across four rhesus macaques and five GI regions. **(B)** Species prevalence in GI, luminal, mucosal and body-control libraries. **(C)** Ten representative LI-enriched clusters; the open red circle marks peak LI prevalence. **(D)** Thirteen matched mucosal/luminal strain clusters that remained significant after Bonferroni correction. **(E)** Paired luminal shotgun and native 2bRAD strain-call counts; colour gives strain-call Jaccard (median strain-call Jaccard 0.605; species Jaccard 0.889). **(F)** Representative Spearman associations between LI strain-cluster composition and bile-acid classes or HUMAnN pathways; BH q values are in Supplementary Table S10. **(G)** Mean dominant-cluster sharing for prespecified rhesus GI comparisons. **(H)** Old/young corresponding-site sharing by GI region and compartment. Exact label sharing describes 0.95 marker-cluster resolution, not SNV/popANI-validated strain identity or transmission.

## Supplementary figures

**Figure S1. ATCC MSA-1002 DNA-input titration uses native BcgI 2bRAD-M reads.** `figures/source_panels/mock_msa1002_titration.{png,pdf}` Native-input species precision and recall are shown for 0.001, 0.01, 0.1 and 1 ng DNA in the original 62-species species-level analysis; precision was 1.0, 0.01-ng recall was 0.95, and recall was complete at 0.1 ng and above. This titration is not directly comparable with the Figure 3 164-genome, 0.95-cluster primary benchmark.

**Figure S2. Layer-1 gate calibration uses the 55-species panel.** `figures/source_panels/gate_calibration.{png,pdf}` Species precision and recall are plotted against the minimum species-marker floor at normal and low depth.

**Figure S3. Cost of unified-database expansion (20- versus 28-species combined tree).** `figures/source_panels/figS_tree_expansion.{png,pdf}` MSA-1002 is unchanged by tree expansion; MSA-1003 has one lower-threshold replicate on the larger tree, confining cost to the approximately 1× marker-depth noise floor.

**Figure S4. inStrain requires a dereplicated reference (MSA-1002 shotgun).** `figures/source_panels/figS_instrain_derep.{png,pdf}` The non-dereplicated 164-genome reference inflates false positives; 98%-ANI dereplication restores precision to 0.97, recall to 0.97 and F1 to 0.99.

**Figure S5. Strain2bScan versus StrainScan on three StrainScan-curated reference sets.** `figures/source_panels/species_expansion.{png,pdf}` Both tools reached precision 1.0 for *A. muciniphila* and *P. copri*; Strain2bScan matched or exceeded recall and was faster and lighter. Strain2bScan completed near-clonal *M. tuberculosis* in 0.89 s, whereas StrainScan did not complete.

**Figure S6. Reference-completeness optimization with Jaccard and containment clustering.** `figures/source_panels/containment_mechanism.{png,pdf} + refqual_figure.{png,pdf}`
**(A)** Mechanism by which an incomplete genome becomes a Jaccard outlier and maximum containment merges it with its complete relative. **(B)** Median precision and recall across 14 resolvable species after reference degradation under Jaccard versus opt-in containment clustering, with near-clonal *M. tuberculosis* inset.

**Figure S7. Per-species rank–rank comparison of 2bRAD and 16S distances.** `figures/source_panels/rankrank_2brad_16s.{png,pdf}` For each of the 15 complete/near-complete genome panels, whole-genome distances are compared with 2bRAD and 16S distances. 2bRAD follows the diagonal; 16S forms discrete bands for many species.
