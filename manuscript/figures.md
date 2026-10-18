# Strain2bScan; Figure legends and file map

Figure numbering follows the two-input structure. Each entry gives the figure file(s) in
`figures/`, the source table(s) in `results/`, and the per-experiment doc in `docs/`.
Multi-panel figures list the panel→file mapping.

**Numbered figure set:** `figures/numbered/Fig01…Fig10` and `FigS1…FigS6` (PNG + PDF) are the
manuscript-ordered figures. Regenerate with `python scripts/number_figures.py` and refresh the
manifest with `python scripts/make_figure_manifest.py`.

## Main figures

**Figure 1. Strain2bScan overview and quality control.** `figures/source_panels/overview.{png,pdf}`
**(A)** Reference construction: type-IIB (2bRAD) digestion into single-copy 25–33 bp tags, within-species clustering at 0.95 Jaccard, marker classification and a cluster × marker database. **(B)** Profiling: one sample digest, species-level marker gating, cross-species marker restriction, cluster detection and abundance estimation for conventional shotgun or native BcgI input. **(C)** Depth-breadth consistency separates a genuine low-depth strain from a partial-marker shadow.

**Figure 2. 2bRAD tags track genome-wide strain distance; 16S does not.** `figures/source_panels/mash_2brad_vs_16s.{png,pdf}`
**(A)** Per-species Spearman correlation between whole-genome strain distance and 2bRAD-tag or 16S distance across 15 complete/near-complete genome panels; 2bRAD median 0.94 versus 16S median 0.36. **(B)** Rank-rank matrices show that 2bRAD preserves genome-wide strain ordering whereas 16S collapses many unrelated pairs.

# Part I — native BcgI 2bRAD-M performance

**Figure 3. The 2bRAD enzyme set is a marker-density and recall control.** `figures/source_panels/enzyme_sweep.{png,pdf}`
In simulated samples from 14 resolvable species of the 15-species pool, increasing the enzyme set from 1 to 14 leaves precision at 1.0 while increasing median recall from 0.50 to 1.00 and strain-specific marker yield from 976 to 12,647. Approximately four enzymes recover most recall; BcgI alone supports native 2bRAD-M libraries. The sweep is an in-silico design experiment, not a native wet-lab library benchmark.

**Figure 4. Native BcgI 2bRAD-M across four ATCC DNA mocks.** `figures/source_panels/fig06_native_2brad_mocks.{png,pdf}`
MSA-1002, MSA-1003, MSA-1005 and MSA-1007 were profiled against a 28-species, 164-genome containment tree. At 90, 95 and 99% human DNA, F1 was 0.976–1.0 and Bray–Curtis similarity 0.745–0.756. A 0.1-ng input retained F1 0.952; 0.01 ng approached the marker-depth noise floor and 0.001 ng was below the detection floor.

**Figure 5. Real saliva contains individual-specific, within-day strain signatures.** `figures/source_panels/saliva_fig07.{png,pdf}`
Thirty-two native BcgI 2bRAD-M saliva libraries from eight subjects sampled four times in one day were profiled against a 19-species oral panel. Leave-one-timepoint-out subject classification reached 100% with strain features versus 78.1% with species features. Within-subject Bray–Curtis distance was 0.327 versus 0.696 between subjects. This is a small within-subject result, not independent strain-identity validation.

**Figure 6. Paired saliva comparison shows directional recovery of low-abundance strains.** `figures/source_panels/saliva_fig08.{png,pdf}`
Shotgun R1 prefix libraries from usable paired samples were digested in silico and compared with matched native BcgI libraries. All 65 shotgun strain-cluster calls were present in native calls; native 2bRAD-M yielded 128–163 additional candidate calls per sample. Additional calls had lower median relative abundance and require independent validation.

# Part II — conventional shotgun metagenome performance

**Figure 7. Accuracy across the complete 2,025-sample single-species simulation pool.** `figures/source_panels/sim2025_summary.{png,pdf}`
Strain2bScan was evaluated on all 2,025 ART error-modelled shotgun simulations from 15 species, followed by in-silico 2bRAD digestion. Median precision, recall and F1 were 1.000; mean precision, recall and F1 were 0.890, 0.846 and 0.862. Mean F1 increased from 0.786 at 0.5× to 0.888 at 10×. Panels show per-species means, depth dependence and different-cluster versus same-cluster design.

**Figure 8. Matched 225-sample simulation comparison.** `figures/source_panels/sim_headtohead.{png,pdf}`
Strain2bScan, StrainScan, StrainGST and inStrain were compared on 225 matched single-species simulations and 12 multi-species communities. Accuracy panels use each tool's own reference/cluster space: 0.95 for Strain2bScan/StrainScan, 0.90 for StrainGST, and full species genome pools for inStrain, so accuracy is not pooled. Panels also show database-build cost, same-environment profiling cost and multi-species profiling time.

**Figure 9. Computational scaling and community throughput.** `figures/source_panels/performance.{png,pdf} + scalability.{png,pdf} + community_throughput.{png,pdf}`
For *C. acnes* shotgun input, Strain2bScan was 8× faster and 11× lighter than StrainScan. Database build and profiling parallelised to 4.6× and 5.8× speed at 16 threads. Community profiling was approximately flat in species count; on a 55-species panel the projected cohort advantage over per-species StrainScan was 121–146×.

**Figure 10. Shotgun DNA mocks separate tools under host contamination.** `figures/source_panels/fig11_wms_toolcompare.{png,pdf}`
Four ATCC shotgun mocks were profiled with Strain2bScan, StrainScan and inStrain, each in its own 0.95 cluster space. In MSA-1002 at 99% host, Strain2bScan retained F1 1.0 and Bray–Curtis dissimilarity 0.302, whereas StrainScan and inStrain had F1 0.095 and 0.333. Clean-mock thresholded precision varied with input composition.

## Supplementary figures

**Figure S1. ATCC MSA-1002 DNA-input titration (native BcgI 2bRAD-M).** `figures/source_panels/mock_msa1002_titration.{png,pdf}`: precision 1.0 and full recall to 0.1 ng input.

**Figure S2. Layer-1 gate calibration on the 55-species panel.** `figures/source_panels/gate_calibration.{png,pdf}`; species precision/recall versus the marker floor at normal and low depth.

**Figure S3. Cost of unified-database expansion (20- versus 28-species combined tree).** `figures/source_panels/figS_tree_expansion.{png,pdf}`. MSA-1002 is unchanged by tree expansion; MSA-1003 has one lower-threshold replicate on the larger tree, confining cost to the approximately 1× marker-depth noise floor.

**Figure S4. inStrain requires a dereplicated reference (MSA-1002 shotgun).** `figures/source_panels/figS_instrain_derep.{png,pdf}`. The non-dereplicated 164-genome reference inflates false positives; 98%-ANI dereplication restores precision to 0.97, recall to 0.97 and F1 to 0.99.

**Figure S5. Strain2bScan versus StrainScan on three StrainScan-curated reference sets.** `figures/source_panels/species_expansion.{png,pdf}`. Both tools reached precision 1.0 for *A. muciniphila* and *P. copri*; Strain2bScan matched or exceeded recall and was faster and lighter. Strain2bScan completed near-clonal *M. tuberculosis* in 0.89 s, whereas StrainScan did not complete.

**Figure S6. Reference-completeness optimization with Jaccard and containment clustering.** `figures/source_panels/containment_mechanism.{png,pdf} + refqual_figure.{png,pdf}`
Fixed simulated reads were profiled while truth-strain references were degraded. Default Jaccard clustering lost precision and recall as completeness fell; opt-in maximum-containment clustering restored accuracy at 90–95% completeness and prevented near-clonal *M. tuberculosis* fragmentation. This is an algorithmic optimization analysis and is therefore supplementary to the primary performance benchmarks.
