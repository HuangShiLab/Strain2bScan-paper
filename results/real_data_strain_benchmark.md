# Strain2bScan strain-level vs species-level classification benchmark

## Datasets

### ECC (Early Childhood Caries)
- 38 saliva BcgI-2bRAD samples: 19 Healthy (H) vs 19 Caries (C)
- Species-level ground truth: `feature_table_species.txt` (970 GTDB species)
- Strain DB: 27 oral species from kraken16 genome pool
- Profile command: `strain2bscan multi-profile --dbs dbs_v2 --enzyme BcgI --min-species-markers 50 --min-species-detect 3`

### Lim_ORPI (clean vs unclean denture)
- 97 BcgI-2bRAD samples: 55 Clean_Denture vs 42 Unclean_Denture
- Species-level ground truth: `Abundance_Stat.all.xls` (1601 species)
- Strain DB: 30 species selected from top abundant/prevalent taxa in the species table
- Profile command: same as ECC

## Results (5-fold stratified CV)

| Dataset | Feature level | Classifier | Accuracy | AUROC | F1 |
|---|---|---|---|---|---|
| ECC | Species (970 taxa) | LR | **0.806** | **0.814** | **0.811** |
| ECC | Species (970 taxa) | RF | 0.750 | 0.799 | 0.780 |
| ECC | Strain (18 clusters) | LR | 0.556 | 0.585 | 0.556 |
| ECC | Strain (18 clusters) | RF | 0.611 | 0.591 | 0.611 |
| Lim_ORPI | Species (1601 taxa) | LR | 0.722 | 0.662 | 0.839 |
| Lim_ORPI | Species (1601 taxa) | RF | **0.722** | **0.723** | **0.839** |
| Lim_ORPI | Strain (9 clusters) | LR | 0.667 | 0.546 | 0.727 |
| Lim_ORPI | Strain (9 clusters) | RF | 0.611 | 0.515 | 0.696 |

## Key observations

1. **Species-level outperforms strain-level on both datasets.**
   - ECC: AUROC 0.81 (species) vs 0.59 (strain)
   - Lim_ORPI: AUROC 0.72 (species) vs 0.52 (strain)

2. **Very few strain calls were produced.**
   - ECC: 36 samples had strain calls, but only 18 distinct strain clusters passed prevalence filter (present in ≥2 samples).
   - Lim_ORPI: only 18 of 97 samples produced any strain calls; 9 distinct clusters passed prevalence filter.

3. **Most species were detected but not strain-resolved.**
   - Typical output per sample: "summary: X/30 species strain-resolved (0–2 strain calls), Y detected-not-resolvable, Z absent".
   - Example Lim_ORPI SF097: 9/30 species strain-resolved (2 strain calls), 9 detected-not-resolvable, 12 absent.

4. **BcgI single-enzyme marker density is the limiting factor.**
   - Even with `--min-support 2 --min-coverage 0.01` and `--layer1 cst`, *Streptococcus mutans* produced no strain calls in Lim_ORPI SF001.
   - Message: "no strain/cluster resolved — insufficient strain-specific 2b tags for this enzyme set".

## Why strain-level underperforms here

- **Marker starvation**: BcgI generates only a few hundred strain-specific tags per species. For many oral species this is below the support/coverage floor needed to distinguish strains in real metagenomes.
- **Database mismatch**: The kraken16 genome pool contains genomes assembled/annotated for other purposes; oral species coverage is uneven and some important taxa have only 1–3 genomes.
- **Feature sparsity**: With only 9–18 strain features, ML models cannot learn a stable decision boundary.
- **Species-level signal is already strong**: Both phenotypes (caries, denture hygiene) appear to be driven more by which species are present than by which strain of a species.

## Recommendations for improving strain-level AUROC

1. **Increase marker density**
   - Use multi-enzyme digital digestion (`--enzyme all` or `recommended`) on the 2bRAD-M reads if the protocol preserves full genomic fragments. This is the most direct fix: on test panels, `all` gives ~10–30× more markers than BcgI alone.
   - Alternatively run on original shotgun/WMS data rather than BcgI-enriched 2bRAD tags.

2. **Use k-mer mode for shotgun data**
   - `--marker-source kmer --sketch-scale 30` produces enough markers for CST internal nodes to carry signal. Only relevant if native shotgun reads are available.

3. **Expand and curate the strain reference panel**
   - Target oral-specific reference collections (e.g., HOMD, expanded GTDB oral representatives) rather than the broad kraken16 pool.
   - Ensure every target species has ≥10 high-quality genomes.

4. **Relax gates carefully**
   - The current default gates already produce many false-negative strain calls. Simply lowering thresholds will likely add false positives without recovering the true strain structure.
   - A sample-adaptive trace-gap (`--trace-gap 10 --trace-floor 1e-4`) may help on mock communities but is risky on open real samples.

5. **Feature engineering**
   - Treat "detected but not strain-resolved" species depth as a feature rather than discarding them. This is essentially species-level information but filtered through strain2bscan's species gate.

## Conclusion

On these two real BcgI-2bRAD oral datasets, **strain-level Strain2bScan profiling does not improve phenotype prediction over species-level abundance**. The bottleneck is BcgI marker density: the enzyme does not generate enough strain-discriminating tags for reliable strain calling in real metagenomes. Multi-enzyme or k-mer marker sources would be needed to test whether strain-level resolution can outperform species-level classification.
