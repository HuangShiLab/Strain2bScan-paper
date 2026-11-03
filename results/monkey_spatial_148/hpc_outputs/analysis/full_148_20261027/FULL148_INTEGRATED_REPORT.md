# Full 148-sample monkey strain-resolved spatial analysis

## Cohort definition

The full cohort contains **148 samples**:

| Group | n | Description |
|---|---:|---|
| GI spatial samples | 128 | 66 luminal + 62 mucosal |
| Body/control samples | 20 | monkey-matched B_* body-site controls |
| **Total** | **148** | |

Four additional vendor files (`O1_Tg_F_Th`, `S3779425`, `Y1_Th_TG_f`, `Y2_Tg_TR`) were excluded because they are not part of the 148-sample GI/body spatial cohort.

## Data modality coverage

| Analysis tier | n |
|---|---:|
| GI 2bRAD only | 62 |
| GI 2bRAD + shotgun + HUMAnN function | 30 |
| GI multimodal 2bRAD + shotgun + function + bile acid | 36 |
| Body/control 2bRAD | 20 |

Thus, all 128 GI samples have native 2bRAD strain profiling. The 66 luminal GI samples also have shotgun/HUMAnN profiles; 36 of them additionally have bile-acid data.

## Profiling method

- Strain2bScan branch: `strainscan-port`
- Panel: 18 monkey-gut species represented by 315 study MAGs and 303 0.95 strain-cluster units
- Marker space: BcgI for both native 2bRAD and shotgun
- Cohort profiled: all 148 samples with identical gates
- Call gates:
  - species panel markers >= 20
  - detected species markers >= 2
  - unique cluster-supporting markers >= 2
  - cluster coverage >= 0.01
  - consistency >= 0.5

## Sample-level strain recovery

| Group | Samples profiled | Samples with >=1 strain call | Positive rate | Median calls/sample |
|---|---:|---:|---:|---:|
| GI spatial | 128 | 113 | 88.3% | 13.5 |
| Luminal GI | 66 | 56 | 84.8% | — |
| Mucosal GI | 62 | 57 | 91.9% | — |
| Body/control | 20 | 2 | 10.0% | 0 |

Only three body-control strain calls were observed, indicating that the main GI strain signals are not explained by widespread body-control contamination.

## Main species-level signals

Top GI species by strain-call prevalence:

| Species | GI prevalence | Luminal prevalence | Mucosal prevalence | Body-control prevalence |
|---|---:|---:|---:|---:|
| *Helicobacter_B_sp003364135* | 76.6% | 72.7% | 80.6% | 10.0% |
| *UBA4372_sp900543815* | 66.4% | 72.7% | 59.7% | 0% |
| *Cryptobacteroides_sp000433355* | 65.6% | 77.3% | 53.2% | 5.0% |
| *Succinivibrio_sp000431835* | 57.8% | 65.2% | 50.0% | 0% |
| *WG_1_sp900539665* | 55.5% | 66.7% | 43.5% | 0% |
| *Oribacterium_sp900768305* | 50.0% | 53.0% | 46.8% | 0% |
| *CAG_873_sp001701165* | 49.2% | 59.1% | 38.7% | 0% |
| *Lachnospira_eligens_A* | 47.7% | 53.0% | 41.9% | 0% |
| *Treponema_D_succinifaciens* | 44.5% | 54.5% | 33.9% | 0% |
| *Treponema_D_berlinense* | 44.5% | 54.5% | 33.9% | 0% |
| *Ligilactobacillus_animalis* | 42.2% | 54.5% | 29.0% | 0% |
| *Cryptobacteroides_sp900546445* | 42.2% | 48.5% | 35.5% | 0% |

The body-control signal is low. No top species or cluster passed the prespecified potential-background criterion.

## Spatial strain structure

Across all GI samples, **89 strain clusters** were classified as LI-enriched using:

- peak-region prevalence >= 0.3,
- region specificity >= 0.8,
- peak region = LI.

Top LI-enriched clusters included:

| Species | Cluster | GI prevalence | Luminal prevalence | Mucosal prevalence | LI peak prevalence | Region specificity |
|---|---|---:|---:|---:|---:|---:|
| *Treponema_D_succinifaciens* | C25 | 32.8% | 40.9% | 24.2% | 75.0% | 1.00 |
| *WG_1_sp900539665* | C22 | 32.8% | 43.9% | 21.0% | 71.4% | 0.94 |
| *Treponema_D_succinifaciens* | C2 | 28.9% | 42.4% | 14.5% | 64.3% | 0.96 |
| *Succinivibrio_sp000431835* | C6 | 26.6% | 36.4% | 16.1% | 58.9% | 0.96 |
| *Cryptobacteroides_sp900546445* | C1 | 25.8% | 36.4% | 14.5% | 57.1% | 0.96 |
| *Treponema_D_succinifaciens* | C15 | 24.2% | 39.4% | 8.1% | 55.4% | 1.00 |
| *UBA636_sp900546285* | C12 | 23.4% | 30.3% | 16.1% | 53.6% | 1.00 |
| *Cryptobacteroides_sp000433355* | C6 | 25.8% | 43.9% | 6.5% | 53.6% | 0.88 |

These clusters were absent from body controls and therefore represent a genuine distal-gut luminal strain reservoir rather than body-site contamination.

## Mucosal versus luminal partitioning

Matched mucosal/luminal pairs were tested within the same monkey, GI region and position. Absence was represented as zero abundance.

- 153 strain clusters had enough positive support for a Wilcoxon paired test.
- 20 were nominal p < 0.05.
- **13 clusters remained significant after Bonferroni correction.**
- All 13 were **luminal-enriched**.

Bonferroni-significant luminal-enriched clusters:

| Species | Cluster | Matched pairs | Mucosal positive | Luminal positive | FDR/Bonferroni p |
|---|---|---:|---:|---:|---:|
| *Treponema_D_succinifaciens* | C25 | 50 | 15 | 27 | 9.6e-4 |
| *UBA4372_sp900543815* | C13 | 50 | 4 | 24 | 2.8e-3 |
| *Treponema_D_succinifaciens* | C6 | 50 | 7 | 23 | 4.1e-3 |
| *Treponema_D_succinifaciens* | C22 | 50 | 8 | 22 | 4.7e-3 |
| *UBA4372_sp900543815* | C5 | 50 | 7 | 32 | 6.5e-3 |
| *Cryptobacteroides_sp900546445* | C13 | 50 | 4 | 24 | 6.7e-3 |
| *Cryptobacteroides_sp900546445* | C3 | 50 | 4 | 27 | 1.3e-2 |
| *WG_1_sp900539665* | C20 | 50 | 16 | 30 | 1.8e-2 |
| *Cryptobacteroides_sp900546445* | C14 | 50 | 4 | 22 | 1.8e-2 |
| *Cryptobacteroides_sp000433355* | C3 | 50 | 5 | 28 | 3.1e-2 |
| *Treponema_D_succinifaciens* | C10 | 50 | 9 | 20 | 3.3e-2 |
| *Treponema_D_succinifaciens* | C12 | 50 | 4 | 21 | 3.3e-2 |
| *Treponema_D_berlinense* | C10 | 50 | 4 | 17 | 4.5e-2 |

This is strong evidence that the distal-gut luminal compartment contains a distinct strain reservoir.

## Integration with bile acids and function

The 148-sample 2bRAD layer was integrated with prior shotgun/HUMAnN/bile-acid overlays:

- 66 luminal GI samples have shotgun strain profiling and HUMAnN function.
- 36 luminal GI samples have targeted bile-acid SRM data.
- In the LI subset, strain composition coupled to bile-acid composition and selected functional pathways.

Key coupling observations from the LI analysis:

1. **Bile acids**  
   LI-luminal strain reservoir abundance was coupled to a bile-acid composition depleted in conjugated/high-mass bile acids and enriched in lower-mass/free bile-acid classes.

2. **SCFA and fermentation**  
   *UBA4372_sp900543815* C4 was positively coupled to pyruvate-to-butyrate and acetyl-CoA-to-butyrate pathways. Other clusters were positively coupled to propionate programs.

3. **Mucus/glycan-associated functions**  
   *Treponema_D_succinifaciens* C2 and *Succinivibrio_sp000431835* C4 were strongly coupled to N-acetylneuraminate and N-acetylglucosamine/mannosamine/neuraminate degradation pathways, consistent with host-glycan-derived sugar utilization.

## Main biological conclusion

The full 148-sample analysis establishes a robust strain-resolved spatial model of the monkey gut:

> The monkey gut is organized at strain-cluster resolution, not merely at species resolution. A large set of 0.95 strain clusters forms a LI-luminal reservoir, is absent from body controls, and is significantly luminal-enriched in matched GI positions. In shotgun/HUMAnN/bile-acid subsets, this reservoir couples to bile-acid transformation, SCFA fermentation programs and host-glycan-associated metabolic functions.

## Key outputs

- Full metadata: `full_148_sample_metadata.tsv`
- Full 148 strain-call long table: `full148_2brad_strain_calls_long.tsv`
- Sample-level summary: `full148_sample_strain_summary.tsv`
- Integrated modality master table: `integrated_sample_master_148.tsv`
- GI/control species summary: `full148_species_GI_control_summary.tsv`
- Cluster region/control summary: `full148_2brad_cluster_GI_control_region_summary.tsv`
- Mucosal/luminal paired tests: `full148_mucosal_luminal_paired_tests.tsv`
- Full report: `FULL148_INTEGRATED_REPORT.md`
