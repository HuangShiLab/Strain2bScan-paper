# Strain-cluster coupling to bile acids, SCFA fermentation and mucus/glycan functions

## Question

Do the distal-gut luminal strain clusters identified by Strain2bScan couple to targeted bile-acid features, short-chain-fatty-acid fermentation pathways, or mucus/glycan utilization functions?

## Data processed

### Strain profiles

- Shotgun Strain2bScan calls from `work/monkey_panel_20261023/profile_shotgun/merged.tsv`.
- Primary strain variable: **within-species strain-cluster composition** (`abundance`), not total strain-call load.
- Large-intestine luminal samples were used for the distal-gut test.
- Clusters detected in at least eight LI samples were tested: **219 strain clusters**.

### Targeted bile acids

Because only vendor mzML files and no quantification sheet were present, SRM chromatograms were parsed directly from mzML.

- 36 biological LI/SI samples passed parsing.
- 7,675 SRM chromatograms were peak-parsed.
- 176 internal-standard chromatograms were used for normalization.
- Peak areas were normalized to internal-standard signal and transformed as `log2(area+1)`.
- Negative-mode transitions were grouped by precursor mass into four interpretable bile-acid feature classes:
  - `LCA_like_Q1_375`
  - `dihydroxy_Q1_391`
  - `cholate_like_Q1_407`
  - `conjugated_or_highmass`
- For coupling tests, bile-acid signals were converted to within-sample composition so that tests did not merely reflect total BA signal.
- Final LI strain/BA overlap: **23 samples**.

### Functional profiles

- HUMAnN `pathabundance-cpm` tables were aggregated across 70 kneaddata inputs.
- 470 unstratified MetaCyc pathways were recovered.
- Selected pathways covered:
  - butyrate/SCFA fermentation,
  - lactate/propionate/mixed-acid fermentation,
  - sialic-acid and N-acetyl-sugar degradation,
  - fucose/rhamnose degradation,
  - glucarate/galactarate degradation.
- Final LI strain/function overlap: **26 samples**.

## Statistical model

For each strain-cluster/feature pair, both variables were residualized against:

1. monkey identity,
2. LI position,

then correlated using Spearman correlation of residuals. Benjamini–Hochberg FDR was calculated separately for bile-acid and function tests.

This asks whether strain composition varies with BA or function composition beyond host background and position.

## Bile-acid coupling

### Primary compositional result

- 876 strain-cluster × BA-category tests.
- 410 had FDR q < 0.05.
- 409 significant associations were negative and one was positive.

The dominant pattern is therefore not a single BA specialist, but a coordinated shift:

> Samples dominated by the LI-luminal strain reservoir have a bile-acid composition depleted in conjugated/high-mass bile acids and relatively enriched in lower-mass/free bile-acid classes.

Key examples:

| Strain cluster | BA class | rho | p | FDR q | Direction |
|---|---|---:|---:|---:|---|
| *Cryptobacteroides_sp000433355* C1 | conjugated/high-mass | -0.943 | 1.8e-11 | 1.0e-8 | negative |
| *Cryptobacteroides_sp000433355* C5 | cholate-like Q1 407 | -0.935 | 6.6e-11 | 1.4e-8 | negative |
| *Treponema_D_succinifaciens* C25 | conjugated/high-mass | -0.923 | 3.6e-10 | 5.2e-8 | negative |
| *Succinivibrio_sp000431835* C25 | cholate-like Q1 407 | -0.909 | 1.9e-9 | 1.4e-7 | negative |
| *UBA4372_sp900543815* C3 | conjugated/high-mass | -0.895 | 8.1e-9 | 2.7e-7 | negative |
| *UBA636_sp900546285* C11 | dihydroxy Q1 391 | +0.648 | 8.2e-4 | 2.5e-3 | positive |

### Interpretation

The distal luminal reservoir is coupled to a bile-acid milieu that is relatively depleted of conjugated/high-mass bile acids. This is compatible with strong luminal bile-acid transformation in the LI, but the current data do not prove that these clusters directly perform the conversions.

### Direct gene-potential caveat

The available MAG scans included broad 3β-hydroxysteroid dehydrogenase PF hits, but did not provide a cluster-resolved BSH/bai operon catalogue. Therefore, bile-acid coupling should be described as **strain-phenotype association**, not direct pathway assignment to a strain cluster.

## SCFA and fermentation coupling

Functional tests used 219 strain clusters × 14 selected pathways = 3,066 tests.

- 1,233 associations had FDR q < 0.05.
- 1,019 were positive and 214 negative.

Key positive associations:

| Strain cluster | Pathway | rho | p | FDR q | Interpretation |
|---|---|---:|---:|---:|---|
| *UBA4372_sp900543815* C4 | pyruvate fermentation to butanoate | 0.926 | 1.2e-11 | 4.0e-10 | strong positive butyrate-pathway coupling |
| *UBA4372_sp900543815* C4 | acetyl-CoA fermentation to butanoate II | 0.905 | 2.2e-10 | 4.9e-9 | strong positive butyrate-pathway coupling |
| *Oribacterium_sp900768305* C2 | pyruvate fermentation to butanoate | 0.955 | 3.7e-14 | 2.9e-12 | strong positive butyrate-pathway coupling |
| *Cryptobacteroides_sp000433355* C5 | pyruvate fermentation to propanoate I | 0.981 | 1.5e-18 | 5.6e-16 | propionate-associated composition |

Key negative association:

| Strain cluster | Pathway | rho | p | FDR q | Interpretation |
|---|---|---:|---:|---:|---|
| *UBA4372_sp900543815* C3 | succinate fermentation to butanoate | -0.951 | 1.0e-13 | 7.1e-12 | alternative butyrate route anti-correlated |

Interpretation:

The LI reservoir is coupled to a mature distal-gut fermentation state. Some dominant clusters—especially *UBA4372* C4 and *Oribacterium* C2—track pyruvate-to-butyrate pathways. Other clusters track propionate programs. This supports strain-level differentiation of central fermentation routes.

## Mucus/glycan and host-glycan-derived sugar coupling

Several LI strain clusters were strongly coupled to host-glycan-associated sugar degradation pathways.

| Strain cluster | Pathway | rho | p | FDR q | Interpretation |
|---|---|---:|---:|---:|---|
| *Treponema_D_succinifaciens* C2 | superpathway of N-acetylneuraminate degradation | 0.964 | 2.2e-15 | 2.6e-13 | sialic-acid utilization coupling |
| *Treponema_D_succinifaciens* C2 | superpathway of N-acetylglucosamine/mannosamine/neuraminate degradation | 0.965 | 2.1e-15 | 2.1e-13 | mucus-glycan-derived sugar utilization coupling |
| *Succinivibrio_sp000431835* C4 | superpathway of N-acetylglucosamine/mannosamine/neuraminate degradation | 0.789 | 2.0e-6 | 1.2e-5 | mucus-glycan-derived sugar utilization coupling |

This is important because sialic acid and N-acetyl sugars are compatible with mucin O-glycan degradation. It suggests that part of the LI luminal reservoir may depend on, or co-occur with, host-glycan-derived sugar processing rather than solely dietary fiber utilization.

## Integrated model

The coupling results refine the strain-resolved spatial model:

1. **Distal-gut luminal strain reservoir**  
   LI-dominant clusters in *Treponema*, *Cryptobacteroides*, *Succinivibrio*, *WG-1*, *UBA4372* and related lineages form a coherent luminal community.

2. **Bile-acid milieu**  
   This reservoir is associated with a bile-acid profile depleted in conjugated/high-mass BA and enriched in lower-mass/free BA classes.

3. **Fermentation programs**  
   Strain clusters partition across butyrate, propionate and mixed-acid fermentation states. *UBA4372* C4 is a central candidate for pyruvate/acetyl-CoA-to-butyrate coupling.

4. **Mucus-glycan interface**  
   LI clusters such as *Treponema_D_succinifaciens* C2 and *Succinivibrio* C4 are coupled to sialic-acid/N-acetyl-sugar degradation, suggesting a luminal community sustained partly by host glycan-derived sugars.

5. **Functional interpretation**  
   These are sample-level strain–function couplings, not proof that a pathway gene belongs to a particular cluster. Follow-up work should map cluster-specific MAG genes to BSH/bai, SCFA and mucin-glycan enzymes.

## Main conclusion

> The distal luminal strain reservoir is not independent of gut chemistry: its abundance is coupled to a free-BA-enriched bile-acid profile, a mature SCFA/fermentation program, and host-glycan-derived sugar utilization pathways. This links strain-level spatial structure to bile-acid, SCFA and mucus-associated metabolic environments.

## Key files

- Bile-acid peak metrics: `bileacid_transition_features.tsv`
- BA category matrix: `bileacid_category_log2_matrix.tsv`
- BA target-group matrix: `bileacid_group_log2_matrix.tsv`
- HUMAnN unstratified pathway matrix: `humann_pathabundance_cpm_unstratified.tsv`
- Selected function matrix: `humann_selected_function_matrix.tsv`
- Absolute-abundance coupling: `coupling_bileacid_categories_li.tsv`, `coupling_functions_li.tsv`
- Compositional coupling:  
  - `coupling_compositional_bileacid_li.tsv`  
  - `coupling_compositional_functions_li.tsv`
- Selected key results: `key_coupling_results.tsv`
- BA-gene potential scan: `strain_cluster_ba_gene_potential.tsv`
