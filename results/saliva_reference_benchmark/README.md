# Controlled saliva reference-panel comparison

## Question

The earlier pre-release Strain2bScan prototype (formerly Strain2bFunc) saliva analysis suggested that host identity could be
recovered more strongly when using a broad GTDB/HOMD-like reference. This run
isolates the reference-panel effect for Strain2bScan by keeping the reads,
BcgI marker space, profiling gates and downstream metrics fixed while changing
only the reference panel.

## Input

* 32 native BcgI 2bRAD saliva libraries (PRJNA1131785).
* Reads were the already-extracted BcgI FASTA files from
  `/lustre1/g/aos_shihuang/Strain2b/data/saliva_data/2bRAD_result/enzyme_result/`.
* The 32 sample aliases are identical to the manuscript saliva time course.
* Profiling used Strain2bScan commit `8fcb0bb`
  (`strainscan-port`) on HKU HPC.

The species list was derived from the archived broad saliva species table using
the prototype default selection rule: abundance >= 0.01% in >= 80% of the 32
samples. This gives 79 species
(`prototype_selected79_species.txt`). The broader unfiltered table has 403
species (`prototype_broad403_species.txt`).

## Panels

| Panel | Reference | Species DBs | Construction |
|---|---|---:|---|
| current oral panel | 19-species oral panel | 19 | existing manuscript panel; 18 species had strain calls |
| GTDB selected79 | GTDB r226 2bRAD-M tags | 157 | 3,072 selected genomes across 198 GTDB species; clustered at 0.95 Jaccard |
| HOMD selected79 | HROM/HOMD nonredundant genomes | 153 | 255 matched HROM species clusters; clustered at 0.95 Jaccard |

All GTDB/HOMD species DBs were built with `strain2bscan cluster --enzyme BcgI
--similarity 0.95`. Profiling used the same gates as the current saliva run:
`--min-species-markers 50 --min-species-detect 5 --min-support 8
--min-coverage 0.1 --min-consistency 0.5`.

## Overall host-discrimination result

Source: `saliva_reference_panel_comparison.tsv`.

| Panel | Level | Features | PERMANOVA R2 | p | LOO accuracy |
|---|---|---:|---:|---:|---:|
| current oral | species | 18 | 0.7553 | 0.0002 | 0.7813 |
| current oral | strain/cluster | 125 | 0.7567 | 0.0002 | 0.9688 |
| HOMD selected79 | species | 62 | 0.7000 | 0.0002 | 0.8125 |
| HOMD selected79 | strain/cluster | 618 | 0.6837 | 0.0002 | 1.0000 |
| GTDB selected79 | species | 76 | 0.7818 | 0.0002 | 0.9688 |
| GTDB selected79 | strain/cluster | 293 | 0.7675 | 0.0002 | 1.0000 |

Interpretation:

* Expanding from the current 19-species oral panel to GTDB selected79 raised
  overall species R2 from 0.7553 to 0.7818 and strain-cluster R2 from 0.7567
  to 0.7675. Thus, reference breadth does improve the Strain2bScan host signal.
* The improvement is modest, not as large as simply comparing a broad species
  table with a narrow strain panel.
* HOMD selected79 produced perfect LOO host classification but lower overall
  PERMANOVA R2, because the much larger feature set added noisy/sparse
  species-specific profiles.

## Prevotella and Rothia

Source: `saliva_reference_perspecies_subject.tsv`.

Key examples:

* Current oral panel: *Rothia mucilaginosa* R2 = 0.8521; only
  *Prevotella melaninogenica* was strain-resolved (R2 = 0.4640).
* GTDB selected79: *Prevotella intermedia* R2 = 0.9197;
  *Prevotella loescheii* R2 = 0.7001; *Rothia aeria* R2 = 0.7221;
  *Rothia dentocariosa* R2 = 0.7292.
* HOMD selected79: multiple additional Prevotella species were resolved,
  including *P. intermedia* R2 = 0.6824 and *P. melaninogenica* R2 = 0.6801.

This confirms that the earlier prototype observation was partly driven by
reference scope: broader oral/GTDB references recover additional
Prevotella/Rothia lineages with strong host-specific strain structure.

## Metrics

* Community matrices use Strain2bScan `sample_fraction`.
* Host effect is PERMANOVA R2 on Bray-Curtis distance with subject as the
  grouping factor (4,999 permutations, seed 42).
* LOO accuracy is leave-one-out 1-nearest-neighbour classification on the same
  Bray-Curtis distance matrix.
* Per-species analyses use within-species abundance (`abundance`) and require
  >= 10 present samples and >= 2 clusters.
* Per-species p values were BH corrected within each panel.

## Limitations

* This isolates reference-panel scope for Strain2bScan. It does not yet run the
  prototype genome-level elastic-net deconvolution on identical inputs.
* GTDB/HOMD species names were matched by normalized nomenclature, including
  GTDB placeholder-suffix fallbacks. Some oral-taxonomy or SGB-style names do
  not map exactly and remain absent.
* Species-gate thresholds were held constant for comparability, although broad
  panels may benefit from separate gate calibration.
* A maximum p value of 0.0002 reflects 4,999 permutations plus the observed
  statistic.

## Provenance

* HKU HPC project: `/lustre1/g/aos_shihuang/projects/strain2bscan_reference_benchmark_20261120`
* GTDB tags: `/lustre1/g/aos_shihuang/databases/GTDB/GTDBr226/2bradm_gtdb_db/db/BcgI.species.fa.gz`
* GTDB mapping: `/lustre1/g/aos_shihuang/databases/GTDB/GTDBr226/2bradm_gtdb_db/db/abfh_classify_with_speciename.txt.gz`
* HROM genomes: `/lustre1/g/aos_shihuang/databases/HROM/data/genome_catalog/HROM_nonredundant_genomes`
* HROM metadata: `/lustre1/g/aos_shihuang/databases/HROM/data/genome_catalog/HROM_representative_genome_metadata.tsv`
* GTDB build job: `4476964` (completed)
* HOMD build/profile job: `4476568` (completed)
