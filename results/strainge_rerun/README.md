# StrainGE rerun provenance

## Scope

This rerun adds the StrainGST component of StrainGE to the existing shotgun
benchmark. It covers:

1. the 225 published matched single-species simulation runs (`_rep1_` only);
2. the 12 multi-species simulation communities;
3. four primary ATCC DNA-mock WMS libraries: MSA-1002 at 99% host DNA,
   MSA-1003, MSA-1005 and MSA-1007.

Native BcgI 2bRAD reads were not submitted to StrainGE because its standard
workflow expects WMS reads. This is therefore a shotgun-input comparison. The
StrainGR variant-calling component was not run; these are strain-reference
identification and abundance results from StrainGST.

## Software

- Tool: StrainGE 1.3.9 (`straingst`/StrainGST)
- Source snapshot used for documentation: Broad Institute `StrainGE`
  commit `86de923dba7f70b25cf1f887a4873fa4eae63eeb`
- Database workflow: k-merize each panel genome; remove near-subsets; cluster
  at Jaccard 0.90; create a species-level StrainGST pan-genome database.
- Sample workflow: k-merize paired reads once, then run StrainGST against each
  relevant species database.
- StrainGST iterations: 5 for single-species simulations, 8 for multi-species
  simulations, and 32 for ATCC DNA mocks.
- Identification calls use StrainGST score >= 0.02. DNA-mock detection metrics
  additionally use relative abundance >= 1e-4, matching the primary mock
  threshold.

## Reference and truth spaces

Each tool is scored in its own reference space. Truth genomes and StrainGST
reported references are mapped through the StrainGE 0.90 Jaccard clusters.
No Strain2bScan or StrainScan cluster IDs are reused for StrainGE.

The simulation panel contained one empty NCBI FASTA
(`Klebsiella_pneumoniae/GCA_057929645.fna`); it is recorded in
`skipped_empty_simulation_genomes.tsv` and was not part of the simulation truth.

## Outputs

- `simulation_database_build.tsv`: per-species StrainGE database build cost.
- `mock_database_build.tsv`: per-species ATCC-panel database build cost.
- `simulation_run_efficiency.tsv`: per-sample sample k-merization and summed
  per-species StrainGST search wall time and peak RSS.
- `mock_run_efficiency.tsv`: the same efficiency fields for the four DNA mocks.
- `simulation_persample.tsv` and derived summaries: strain-reference precision,
  recall and F1 in the StrainGE 0.90-reference space.
- `mock_persample.tsv`: DNA-mock precision, recall, F1 at 1e-4, and abundance
  profile similarity.
- `manuscript_summary.tsv`: source-derived values used in manuscript Table 9.

## Headline results

- 225/225 matched single-species simulations completed. Median precision,
  recall and F1 were 1.000 in the StrainGE 0.90-reference space; median wall
  time was 18.66 s/sample and median peak RSS was 1.30 GB.
- All 12 multi-species communities completed. Median precision, recall and F1
  were 0.901, 1.000 and 0.941; median wall time was 236.98 s/sample and median
  peak RSS was 7.90 GB.
- On the same matched simulations, Strain2bScan medians were 0.35 s and
  0.051 GB (single-species) and 4.33 s and 0.654 GB (multi-species). These
  runtime comparisons are in `runtime_comparison_summary.tsv`; accuracy is not
  pooled because StrainGE and Strain2bScan use different cluster spaces.
- All four ATCC WMS mocks completed. F1 was 0.700, 0.542, 0.476 and 0.600 for
  MSA-1002 99% host, MSA-1003, MSA-1005 and MSA-1007, respectively. The wall
  time includes one sample k-merization plus 28 species-database searches
  (1007.22-2158.75 s); peak RSS was 14.36-21.52 GB.
