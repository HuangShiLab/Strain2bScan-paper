# Cohort-specific isolate panel validation

Use PRJNA1191225 isolate reads to assemble study-specific genomes, build an E. coli / Bifidobacterium panel, validate isolate reads, and reprofile P08 W1-W3 metagenomes.

## Assemblies

| assembly | contigs | total_bp | largest_bp | N50 |
|---|---|---|---|---|
| Bifidobacterium_bifidum_LHCA82 | 180 | 2299850 | 170675 | 75562 |
| Bifidobacterium_breve_LHCA81 | 90 | 2393363 | 286754 | 228295 |
| Bifidobacterium_longum_LHCA43 | 164 | 2355114 | 172561 | 55550 |
| Escherichia_coli_LHCA45 | 185 | 4962385 | 640837 | 239924 |
| Escherichia_coli_LHCA56 | 448 | 5147585 | 519210 | 222606 |
| Escherichia_coli_LHCA72 | 457 | 5040433 | 365311 | 159193 |

## Panel

- E. coli: 3 assemblies clustered into 3 resolvable clusters; 47,410 total and 19,122 unique markers.
- B. bifidum LHCA82: 17,985 markers.
- B. breve LHCA81: 19,214 markers.
- B. longum LHCA43: 18,426 markers.

## Isolate self-validation: 6/6 correct

| sample | expected | observed | correct | coverage | support | depth | sample_fraction |
|---|---|---|---|---|---|---|---|
| LHCA45 | Escherichia_coli|C0 | Escherichia_coli|C0 | True | 1.0000 | 7842 | 25.1572 | 0.964047 |
| LHCA56 | Escherichia_coli|C1 | Escherichia_coli|C1 | True | 1.0000 | 5533 | 25.6356 | 0.796484 |
| LHCA72 | Escherichia_coli|C2 | Escherichia_coli|C2 | True | 0.9998 | 5746 | 16.3273 | 0.816417 |
| B_longum_LHCA43 | Bifidobacterium_longum|Bifidobacterium_longum_LHCA43 | Bifidobacterium_longum|Bifidobacterium_longum_LHCA43 | True | 0.9995 | 15949 | 9.7554 | 0.799987 |
| B_breve_LHCA81 | Bifidobacterium_breve|Bifidobacterium_breve_LHCA81 | Bifidobacterium_breve|Bifidobacterium_breve_LHCA81 | True | 1.0000 | 17029 | 24.7082 | 0.835125 |
| B_bifidum_LHCA82 | Bifidobacterium_bifidum|Bifidobacterium_bifidum_LHCA82 | Bifidobacterium_bifidum|Bifidobacterium_bifidum_LHCA82 | True | 1.0000 | 17311 | 22.8871 | 0.812597 |

Runtime: 9.53 s for six isolate read sets; peak RSS 55.7 MiB.

## P08 longitudinal result with cohort panel

| sample | species | cluster | abundance | coverage | support | depth | global_abundance | sample_fraction | n_markers |
|---|---|---|---|---|---|---|---|---|---|
| P08_W1 | Bifidobacterium_bifidum | Bifidobacterium_bifidum_LHCA82 | 1.000000 | 0.7638 | 12915 | 236.6885 | 1.000000 | 0.645624 | 17985 |
| P08_W2 | Bifidobacterium_bifidum | Bifidobacterium_bifidum_LHCA82 | 1.000000 | 0.7734 | 12938 | 213.8555 | 0.999865 | 0.583400 | 17985 |
| P08_W2 | Escherichia_coli | C0|C1|C2 | 1.000000 | 0.0258 | 492 | 0.0289 | 0.000135 | 0.000140 | 31877 |
| P08_W3 | Bifidobacterium_bifidum | Bifidobacterium_bifidum_LHCA82 | 1.000000 | 0.7923 | 13019 | 238.1707 | 1.000000 | 0.609940 | 17985 |

Key result:
- Bifidobacterium bifidum LHCA82 persists across W1-W3 with breadth 0.764/0.773/0.792 and depth 237x/214x/238x.
- Within-panel abundance remains high: 64.6%, 58.3%, 61.0%.
- E. coli appears only at W2, at depth 0.029x and fraction 1.40e-4; the panel reports unresolved C0|C1|C2 instead of falsely selecting one isolate.

## Generic MSA panel vs cohort panel

| week | generic_panel_calls | cohort_panel_calls |
|---|---|---|
| W1 | Staphylococcus_aureus|C3@frac=6.50e-05 | Bifidobacterium_bifidum|Bifidobacterium_bifidum_LHCA82@frac=6.46e-01 |
| W2 | Escherichia_coli|C0@frac=6.10e-05;Escherichia_coli|C2@frac=5.70e-05 | Bifidobacterium_bifidum|Bifidobacterium_bifidum_LHCA82@frac=5.83e-01;Escherichia_coli|C0|C1|C2@frac=1.40e-04 |
| W3 | Enterococcus_faecalis|C5@frac=1.20e-04;Enterococcus_faecalis|C1@frac=9.00e-05;Enterococcus_faecalis|C3@frac=7.10e-05;Enterococcus_faecalis|C0@frac=5.00e-05;Enterococcus_faecalis|C4@frac=4.40e-05;Enterococcus_faecalis|C2@frac=3.80e-05 | Bifidobacterium_bifidum|Bifidobacterium_bifidum_LHCA82@frac=6.10e-01 |

The generic panel missed the dominant B. bifidum strain because its Bifidobacterium DB was B. adolescentis. This demonstrates why cohort-specific panels matter.

## Efficiency

- Isolate validation: 6 samples in 9.53 s, peak 55.7 MiB.
- P08 cohort-panel run: 3 samples in 52.49 s, peak 77.5 MiB.

## Files

isolate_assembly_stats.tsv; isolate_validation.scorecard.tsv; P08_cohortpanel_calls.tsv; P08_generic_vs_cohort.tsv; cohort_panel_run_summary.tsv
