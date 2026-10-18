# Strain2bScan manuscript structure after Shi Huang's 2026-10-18 edit

## Target narrative

The Results section now follows two input modes. Within each input mode, evidence moves from
in-silico simulation to controlled DNA mock to real metagenome.

### Shared foundation

1. **Figure 1** — Strain2bScan overview and the two input modes.
2. **Figure 2** — 2bRAD markers retain genome-wide strain relationships; 16S does not.

### Part I — native BcgI 2bRAD-M performance

3. **Figure 3** — in-silico enzyme/marker-density simulation across 14 resolvable species from the
   15-species pool. This is a design experiment; BcgI is the single-enzyme condition used for native
   2bRAD-M libraries.
4. **Figure 4** — DNA mocks MSA-1002, MSA-1003, MSA-1005 and MSA-1007. The 99.9%-host row has been
   removed from the primary figure.
5. **Figure 5** — real saliva time-course: 32 native 2bRAD-M libraries from eight subjects.
6. **Figure 6** — real paired saliva comparison: native 2bRAD-M versus shotgun WMS prefix libraries.
   The exploratory four-library oral cohort remains in Table S3.

### Part II — conventional shotgun / in-silico digestion performance

7. **Figure 7** — complete 2,025-sample single-species simulation pool. This establishes the
   Strain2bScan accuracy envelope; the archived provenance is explicitly disclosed.
8. **Figure 8** — matched 225-sample simulation subset used for the current-source multi-tool
   comparison. Strain2bScan and StrainScan are scored in their own 0.95 spaces; StrainGST is shown
   separately in its own 0.90 space. Figure 8 now includes StrainGST accuracy and runtime.
9. **Figure 9** — computational scaling and community throughput.
10. **Figure 10** — ATCC shotgun DNA mocks with Strain2bScan, StrainScan and inStrain.
11. **Tables 4–7** — real shotgun metagenomes and isolate WGS application cohorts.

### Supplementary

- **Figure S6** — former main Figure 4: Jaccard versus containment reference-completeness
  optimization. This is now framed as an algorithmic optimization analysis rather than a primary
  performance benchmark.

## Comment-to-action status

| ShiEdit comment | Status |
|---|---|
| C5: add simulation/mock comparisons to other tools | Partially addressed. StrainGST is already complete on 225 simulations, 12 communities and four shotgun mocks, and is now shown in Figure 8. inStrain has only six exploratory single-species rows and three community rows in `results/instrain_sim.tsv`; it is not yet a complete 225-sample comparison. StrainPhlAn has not been run. |
| C82: use 32 paired shotgun metagenomes rather than four | NCBI confirms 32 paired WMS saliva runs (`WMS_S9_*`, `WMS_S11_*`, `WMS_S13_*`, `WMS_S17_*`) plus controls. The current manuscript correctly says that only four libraries were analyzed because only prefix FASTQs were available. A full 32-library rerun is a required new benchmark. |
| C94: clarify whether Figure 3 was shotgun or 2bRAD | Clarified. Simulation FASTQs are ART error-modelled shotgun reads followed by in-silico 2bRAD digestion; they are not native BcgI 2bRAD libraries. |
| C97: replace three-species Figure 3 with the 15-species pool | Main-text Figure 3 now uses the enzyme simulation from the 15-species pool. The complete 2,025-sample Strain2bScan simulation is new Figure 7. The n=2025 table is archived, so a current-source n=2025 rerun is still needed before it can be treated as final primary evidence. |
| C98: move Jaccard/containment work to supplement | Done. Former Figure 4 is now Figure S6. |
| C104/C112: include more tools in Figures 9/10 | Figure 8 now includes StrainGST. Figure 9 remains an efficiency benchmark and Figure 10 already includes inStrain. A complete inStrain/StrainPhlAn simulation benchmark still requires new runs. |

## Required benchmark extensions before final submission

1. **Current-source n=2025 rerun.** Regenerate or restore the 2,025 FASTQs and rerun the current
   Strain2bScan commit. The existing n=2025 table is historical and differs from current-source
   output on 101 of the 225 overlapping runs.
2. **Thirty-two paired saliva shotgun runs.** Replace the four-prefix analysis with a uniform
   analysis of all 32 WMS saliva runs (or a prespecified equal-depth prefix rule applied to all 32).
3. **Complete inStrain simulation benchmark.** Extend beyond the existing nine exploratory rows and
   define reference/dereplication and cluster-space scoring before comparing accuracy.
4. **StrainPhlAn benchmark.** Run StrainPhlAn on the same simulation reads or explain that it
   estimates strain types rather than the same cluster abundance task; do not force its output into
   P/R/F1 without a transparent mapping.
5. **Figure 5/enzyme-sweep provenance.** Preserve or regenerate the exact number of samples and
   strains behind each species × enzyme aggregate; the current table contains 70 aggregate rows but
   not per-sample counts.
