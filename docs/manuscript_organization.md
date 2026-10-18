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
6. **Figure 6** — real paired saliva comparison: native 2bRAD-M versus shotgun WMS prefix libraries. A complete canonical 32-sample WMS rebuild is now available for follow-up benchmarking. The exploratory four-library oral cohort remains in Table S3.

### Part II — conventional shotgun / in-silico digestion performance

7. **Figure 7** — complete current-source 2,025-sample single-species simulation pool. This establishes the Strain2bScan accuracy envelope.
8. **Figure 8** — matched 225-sample simulation subset used for the current-source multi-tool comparison. Strain2bScan and StrainScan use 0.95 spaces, StrainGST uses a 0.90 space, and inStrain is genome-resolved then mapped to truth 0.95 clusters. Accuracy is not pooled across reference spaces.
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
| C5: add simulation/mock comparisons to other tools | Complete for executable direct comparators. StrainGST is complete on 225 simulations, 12 communities and four shotgun mocks. inStrain is complete on all 225 matched simulations (`results/instrain_sim225/`). Both are shown in Figure 8. StrainPhlAn was not forced into the P/R/F1 framework because it infers strain phylogeny/types rather than cluster abundance calls. |
| C82: use 32 paired shotgun metagenomes rather than four | Complete. All 32 paired WMS samples now have canonical gzipped mate files in `/lustre1/g/aos_shihuang/Strain2b/data/saliva_data/WMS_data/complete_32_fqgz`; `results/wms32_rebuild/` records the 64-file manifest and completion evidence. |
| C94: clarify whether Figure 3 was shotgun or 2bRAD | Clarified. Simulation FASTQs are ART error-modelled shotgun reads followed by in-silico 2bRAD digestion; they are not native BcgI 2bRAD libraries. |
| C97: replace three-species Figure 3 with the 15-species pool | Main-text Figure 3 now uses the enzyme simulation from the 15-species pool. The complete current-source 2,025-sample Strain2bScan rerun is Figure 7 and `results/sim2025_current/`. |
| C98: move Jaccard/containment work to supplement | Done. Former Figure 4 is now Figure S6. |
| C104/C112: include more tools in Figures 9/10 | Figure 8 now includes StrainGST and the complete inStrain 225-sample benchmark. Figure 9 remains an efficiency benchmark and Figure 10 already includes inStrain. StrainPhlAn requires a separate phylogeny-placement metric rather than direct P/R/F1 pooling. |

## Required benchmark extensions before final submission

1. **Current-source n=2025 rerun.** Complete. `results/sim2025_current/` contains all 2,025 runs at source commit `6ae760c2f3bd2f39efc0e6aa06c53a84366cb16e`; it supersedes the archived table.
2. **Thirty-two paired saliva shotgun runs.** Complete. All 64 canonical mate files are present in `/lustre1/g/aos_shihuang/Strain2b/data/saliva_data/WMS_data/complete_32_fqgz` and recorded in `results/wms32_rebuild/`.
3. **Complete inStrain simulation benchmark.** Complete. `results/instrain_sim225/persample.tsv` contains all 225 matched runs; summaries are in `by_depth.tsv` and `by_species.tsv`.
4. **StrainPhlAn scope.** Deferred with an explicit scope statement: it infers strain phylogeny/types and does not emit cluster abundance calls directly comparable to the P/R/F1 benchmark.
5. **Figure 5/enzyme-sweep provenance.** Preserve or regenerate the exact number of samples and
   strains behind each species × enzyme aggregate; the current table contains 70 aggregate rows but
   not per-sample counts.
