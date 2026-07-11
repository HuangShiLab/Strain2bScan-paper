# Strain2bScan — Figure legends and file map

Figure numbering follows `docs/manuscript_organization.md` (two-input-mode structure). Each entry gives
the figure file(s) in `figures/`, the source table(s) in `results/`, and the per-experiment doc in
`docs/`. Multi-panel figures list the panel→file mapping.

**Numbered figure set:** `figures/numbered/Fig01…Fig10` and `FigS1…FigS3` (PNG + PDF) are the
manuscript-ordered figures, one file per figure number (multi-panel figures are montaged vertically,
one source panel per row). Regenerate with `python scripts/number_figures.py` after any panel figure
changes.

## Main figures

**Figure 1. Strain2bScan overview and the two input modes.** `figures/overview.{png,pdf}`
Pipeline schematic. (A) Reference construction: type-IIB (2bRAD) digestion of reference genomes into
single-copy 32–38 bp tags → within-species clustering at 0.95 Jaccard (MinHash-accelerated) → a
cluster × marker database of species-core, cluster-specific and strain-specific markers. (B) Profiling:
a sample is digested once into canonical markers, gated on species-specific markers (Layer-1), and
strains are detected and quantified within each present species from unique markers (Layer-2). The two
input modes — in-silico digestion of conventional shotgun, or native BcgI 2bRAD libraries whose reads
are the tags — enter the same tag space, shared with the Fast2bRAD-M species layer.

**Figure 2. 2bRAD tags track genome-wide strain distance; 16S does not.** `figures/mash_2brad_vs_16s.{png,pdf}`
Per-species Spearman correlation between whole-genome (bottom-3000 21-mer MinHash) between-strain
distance and the corresponding 2bRAD-tag (blue) or 16S rRNA (red) distance, across 15 species
(complete/near-complete genomes only, CheckM ≥97 %/≤5 %; 14–50 genomes per species, shown at right).
Error bars are 95 % CIs over genome subsamples. Genus and species on separate label lines. 2bRAD median
0.94 vs 16S median 0.36. Source: `results/mash_2brad_vs_16s_recomputed.tsv`, `data/genome_qc_16s_panel.tsv`;
doc: `docs/motivation_16s.md`.

**Figure 3. Accurate strain profiling and depth sensitivity.** (A) `figures/cross_species.{png,pdf}`
— precision/recall/abundance error (Bray–Curtis) on real reference panels with simulated mixtures for
*C. acnes*, *S. aureus*, *S. epidermidis* (precision 1.0 throughout). (B) `figures/depth_sensitivity.{png,pdf}`
— detection of a single *C. acnes* strain across a coverage ladder; onset at 0.5×, matching StrainScan.
Sources: `results/cross_species.tsv`, `results/depth_sensitivity.tsv`; docs: `docs/cross_species.md`,
`docs/depth_sensitivity.md`.

**Figure 4. Robustness to reference panel size and reference quality.** (A) `figures/panelsize_prevotella.{png,pdf}`
— *P. copri* panels of 40/80/112 genomes: precision 1.0 at every size; 112→51 clusters matches
StrainScan's own DB. (B) `figures/refqual_figure.{png,pdf}` — reference-degradation gradient: precision
1.0 to 70 % completeness, collapse at 50 %/10 %/400-contigs. Sources: `results/panelsize_prevotella.tsv`,
`results/refqual_degradation.tsv`; doc: `docs/refgenome_quality.md`.

**Figure 5. The 2bRAD enzyme set is a resolution/cost knob enabling native BcgI operation.**
`figures/enzyme_sweep.{png,pdf}` — on *C. acnes*, precision/recall, strain-specific marker count and
build time vs number of enzymes (1→14). BcgI-alone gives precision 1.0; ~4 enzymes is the sweet spot;
cluster count is invariant (16). Single-enzyme BcgI is what enables native 2bRAD-M libraries. Source:
`results/enzyme_sweep.tsv`; doc: `docs/enzyme_sweep.md`.

**Figure 6. DNA mock, host-contamination series: native 2bRAD holds strain recovery where shotgun collapses.**
`figures/mock_hostcontam.{png,pdf}` — ATCC MSA-1002 mock across 0/90/99 % human DNA, native BcgI 2bRAD
vs in-silico-digested shotgun (WMS), same 62-species BcgI DB. (A) Species recall: precision 1.0 for both
at every level; 2bRAD 20/20 at 99 % host vs shotgun 12/20. (B) Usable BcgI marker yield: 2bRAD flat and
high; shotgun 96k→53k→33k as host rises. Source: `results/mock_hostcontam.tsv`; doc: `docs/mock_hostcontam.md`.

**Figure 7. Real saliva: strain profiles discriminate individuals and are temporally stable.**
(A–C) `figures/saliva_individual_discrimination.{png,pdf}` — PCoA (species vs strain) coloured by subject,
and per-species strain-level subject R² (8 subjects × 4 timepoints; strain R² 0.833 > species 0.822,
p 2e-4; *Rothia mucilaginosa* R² 0.921; 17/18 species significant). (D–E) `figures/saliva_temporal_ml.{png,pdf}`
— within- vs between-subject strain distance (0.19 vs 0.44, p 5e-25) and leave-one-timepoint-out host-ID
accuracy (100 % strain vs 94 % species). Sources: `results/saliva_permanova.tsv`,
`results/saliva_perspecies_subject.tsv`, `results/saliva_temporal_ml.tsv`; docs:
`docs/saliva_individual_discrimination.md`, `docs/saliva_temporal_ml.md`.

**Figure 8. Native 2bRAD confirms all shotgun strains and recovers the low-abundance ones shotgun misses.**
`figures/saliva_concordance.{png,pdf}` — paired shotgun↔2bRAD on the same saliva samples. (A) Per sample,
strains detected: shotgun-confirmed (100 % of shotgun calls, 65/65, are in native 2bRAD) plus 2bRAD-only
(128–163/sample). (B) Community relative abundance of shared vs 2bRAD-only strains (2bRAD-only
significantly lower, p 1.2e-23). Source: `results/saliva_concordance.tsv`; doc: `docs/saliva_concordance.md`.

**Figure 9. Fast, light, and scalable to whole communities.** (A) `figures/performance.{png,pdf}` —
per-sample time/memory vs StrainScan (~8×/~11×). (B) `figures/scalability.{png,pdf}` — thread scaling of
build and profile (8.5×/6.5× to 16 threads). (C) `figures/community_throughput.{png,pdf}` — 55-species
community: per-sample cost flat in #species, ~121–146× faster than a per-species tool run once per species.
Sources: `results/headtohead_performance.tsv`, `results/parallel_and_build_scaling.tsv`,
`results/multispecies_*.tsv`, `results/community_throughput.tsv`; docs: `docs/scalability.md`, `docs/multispecies.md`.

**Figure 10. Matches or exceeds StrainScan on its own databases.** `figures/species_expansion.{png,pdf}`
— head-to-head on StrainScan's reference sets: *A. muciniphila*, *P. copri* (precision 1.0 both; recall
0.93 vs 0.24, 0.94 vs 0.90; ~17–23× faster, ~15–24× lighter) and near-clonal *M. tuberculosis*, where
StrainScan did not complete (>3.3 h, >25 GB) and Strain2bScan finished in 0.89 s. Sources:
`results/species_expansion.tsv`, `results/headtohead_performance.tsv`; docs: `docs/species_expansion.md`,
`docs/headtohead_strainscan.md`.

## Supplementary figures

**Figure S1. Per-species rank–rank scatter of between-strain distance (2bRAD and 16S vs whole-genome).**
`figures/mash_2brad_vs_16s_scatter.{png,pdf}` — pair-level view behind Fig 2; 2bRAD hugs the diagonal in
every species, 16S forms flat rank-bands (most extreme in *P. dorei*, *M. tuberculosis*). Data:
`results/pairdist/*.tsv`.

**Figure S2. ATCC MSA-1002 DNA-input titration (native 2bRAD).** `figures/mock_msa1002_titration.{png,pdf}`
— precision 1.0 and full recall to 0.1 ng input. Source: `results/mock_msa1002_titration.tsv`; doc:
`docs/mock_msa1002.md`.

**Figure S3. Layer-1 gate calibration on the 55-species panel.** `figures/gate_calibration.{png,pdf}` —
species precision/recall vs the marker floor at normal and low depth. Source:
`results/gate_calibration_*.tsv`; doc: `docs/gate_calibration.md`.

**Table S4. Exploratory clinical oral cohort profiling.** `results/clinical_exploratory.tsv` — 4 clinical
oral 2bRAD samples, 15–17 species and 115–158 strain calls each, ~2 s/sample (no case/control labels
public). Doc: `docs/clinical_oral_exploratory.md`.

**Table S5. Genome quality-control for the 16S/2bRAD motivation panel.** `data/genome_qc_16s_panel.tsv` —
assembly level, CheckM completeness/contamination, contig count and length per genome; high-quality flag.
