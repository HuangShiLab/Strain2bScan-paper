# Strain2bScan — Abstract (draft)

**Working title:** *Strain2bScan: strain-level profiling from 2bRAD-reduced markers for low-biomass
microbiomes and community-scale metagenomes.*

## Abstract

**Background.** Strain-level variation drives clinically important microbial phenotypes, but 16S
surveys cannot resolve strains and shotgun k-mer methods, while accurate, are held back from two
increasingly important settings: they are computationally heavy across the many-species, many-sample
cohorts of modern studies, and they fail on low-biomass, host-contaminated samples (saliva, tumour/FFPE,
skin) where most DNA is human. Reduced-representation 2bRAD sequencing samples a sparse, reproducible
subset of the genome, is robust to low input and host contamination, and has been used for species-level
profiling — but not for strain resolution.

**Results.** We present **Strain2bScan**, a dependency-free Rust strain profiler that ports the
StrainScan resolution framework (within-species clustering plus unique-marker scoring) onto **2bRAD
markers** — the 32–38 bp tags released by type-IIB restriction digestion — and uniquely accepts **two
input modes**. As a shared foundation it is accurate (**precision 1.0** across species, high recall, low
abundance error), detects strains to **0.5× coverage matching StrainScan**, is robust to reference panel
size and to reference degradation to 70 % completeness, and — unlike 16S, whose between-strain distances
are uncorrelated with genome divergence (median Spearman 0.36 across 15 species) — its tags track
whole-genome strain distance (median 0.94). **(1) Native 2bRAD-M for low-biomass microbiomes:** on the
ATCC MSA-1002 mock, native 2bRAD holds precision 1.0 and **20/20 strain recall at 99 % human DNA** where
in-silico-digested shotgun drops to 12/20; on real saliva (8 subjects × 4 timepoints) strain-level
profiles discriminate individuals better than species-level (PERMANOVA **R² 0.83 > 0.82**; leave-one-
timepoint-out host-ID **100 %**), are temporally stable, and recover **128–163 low-abundance strains per
sample that host-limited shotgun misses**. **(2) Conventional metagenomes at scale:** the sample is
digested once and matched against every species database, so per-sample cost is independent of the
number of species — **~8× faster and ~11× lighter** than StrainScan per sample and **~130–146× faster on
a 55-species community**, while matching or exceeding StrainScan's recall on its own databases and
completing near-clonal *Mycobacterium tuberculosis* in ~1 s where StrainScan does not (>3 h, >25 GB). The
two modes agree, so one tool spans both regimes.

**Conclusions.** Strain2bScan delivers accurate, genome-resolved strain profiling from a sparse marker
set: it uniquely enables strain-level analysis of low-biomass, high-host 2bRAD-M data, and scales
conventional-metagenome strain profiling to communities of many species across many samples.

**Availability.** Rust source: https://github.com/HuangShiLab/Strain2bScan · reproducible benchmarks and
figures: https://github.com/HuangShiLab/Strain2bScan-paper.

---

### Author-facing notes (delete before submission)
- **Framing = two input modes**: (1) native 2bRAD-M → low-biomass/high-host strain analysis;
  (2) in-silico-digested shotgun → community-scale throughput vs StrainScan. Shared foundation
  (accuracy, robustness, 2bRAD-vs-16S motivation) precedes the two pillars. Full section→evidence map in
  `docs/manuscript_organization.md`.
- Figure map (10 main): 1 overview, 2 mash_2brad_vs_16s (A bars + B 3×5 rank–rank scatter, combined),
  3 cross_species(+depth), 4 robustness (panelsize+refqual), 5 enzyme_sweep, 6 mock_hostcontam,
  7 saliva (individual+temporal_ml), 8 saliva_concordance, 9 efficiency (performance+scalability+
  community_throughput), 10 species_expansion. Supp: S1 mock_msa1002_titration, S2 gate_calibration,
  Table S3 clinical_exploratory, Table S4 genome_qc_16s_panel. (Former Fig S1 scatter is now Fig 2B.)
- All numbers regenerated with the corrected-enzyme binary; 2bRAD-native results on real error-containing
  reads; simulated benchmarks are closed-world (stated in Discussion).
