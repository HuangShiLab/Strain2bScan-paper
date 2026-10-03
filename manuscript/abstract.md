# Strain2bScan; Abstract (draft)

**Working title:** *Strain2bScan: strain-level profiling from 2bRAD-reduced markers for low-biomass
microbiomes and community-scale metagenomes.*

## Abstract

**Background.** Strain-level variation drives clinically important microbial phenotypes. Although 16S
surveys resolve species, they generally cannot resolve strains. Shotgun k-mer methods recover strain
information but face two practical barriers: per-sample cost grows when many species must be queried, and
host-dominated low-biomass specimens contain too little microbial DNA for robust detection. 2bRAD sequencing
samples a sparse, reproducible fraction of the genome, is compatible with low-input and host-contaminated
specimens, and has previously been used mainly for species-level profiling.

**Results.** We present Strain2bScan, a Rust strain profiler that implements within-species clustering and
unique-marker scoring on 2bRAD markers, the 25–33 bp tags produced by type-IIB restriction digestion. It
accepts both native 2bRAD-M libraries and in-silico-digested shotgun metagenomes. 2bRAD-marker distances
tracked whole-genome strain distances (median Spearman ρ = 0.94), whereas 16S distances were much weaker
(median ρ = 0.36). With complete references, precision was 1.0 across simulated single-species benchmarks,
detection reached 0.5× coverage, and accuracy depended strongly on reference completeness. In native BcgI
2bRAD mode, Strain2bScan retained precision 1.0, full strain recall and accurate abundance at 99% human DNA.
On real saliva, leave-one-timepoint-out host identification reached 100% with strain features versus 78.1%
with species features; strain and species PERMANOVA R² values were comparable (0.757 versus 0.755), so the
advantage was strongest for subject identification and low-abundance strain recovery. Native 2bRAD confirmed
65/65 shotgun strain calls and recovered 128–163 additional low-abundance strains per sample. For
conventional shotgun input, digest-once profiling reduced cost relative to per-species k-mer querying:
Strain2bScan was approximately 8× faster and 11× lighter per sample and 130–146× faster on a 55-species
community. On a matched 15-species benchmark it matched StrainScan precision (1.0), improved median recall
(0.80 versus 0.67), built databases 249–614× faster and 43–138× lighter, and completed *Klebsiella
pneumoniae*, which StrainScan could not build. On host-contaminated shotgun mocks, it was the only tool
tested that preserved both detection and abundance at 99% human DNA. In public WGS cohorts, 22 libraries
were profiled in 100.6 s using 343 MiB peak RSS; a cohort-specific isolate panel recovered all six expected
isolate self-calls and tracked a persistent *Bifidobacterium bifidum* strain across three weekly infant
stool samples.

**Conclusions.** Strain2bScan makes genome-resolved strain profiling practical for both native
low-biomass/high-host 2bRAD-M data and community-scale shotgun cohorts. Its accuracy depends on reference
completeness and panel design, so niche-appropriate or cohort-specific panels are essential for
biological interpretation.

**Availability.** Rust source: https://github.com/HuangShiLab/Strain2bScan. Public reads: PRJNA1131785,
PRJNA288562, PRJNA1517970, PRJNA1191223 and PRJNA1191225. Derived outputs and reproducibility artifacts:
https://github.com/HuangShiLab/Strain2bScan-paper.

---

### Author-facing notes (delete before submission)
- **Framing = two input modes**: (1) native 2bRAD-M → low-biomass/high-host strain analysis;
  (2) in-silico-digested shotgun → community-scale throughput vs StrainScan. Shared foundation
  (accuracy, robustness, 2bRAD-vs-16S motivation) precedes the two pillars. Full section→evidence map in
  `docs/manuscript_organization.md`.
- Figure map (12 main): 1 overview, 2 mash_2brad_vs_16s (A bars + B 3×5 rank–rank scatter, combined),
  3 cross_species(+depth), 4 reference-genome completeness (refqual, all 15 species), 5 enzyme_sweep,
  6 native-2bRAD strain-ID + abundance across four ATCC mocks (fig6_2brad),
  7 saliva (individual+temporal_ml), 8 saliva_concordance, 9 efficiency (performance+scalability+
  community_throughput), 10 species_expansion, 11 sim_headtohead (systematic 15-species head-to-head:
  accuracy + build/profile/multi cost), 12 shotgun strain-ID vs StrainScan/inStrain, incl. high host
  contamination (fig12_wms_toolcompare). Supp: S-tree DB-expansion cost (20 vs 28 species),
  S-inStrain dereplication control, S1 mock_msa1002_titration, S2 gate_calibration,
  Table 1–3 sim head-to-head (`manuscript/tables.md`), Table S3 clinical_exploratory,
  Table S4 genome_qc_16s_panel. (Former Fig S1 scatter is now Fig 2B.)
- All numbers regenerated with the corrected-enzyme binary; 2bRAD-native results on real error-containing
  reads; simulated benchmarks are closed-world (stated in Discussion).
