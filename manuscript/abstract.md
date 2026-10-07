# Strain2bScan: strain-level profiling from 2bRAD-reduced markers for low-biomass microbiomes and community-scale metagenomes

## Abstract

Strain-level variation underlies clinically important microbial phenotypes, but 16S rRNA surveys usually
resolve species rather than strains. Shotgun k-mer profilers recover strain signal, yet repeated
per-species queries scale poorly across communities and host-dominated low-biomass specimens often contain
too little microbial DNA. 2bRAD sequencing samples a reproducible genomic fraction before host DNA dominates
a library, but has mainly supported species-level profiling. We present Strain2bScan, a Rust profiler that
clusters conspecific reference genomes and scores unique 25–33 bp type-IIB restriction markers. It accepts
native BcgI 2bRAD-M libraries and in-silico-digested shotgun reads in one marker framework. Across 15
species, 2bRAD-marker distances tracked whole-genome strain distances (median Spearman ρ = 0.94), whereas
16S distances were weaker (median ρ = 0.36). With complete references, precision was 1.0 across simulated
single-species benchmarks, detection reached 0.5× coverage, and accuracy depended strongly on reference
completeness. In the primary native BcgI 2bRAD-M run at 99% human DNA, Strain2bScan retained F1 = 1.0 at an
abundance threshold of 10⁻⁴ and had the lowest Bray–Curtis dissimilarity among tested tools. In 8 saliva
subjects, leave-one-timepoint-out host identification reached 100% with strain features versus 78.1% with
species features, although PERMANOVA R² was similar (0.757 versus 0.755). In three usable paired samples,
all 65 shotgun strain-cluster calls were also present in native BcgI 2bRAD-M, which yielded 128–163
additional candidate strain-cluster calls per sample. In single-panel testing, Strain2bScan was
approximately 8× faster and 11× lighter than StrainScan; a projected 55-species comparison was 121–146×
lower in cost than per-species querying. In a run-level 15-species StrainScan rerun, species-median
precision was 1.0 for both tools, Strain2bScan had higher median recall/F1 (0.733/0.844 versus
0.667/0.800), and the 95% intervals for paired accuracy differences included zero; archived timing runs
showed 249–614× faster and 43–138× lighter database construction. In the primary host-contamination
comparison, Strain2bScan was the only tested tool that preserved both detection and abundance at 99% human
DNA. In
exploratory public WGS subsets, 22 libraries were profiled in 100.6 s using 343 MiB peak RSS; a
cohort-specific isolate panel showed the expected 6/6 panel/read compatibility result and detected a
persistent *Bifidobacterium bifidum* strain-resolved cluster across three weekly infant stool samples.
These were feasibility and concordance analyses, not independent strain-accuracy validations. Accuracy
therefore depends on reference completeness and niche-matched panel design.

**Availability.** Rust source: https://github.com/HuangShiLab/Strain2bScan. Public reads: PRJNA1131785,
PRJNA288562, PRJNA1517970, PRJNA1191223 and PRJNA1191225. Derived outputs:
https://github.com/HuangShiLab/Strain2bScan-paper.
