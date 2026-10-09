# Strain2bScan: strain-level profiling from 2bRAD-reduced markers for low-biomass microbiomes and community-scale metagenomes

## Abstract

Strain-level variation underlies clinically important microbial phenotypes, but 16S rRNA surveys usually
resolve species rather than strains, whereas shotgun profilers can scale poorly and lose signal in
host-dominated, low-biomass samples. We present Strain2bScan, a Rust profiler that clusters conspecific
reference genomes and scores unique 25–33 bp type-IIB restriction markers. It profiles native BcgI
2bRAD-M libraries and in-silico-digested shotgun reads in one marker framework. Across 15 species,
2bRAD-marker distances tracked whole-genome strain distances more closely than 16S distances (median
Spearman ρ = 0.94 versus 0.36). Native BcgI profiling retained F1 = 1.0 at an abundance threshold of 10⁻⁴
in the 99%-human-DNA mock. In the primary shotgun comparison, Strain2bScan was the only tested tool that
preserved detection and abundance at 99% host. In a
run-level StrainScan rerun, accuracy was comparable: species-median precision was 1.0 for both tools,
median recall/F1 were 0.733/0.844 versus 0.667/0.800, and paired-difference intervals included zero.
Archived cross-environment runs reduced database-construction cost by 249–614× in wall time and 43–138×
in memory. A
supplementary StrainGST rerun used a different 0.90-reference space; StrainGR was not run. In eight
saliva subjects, a leave-one-timepoint-out classifier reached 100% with strain features versus 78.1%
with species features; the paired accuracy-difference interval included zero. This was a small
within-subject classification result, not independent strain-accuracy validation. Accuracy depends on
reference completeness and niche-matched panel design.

**Availability.** Rust source: https://github.com/HuangShiLab/Strain2bScan. Public reads: PRJNA1131785,
PRJNA288562, PRJNA1517970, PRJNA1191223 and PRJNA1191225. Derived outputs:
https://github.com/HuangShiLab/Strain2bScan-paper.
