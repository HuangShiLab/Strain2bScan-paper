# Strain2bScan: strain-level 2bRAD-marker profiling of microbiomes

## Abstract

Strain-level variation underlies clinically important microbial phenotypes, but 16S rRNA surveys resolve
species rather than strains, whereas commonly used full-k-mer strain profilers can incur repeated per-species
search costs and lose signal in host-rich samples. We present Strain2bScan, a Rust profiler that converts
type-IIB 2bRAD restriction markers into a strain-resolved, digest-once search framework for native BcgI
2bRAD-M libraries and conventional shotgun reads. Across 15 species, 2bRAD markers preserved genome-wide
strain relationships that 16S markers did not (median Spearman ρ = 0.94 versus 0.36). Native BcgI profiling
retained F1 = 1.0 at a 10⁻⁴ abundance threshold in a 99% human-DNA mock, where Strain2bScan was the only
tested tool preserving both detection and abundance. In 225 matched simulations, it matched StrainScan
species-median precision and had higher failure-as-absence recall (+0.127; 95% CI 0.007–0.248). Database
construction used 249–614× less wall time and 43–138× less memory than StrainScan. In saliva, strain
features supported within-subject discrimination (100% versus 78.1% leave-one-timepoint-out
classification). Accuracy depended on reference completeness and niche-matched panel design.
Strain2bScan therefore extends reduced-representation sequencing from species profiles to scalable strain
resolution.

**Availability.** Rust source: https://github.com/HuangShiLab/Strain2bScan. Public reads: PRJNA1131785,
PRJNA288562, PRJNA1517970, PRJNA1191223 and PRJNA1191225. Derived outputs:
https://github.com/HuangShiLab/Strain2bScan-paper.
