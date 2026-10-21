# Strain2bScan: strain-level 2bRAD-marker profiling of microbiomes

## Abstract

Strain-level variation underlies clinically important microbial phenotypes, but 16S rRNA surveys resolve
species rather than strains, whereas full-k-mer strain profilers can incur repeated per-species search
costs and lose signal in host-rich samples. We present Strain2bScan, a Rust profiler that converts
type-IIB 2bRAD restriction markers into a strain-cluster-resolved, digest-once search framework for
native BcgI 2bRAD-M libraries and conventional shotgun reads. Across 15 species, 2bRAD markers preserved
genome-wide strain relationships that 16S markers did not (median Spearman ρ = 0.94 versus 0.36). Native
BcgI profiling retained F1 = 1.0 at a 10⁻⁴ abundance threshold in the MSA-1002 99% human-DNA mock, where
the primary shotgun comparators did not retain both detection and abundance. In 225 matched simulations,
Strain2bScan matched StrainScan species-median precision and had higher failure-as-absence recall
(+0.127; 95% CI 0.007–0.248), although multi-species F1 was lower. In 32 gate-matched saliva pairs,
full-depth WMS and native profiles shared 495 strain-cluster calls, with additional calls unique to each
modality. Same-container profiling used 4–33× less wall time and 5–39× less memory than StrainScan; an
archived cross-platform database benchmark showed larger differences. Accuracy depended on reference
completeness and niche-matched panel design. Strain2bScan therefore extends reduced-representation
sequencing from species profiles to scalable strain-cluster resolution.

**Availability.** Rust source: https://github.com/HuangShiLab/Strain2bScan. Public reads: PRJNA1131785,
PRJNA288562, PRJNA1517970, PRJNA1191223 and PRJNA1191225. Derived outputs:
https://github.com/HuangShiLab/Strain2bScan-paper.
