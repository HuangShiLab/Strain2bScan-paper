# Strain2bScan — Abstract (draft)

**Working title:** *Strain2bScan: fast, scalable strain-level metagenomic profiling on
2bRAD-reduced markers.*

## Abstract

**Background.** Strain-level microbiome profiling resolves clinically and ecologically
important within-species variation, but existing k-mer methods are computationally heavy at
the scale of modern cohorts, where hundreds of species must be resolved across hundreds of
samples. Reduced-representation 2bRAD sequencing samples a sparse, reproducible subset of the
genome and has been used for species-level profiling, but not for strain resolution.

**Results.** We present Strain2bScan, a Rust strain profiler that ports the StrainScan
resolution framework (within-species clustering plus unique-marker scoring) onto **2bRAD
markers** — the 32–38 bp tags released by type-IIB restriction digestion, a ~1–2% subset of
the genome. Strain2bScan accepts both BcgI 2bRAD experimental libraries and in-silico
digestion of conventional metagenomes with up to 16 enzymes. On a real *Cutibacterium acnes*
benchmark it profiles **~14× faster and with ~7× less peak memory** than StrainScan at
comparable accuracy for sufficiently distinct, sufficiently covered strains; because the
sample is digested once and matched against every per-species database, its per-sample cost is
**independent of the number of species and linear in the number of samples**, and MinHash-based
clustering scales to large panels without changing the result. We map the operating envelope:
accuracy is high at ≥~5× per-strain depth and for species with high intra-species diversity
(*C. acnes*), and declines below ~1× depth and for low-diversity species whose strains are
near-identical (*Staphylococcus epidermidis*); the number of enzymes is a tunable
efficiency-versus-resolution knob (≈8 enzymes optimal). Strain2bScan reports honestly when a
species is not strain-resolvable and includes assembly-quality filtering that counters the
completeness-driven bias of Jaccard clustering.

**Conclusions.** Strain2bScan trades a controlled amount of low-depth and low-diversity
sensitivity for large gains in speed, memory, and scalability, making genome-resolved
strain surveillance across many species and samples practical, and uniquely enabling
strain-level analysis of BcgI 2bRAD experimental data.

**Availability.** Rust source: https://github.com/HuangShiLab/Strain2bScan · reproducible
benchmarks and figures: https://github.com/HuangShiLab/Strain2bScan-paper.

---

### Author-facing notes (delete before submission)
- Keep the honest framing: **efficiency + envelope + 2bRAD capability**, not "more accurate."
- Numbers to finalize against the last re-verified run: 14× / 7× (Fig 1); ≈8-enzyme optimum
  (Fig 4); resolvability→precision (Fig 5); depth onset 0.5× vs 1×/5× (Fig 3).
- Add one sentence on the StrainScan same-panel head-to-head once the Linux run is in (step 3).
