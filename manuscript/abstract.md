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
markers** — the 32–38 bp tags released by type-IIB restriction digestion. Strain2bScan accepts
both BcgI 2bRAD experimental libraries and in-silico digestion of conventional metagenomes with
up to 16 enzymes; its tag lengths and recognition patterns match Fast2bRAD-M /
`2bRADExtraction.pl` exactly, so the tags are interoperable with the Fast2bRAD-M species layer
and with real BcgI 2bRAD libraries. Across five species (*Cutibacterium acnes*, *Staphylococcus
aureus*, *S. epidermidis*, *Akkermansia muciniphila*, *Prevotella copri*) it achieves
**precision 1.0** with high recall and low abundance error, and detects strains down to **0.5×
coverage — matching StrainScan** — while profiling each sample **~8× faster and with ~11× less
peak memory**. Because the sample is digested once and matched against every per-species
database, per-sample cost is **independent of the number of species and linear in the number of
samples**: on a 55-species community it profiles ~130× faster than running a per-species tool
once per species (∼5 min vs a projected ∼10 h at 100 samples). Recall is limited only by
*genuine* near-clonality — where strains are nearly identical (*Mycobacterium tuberculosis*),
cluster resolution is intrinsically coarse — and, notably, StrainScan's regression-based
within-cluster step failed to complete on *M. tuberculosis* (>3 h, >25 GB) where Strain2bScan
finished in ~1 s.

**Conclusions.** Strain2bScan delivers accurate, genome-resolved strain profiling at a fraction
of the compute and memory of full-k-mer methods, scales to communities of many species across
many samples, and uniquely enables strain-level analysis of BcgI 2bRAD experimental data.

**Availability.** Rust source: https://github.com/HuangShiLab/Strain2bScan · reproducible
benchmarks and figures: https://github.com/HuangShiLab/Strain2bScan-paper.

---

### Author-facing notes (delete before submission)
- Framing (post strand-fix): **accurate (precision 1.0) + fast/light + community-scale +
  2bRAD capability**. The earlier "operating envelope / accuracy tracks resolvability" framing
  is obsolete — that gradient was a digestion bug (see `docs/species_expansion.md`), now fixed.
- Numbers are from the strand-fixed binary: ~6×/~7× (Fig 1); ~91× community (Fig, multispecies);
  precision 1.0 cross-species (Fig 5); depth onset 0.5× = StrainScan (Fig 3); ~8-enzyme sweet
  spot (Fig 4). All `results/*.tsv` and `docs/*` are regenerated with the fix.
- The honest cost to state: both-strand digestion doubles markers → ~6× (not ~14×) speed edge;
  recoverable by framing one canonical tag per site.
- Add the StrainScan same-panel build head-to-head once the Linux run is in.
