# Strain2bScan manuscript — re-organized around the two input modes

**Thesis.** Strain2bScan is one strain-profiling engine (2bRAD-reduced markers + StrainScan-style
clustering/unique-marker scoring) that serves **two input modes**, each a distinct contribution:

1. **Native 2bRAD-M libraries → strain-level analysis of low-biomass / high-host microbiomes.**
   The wet-lab 2bRAD reduction enriches informative markers *before* host DNA swamps the library, so
   Strain2bScan resolves strains where shotgun cannot. New algorithm + software make 2bRAD, until now
   a species-level method, strain-resolved.
2. **Conventional shotgun metagenomes (in-silico digestion) → high-throughput, community-scale strain
   profiling.** The sample is digested once and matched against every species DB, so cost is flat in
   #species and linear in #samples — enabling cohort-scale analysis where per-species k-mer tools
   (StrainScan) do not scale.

A shared foundation (concept, why-2bRAD-not-16S, core accuracy) precedes the two pillars.

---

## Suggested title / abstract thesis
> *Strain2bScan: strain-level profiling from 2bRAD-reduced markers for low-biomass microbiomes and
> community-scale metagenomes.*

Abstract arc: strain variation matters → 16S can't resolve it and shotgun is heavy at scale and fails
on low-biomass/high-host samples → 2bRAD reduced-representation solves low-biomass but was species-only
→ **Strain2bScan** ports strain resolution onto 2bRAD tags and accepts **both** native 2bRAD and
in-silico-digested shotgun → (1) precision 1.0, robust, depth-matched to StrainScan; (2) low-biomass
pillar: 20/20 recall at 99 % host vs shotgun 12/20; real saliva strain-level individual discrimination
(R² 0.83 > species 0.82; 100 % host-ID) and recovery of low-abundance strains shotgun misses;
(3) scale pillar: ~8×/~11× per sample, ~130–146× at community scale, completes near-clonal *M. tb*
where StrainScan does not.

---

## Introduction — beats
1. Strain-level variation is clinically/ecologically decisive (AMR, virulence, individual-specific
   strains); species-level is insufficient.
2. **16S resolves species, not strains** — quantitative motivation (Fig 2).
3. Shotgun + k-mer strain tools (StrainScan/StrainGE/inStrain) are (a) **heavy at cohort scale** and
   (b) **fail on low-biomass / high-host material** (saliva, FFPE, skin), where most DNA is human.
4. **2bRAD reduced representation** samples a sparse reproducible genome subset at ~1 % cost, robust to
   low input / host contamination / degradation — but has only ever been *species*-level (2bRAD-M /
   Fast2bRAD-M).
5. **Strain2bScan** ports strain resolution onto 2bRAD tags, tag-compatible with Fast2bRAD-M, and
   uniquely accepts **both** native 2bRAD libraries and in-silico-digested shotgun — bridging the
   low-biomass and the large-scale regimes. State the two contributions explicitly.

---

## Results — Part 0: Concept, motivation, and core accuracy (shared foundation)

### §1  Strain2bScan overview and the two data modes — **Fig 1** (`overview.png`)
Claim: one engine, two inputs; digest → cluster (0.95 Jaccard, MinHash) → cluster×marker DB →
two-layer detect/quantify; tags identical to Fast2bRAD-M. Evidence: schematic. Status ✅.

### §2  2bRAD captures strain signal that 16S cannot — **Fig 2** (`mash_2brad_vs_16s.png`; scatter → **Supp S1** `mash_2brad_vs_16s_scatter.png`)
Claim: across 15 species (complete/near-complete genomes, 95 % CIs), between-strain distance tracks
whole-genome distance far better for 2bRAD than 16S — **median Spearman 0.94 vs 0.36**; several 16S CIs
overlap zero (*M. tuberculosis* −0.04, *L. plantarum* 0.02, *P. dorei* −0.13). Evidence:
`results/mash_2brad_vs_16s_recomputed.tsv`, `results/pairdist/*`, `docs/motivation_16s.md`. Status ✅.

### §3  Accurate, robust strain profiling — **Fig 3** (accuracy) + **Fig 4** (robustness)
The algorithm is validated on real reference panels + simulated mixtures (applies to *both* inputs):
- **Cross-species precision 1.0** with high recall / low abundance error — `cross_species.png`,
  `results/cross_species.tsv`, `docs/cross_species.md`. ✅
- **Detection to 0.5× coverage = StrainScan** — `depth_sensitivity.png`, `docs/depth_sensitivity.md`. ✅
- **Robust to reference panel size** (P. copri 40/80/112 → clusters match StrainScan) and **reference
  quality** (precision 1.0 to 70 % completeness) — `panelsize_prevotella.png` + `refqual_figure.png`,
  `docs/refgenome_quality.md`. ✅
  (Combine accuracy panels as Fig 3, robustness panels as Fig 4, or merge into one multi-panel figure.)

---

## Results — Part I: Native 2bRAD-M for low-biomass / high-host microbiomes (**Input mode 1**)

### §4  A tunable enzyme set enables native BcgI 2bRAD — **Fig 5** (`enzyme_sweep.png`)
Claim: the type-IIB enzyme set is a resolution/cost knob; **BcgI alone** already gives precision 1.0,
~4 enzymes is the sweet spot. Single-enzyme BcgI operation is precisely what lets **native BcgI 2bRAD
libraries** be profiled directly (the bridge into this pillar). Evidence: `enzyme_sweep.tsv`,
`docs/enzyme_sweep.md`. ✅

### §5  DNA mock: 2bRAD holds strain recovery under host contamination where shotgun collapses — **Fig 6** (`mock_hostcontam.png`; titration → **Supp S2** `mock_msa1002_titration.png`)
Claim (the low-biomass payoff, real DNA, known truth):
- ATCC MSA-1002 across **0/90/99 % human DNA**, native 2bRAD vs in-silico-digested shotgun: **precision
  1.0 for both at every level**; native 2bRAD keeps **20/20 recall at 99 % host** vs shotgun **12/20**,
  because 2bRAD preserves ~10× more usable markers under host load. `results/mock_hostcontam.tsv`,
  `docs/mock_hostcontam.md`. ✅
- MSA-1002 DNA-input titration: precision 1.0, full recall to **0.1 ng**. `mock_msa1002_titration.tsv`,
  `docs/mock_msa1002.md`. ✅

### §6  Real saliva case study: individual-specific, temporally stable strain signatures — **Fig 7** (multi-panel: `saliva_individual_discrimination.png`, `saliva_temporal_ml.png`) + **Fig 8** concordance (`saliva_concordance.png`)
Native BcgI 2bRAD saliva, 8 subjects × 4 timepoints, oral-commensal panel. The biological headline:
- **Individual discrimination**: subject strain-level PERMANOVA **R² 0.833 > species 0.822**
  (p 2e-4), LOO 1-NN 90.6 %; per-species ***Rothia mucilaginosa* R² 0.921**, 17/18 species significant.
  `saliva_permanova.tsv`, `saliva_perspecies_subject.tsv`, `docs/saliva_individual_discrimination.md`. ✅
- **Temporal stability + host-ID**: within-subject distance 0.19 ≪ between 0.44 (p 5e-25);
  leave-one-timepoint-out subject ID **100 % strain vs 94 % species**. `saliva_temporal_ml.tsv`,
  `docs/saliva_temporal_ml.md`. ✅
- **Paired shotgun↔2bRAD concordance** (self-validation + 2bRAD sensitivity): 100 % of shotgun strain
  calls confirmed by native 2bRAD, and native 2bRAD additionally recovers **128–163 low-abundance
  strains/sample that host-limited shotgun misses** (median rel-ab 0.003 vs 0.010, p 1.2e-23).
  `saliva_concordance.tsv`, `docs/saliva_concordance.md`. ✅
- Clinical oral cohort profiles cleanly (15–17 species, ~2 s/sample) — **Supp / brief mention**; no
  case/control labels public. `clinical_exploratory.tsv`, `docs/clinical_oral_exploratory.md`. ✅

---

## Results — Part II: Conventional metagenomes at community scale (**Input mode 2**)

### §7  Fast, light, and scalable to whole communities — **Fig 9** (`performance.png` + `scalability.png` + `community_throughput.png`)
Claim (the throughput payoff on shotgun input): per sample **~8× faster / ~11× lighter** than
StrainScan; both build and profile parallelise (8.5×/6.5× to 16 threads); the sample is digested
**once**, so per-sample cost is **flat in #species** — on a **55-species** community **~121–146×**
faster than a per-species tool run once per species (~132× at 100 samples). Evidence:
`headtohead_performance.tsv`, `parallel_and_build_scaling.tsv`, `multispecies_*.tsv`,
`community_throughput.tsv`, `docs/scalability.md`, `docs/multispecies.md`. ✅

### §8  Matches or exceeds StrainScan on its own databases — **Fig 10** (`species_expansion.png`)
Claim (head-to-head on shotgun input): on StrainScan's own reference sets both reach **precision 1.0**;
Strain2bScan **matches/exceeds recall** (0.93 vs 0.24; 0.94 vs 0.90) at **~17–23× faster / ~15–24×
lighter**, and **completes near-clonal *M. tuberculosis*** (~1 s) where **StrainScan does not** (>3.3 h,
>25 GB). Evidence: `species_expansion.tsv`, `headtohead_performance.tsv`,
`docs/species_expansion.md`, `docs/headtohead_strainscan.md`. ✅

The validity of the in-silico-digestion (shotgun) mode is established independently in §5 (mock shotgun
recovers 20/20 at 0 % host) and §6 (shotgun calls ⊆ native-2bRAD calls) — i.e. the two modes agree.

---

## Discussion — beats
- Two-mode positioning: 2bRAD-native for **low-biomass/high-host** (saliva, FFPE, skin, ancient/degraded);
  in-silico shotgun for **cohort-scale** where k-mer tools don't scale. One tool, one tag space.
- Interoperability with Fast2bRAD-M species layer (shared tags) → species+strain in one pipeline.
- Limitations: recall bounded by **genuine near-clonality** (*M. tb*); host-limited shotgun can't reach
  the low-abundance tail (which is exactly the 2bRAD advantage); strain resolution needs a
  **niche-appropriate, genome-rich reference panel** (the generic-vs-oral-panel lesson);
  all non-real benchmarks simulated/closed-world.
- Future: oral-cancer case/control (needs labels), FFPE/degraded, deeper multi-tool comparison
  (sylph/StrainGE), 16S species-layer half.

---

## Main-figure sequence (10) and supplement
| # | Section | File(s) | Pillar |
|---|---|---|---|
| 1 | Overview / two modes | overview | shared |
| 2 | 2bRAD vs 16S motivation | mash_2brad_vs_16s | shared |
| 3 | Cross-species accuracy (+depth) | cross_species (+depth_sensitivity) | shared |
| 4 | Robustness (panel size + ref quality) | panelsize_prevotella + refqual_figure | shared |
| 5 | Enzyme knob → native BcgI | enzyme_sweep | **I (2bRAD)** |
| 6 | DNA mock, host-contamination | mock_hostcontam | **I (2bRAD)** |
| 7 | Saliva: discrimination + temporal/ML | saliva_individual_discrimination + saliva_temporal_ml | **I (2bRAD)** |
| 8 | Saliva: shotgun↔2bRAD concordance | saliva_concordance | **I (2bRAD)** |
| 9 | Efficiency & community scale | performance + scalability + community_throughput | **II (metagenome)** |
| 10 | Head-to-head vs StrainScan | species_expansion | **II (metagenome)** |

**Supplement:** S1 rank–rank scatter (mash_2brad_vs_16s_scatter); S2 MSA-1002 titration
(mock_msa1002_titration); S3 gate calibration (gate_calibration); S4 clinical exploratory; S5 genome-QC
table (genome_qc_16s_panel.tsv); per-experiment design in `docs/benchmark_datasets.md`.

(Figs 7+8 may merge into one 6-panel saliva figure; Figs 3+4 may merge; that gives an 8-main-figure paper.)
