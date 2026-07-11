# Results & figure organization — following StrainScan's logic flow

Reference: Liao et al., *High-resolution strain-level microbiome composition analysis from short
reads*, **Microbiome** 11:183 (2023) — the StrainScan paper (DOI 10.1186/s40168-023-01615-w).

## StrainScan's figure logic flow (the template)

StrainScan organizes 11 figures as **concept → method → simulated (escalating difficulty, vs
many tools) → real-data case studies**:

| Fig | Role | Content |
|---|---|---|
| 1 | **Motivation** | why strain-level > cluster-level (vs StrainGE/StrainEst); gene-content differences *within* a cluster |
| 2 | **Overview/schematic** | the pipeline (genome clustering → CST) |
| 3–5 | **Method detail** | k-mer assignment, cluster search, within-cluster strain ID |
| 6 | **Simulated accuracy #1 — single strain × depth** | F1 of **7 tools** across sequencing depths + **running time** |
| 7 | **Simulated accuracy #2 — co-existing strains** | F1 of 5 tools on 2/3/5-strain mixtures of varying similarity + running time |
| 8 | **Simulated abundance accuracy** | Jensen–Shannon divergence (truth vs predicted abundance), 5 tools |
| 9 | **Robustness — novel/distant strains** | Mash distance vs ground truth, binned by how close the true strain is to the DB (simulated + real *E. coli*) |
| 10 | **Real case study #1** | *S. epidermidis* antibiotic-resistance **dynamics** (real metagenomes, ± antibiotic; abundance by 5 tools + gene-content validation) |
| 11 | **Real case study #2** | *C. acnes* skin **diversity** (real metagenomes, sites × individuals × timepoints + phylogenetic tree) |

Plus text/supp case studies: pathogen detection, low-depth *C. difficile*.

**The takeaways of the template:**
1. Lead with **why** (motivation) + **how** (schematic), then benchmark.
2. Simulated benchmarks **escalate in difficulty**: single → co-existing → abundance → open-world
   (novel strain), each with **multi-tool comparison** and a **running-time** panel.
3. Close with **real-data case studies**, each a biological story with an **orthogonal validation**
   (gene content, phylogeny, a known resistant/pathogenic strain, temporal dynamics).

## Proposed Strain2bScan Results structure (parallel flow, our differentiators foregrounded)

Our unique angles vs StrainScan: (i) **efficiency + community-scale** (StrainScan scales poorly —
our strength), (ii) **2bRAD-native + Fast2bRAD-M interoperability**, (iii) the **enzyme knob**.
So we keep the flow but foreground efficiency and add the 2bRAD story.

| Fig | Role (StrainScan analog) | Our content | Status |
|---|---|---|---|
| **1** | Concept + overview (their 1+2) | Schematic: 2bRAD-reduced markers + clustering + strain ID; the **two data modes** (shotgun in-silico digest / real BcgI 2bRAD); ~50–100× sparser than full k-mers; **tags interoperable with Fast2bRAD-M** | **TO MAKE** (no schematic yet) |
| **2** | Simulated accuracy, co-existing strains (their 7+8) | Cross-species strain mixtures: **precision / recall / abundance (1−Bray–Curtis)** for the 5–6 species; all precision 1.0 | ✅ `cross_species` (extend to add species-expansion species) |
| **3** | Accuracy vs depth (their 6A) | Detection vs 0.1–5× depth, Strain2bScan vs StrainScan (matches at 0.5×) | ✅ `depth_sensitivity` |
| **4** | (Strain2bScan-specific) | **The 2bRAD enzyme knob**: BcgI-alone works, ~4-enzyme sweet spot, markers/cost vs #enzymes | ✅ `enzyme_sweep` |
| **5** | Robustness (their 9) | Reference **panel size** (P. copri 40/80/112 → clusters match StrainScan) + reference **quality** degradation | ✅ `panelsize_prevotella` + `refqual_figure` (combine as A/B) |
| **6** | **Efficiency & scalability (their 6B/7B, expanded — our headline)** | Per-sample speed/memory (~8×/~11×) + parallelism + **community-scale ~130–146× on 55 species × hundreds of samples** (flat in species, linear in samples) | ✅ `performance` + `scalability` + `community_throughput` (multi-panel) |
| **7** | Head-to-head vs StrainScan (their tool comparison) | 3 species on **StrainScan's own DBs**: matched precision, ~17–23× faster, **StrainScan DNF on near-clonal *M. tuberculosis*** | ✅ `species_expansion` |
| **8** | **Real-data case study (their 10+11)** | **Saliva**: strain-level individual discrimination (PERMANOVA R² **>** species-level) + **paired shotgun↔2bRAD concordance** + temporal stability; reproduces strain2bfunc at ~78 MB/~1 s | ✅ **DONE** — `docs/saliva_individual_discrimination.md`; subject strain R²=0.833 > species 0.822 (p=2e-4), 1-NN 90.6%, timepoint ns, *Rothia mucilaginosa* R²=0.921. Shotgun↔2bRAD concordance still open. |

Additional real-data figure completed since this table:
- **Fig 9 mock host-contamination** (`docs/mock_hostcontam.md`): native 2bRAD vs shotgun on the
  ATCC MSA-1002 mock across 0/90/99% human DNA — precision 1.0 throughout for both; 2bRAD holds
  20/20 recall at 99% host vs shotgun 12/20 (2bRAD preserves ~10× more usable markers).

Supplementary: parallel-speedup + MinHash-vs-exact validation, gate-sweep accuracy table,
CheckM2 validation of degraded refs, per-sample tables, BcgI-2bRAD verification (`benchmark #9`).

## Gaps to fill before this reads as a complete Microbiome-style paper

1. **Fig 1 schematic** — we have no overview/pipeline figure (StrainScan's Fig 2 is central). Make
   one: genomes → 2bRAD digest (16 enzymes) → single-copy tags → cluster (0.95 Jaccard, MinHash) →
   marker classes → sample digest (shotgun **or** BcgI 2bRAD) → detect + quantify. Highlight the
   two data modes and Fast2bRAD-M interoperability.
2. **Real-data case study (Fig 8)** — the biggest gap; the saliva paired shotgun+2bRAD dataset fills
   it and gives the biological headline (individual discrimination) + self-validation.
3. **Multi-tool comparison (optional)** — StrainScan benchmarks 5–7 tools. We compare mainly vs
   StrainScan. Either add ≥1–2 more (StrainGE/StrainEst/sylph) or explicitly position our claim as
   **efficiency + 2bRAD-native**, not "most accurate of all tools." Recommend adding sylph and
   StrainGE at least on the cross-species set for credibility.
4. **Abundance metric alignment** — StrainScan reports JSD; we report Bray–Curtis. Consider adding
   JSD (or L1) for direct comparability, or keep Bray–Curtis and note the mapping.

## One-line summary of the flow to adopt
**Concept/schematic → simulated accuracy (co-existing strains, across species) → depth → the
enzyme knob → robustness (panel size + quality) → efficiency & community scale (our headline) →
head-to-head vs StrainScan → real saliva case study (individual discrimination + paired
shotgun/2bRAD).**

## Additions from prior Strain2bfunc data (PPTX 2026-07-10)

Three prior results from the Strain2bfunc project (`Strain2bFunc_20260710--Shi.pptx`) strengthen
the paper. All were generated with **Strain2bfunc** (the earlier heavy Python tool). The
2bRAD-representation result is tool-agnostic and reusable as-is; the two data-driven results should
be **re-run with Strain2bScan** to show it reproduces them at a fraction of the compute.

**New motivation figure — 2bRAD captures strain-level signal that 16S cannot. ✅ DONE (all 15
species, recomputed with Strain2bScan).** `docs/motivation_16s.md`,
`results/mash_2brad_vs_16s_recomputed.tsv`, `figures/mash_2brad_vs_16s.*`. Between-strain distance
vs whole-genome (Spearman): **2bRAD median 0.90 (0.59–0.99) vs 16S median 0.36 (−0.14 to 0.78)** —
e.g. *M. tuberculosis* 0.90 vs −0.14, *S. aureus* 0.98 vs 0.09, *L. plantarum* 0.87 vs 0.02. WGS =
bottom-3000 21-mer MinHash; 2bRAD = Strain2bScan BcgI tags; 16S = longest barrnap gene, 21-mer
Jaccard. Reproduces the slide-deck finding. This is the "why 2bRAD for strains" motivation
(StrainScan's Fig 1 analog).

**New figure — DNA mock community, real DNA with known truth (low-biomass niche). ✅ DONE (2bRAD
side) with Strain2bScan.** Raw ATCC MSA-1002 native BcgI 2bRAD reads obtained from Figshare
(12272360) and profiled with Strain2bScan against a 62-species BcgI panel (20 mock + 42 decoys):
**precision 1.0 at every DNA input, recall 20/20 at ≥0.1 ng and 19/20 at 0.01 ng** — reproducing the
prior Strain2bfunc result on real, error-containing data and demonstrating the low-biomass niche.
See `docs/mock_msa1002.md`, `results/mock_msa1002_titration.tsv`, `figures/mock_msa1002_titration.*`.
Still open: cross-species abundance even-ness; the WMS/inStrain contrast (needs the 2.6 GB shotgun +
inStrain); the 90/95/99 % host-contamination series (not in this archive); MSA-1003. The original
PPTX numbers for reference — Strain2bfunc-2bRAD precision/recall/F1 = 1.0 down to 0.001 ng & 99 %
human DNA vs inStrain-WMS precision 0.4–0.5.

**Saliva case study — now with a target result.** The saliva chapter (8 subjects × 4 timepoints =
32) already has a Strain2bfunc result to reproduce: strain-level Adonis R² up to **0.87** (*Rothia
mucilaginosa*) > **species-level 0.72**, and a host-ID ML model at **100% accuracy from 65
strain-level features (6 species)** vs 96.9% from 1530 species-level features (slides 31–32).

**Framing decision (needed).** These are Strain2bfunc outputs. To carry the tool-specific ones
(mock, saliva) into a *Strain2bScan* paper, either (a) **re-run with Strain2bScan** on the raw reads
(best — shows reproduction at lower compute; needs the reads on the Mac Studio/HPC), or (b) present
as prior validation of the 2bRAD approach with attribution. The Mash-distance figure works directly
either way.

**Proposed updated main-figure sequence:** Fig 1 schematic → **Fig 2 2bRAD-vs-16S motivation** →
Fig 3 cross-species accuracy → Fig 4 depth → Fig 5 enzyme knob → Fig 6 robustness → Fig 7 efficiency
& community scale → Fig 8 head-to-head vs StrainScan → **Fig 9 DNA mock community (real DNA, vs
inStrain, low-biomass/high-host)** → Fig 10 saliva case study.

The low-biomass / high-host-contamination angle (slides 10–11: saliva ≥90% human DNA, FFPE <150 bp,
skin <1 ng) is a distinct 2bRAD selling point worth foregrounding in the Introduction, and the mock
figure is its quantitative payoff.
