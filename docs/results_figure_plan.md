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
| **8** | **Real-data case study (their 10+11)** | **Saliva**: strain-level individual discrimination (PERMANOVA R² **>** species-level) + **paired shotgun↔2bRAD concordance** + temporal stability; reproduces strain2bfunc at ~78 MB/~1 s | **TO DO** (needs the saliva data; see `docs/real_metagenome_plan.md`) |

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
