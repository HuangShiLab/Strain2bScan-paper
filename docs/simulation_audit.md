# Simulation audit — current setup vs requirements (Figs 3–5, 9–10)

## 1. Summary of the CURRENT simulation setup

### 1. Sample type — **simulated SHOTGUN**, not 2bRAD
All simulated benchmarks use **shotgun/WGS-style reads**: random 150 bp fragments drawn uniformly from
each reference genome, 50 % reverse-complemented, **error-free**, quality all `I` (`scripts/sim_strain_mock.py`,
`sim_depth.py`, `sim_samples.py`, `run_species_expansion.py`). Strain2bScan then applies **in-silico
BcgI/multi-enzyme digestion** to those reads. Native 2bRAD is tested only on **real** data (mock Fig 6,
saliva Fig 7–8), never simulated. So Figs 3–5, 9–10 are all simulated shotgun.

### 2. Single-species (multi-strain) communities — current
Per-sample design (`sim_strain_mock.py`, and the identical inline sim in `run_species_expansion.py`):
- **Strains mixed:** random **2–5** strains per sample (uniform over {2,3,4,5}).
- **Genome IDs:** drawn from each species' pinned panel — `data/accessions/cae_panel_64.txt`,
  `saureus_panel.txt`, `sepidermidis_panel.txt`, `prevotella_copri_40x.txt`/`_112.txt`,
  `akkermansia_muciniphila_40x.txt`, `mycobacterium_tuberculosis_40x.txt` (Fig 10 uses subsets of
  StrainScan's own DB accession lists). Species covered: **C. acnes, S. aureus, S. epidermidis** (Fig 3/5),
  **P. copri** (Fig 4), **A. muciniphila, P. copri, M. tuberculosis** (Fig 10) — **6 species total**.
- **Ground-truth abundance:** per-strain **depth ~ log-normal(µ=1.0, σ=0.6), floored at 1×** (uneven);
  truth written as relative abundance = depth/Σdepth to `truth/sample*.truth.tsv`.
- **Samples:** 5 per species.

### 3. Multi-species (multi-strain) community — current (Fig 9)
`sim_samples.py`: **55 species × ~4 strains = 218 genomes** (`data/accessions/multispecies_55x4.tsv`),
**30 samples**, each mixing **12 species**; per species ~40 % get 2 strains else 1; per-strain
**depth ~ log-normal(µ=0.7, σ=0.5), floored at 1×**. Truth = species/strain/depth per sample. Depth
scaling (Fig 9C) reuses/cycles these samples across increasing sample counts.

### 4. Sequencing depth — current
- **Single-species (Fig 3A/4/5/10):** one draw of log-normal per-strain depth per sample; **not** a
  controlled depth ladder on a fixed community.
- **Depth sensitivity (Fig 3B):** a **single strain** (`GCF_009737125.1`) at **0.1× / 0.5× / 1× / 5×** —
  the only controlled depth series, and it is single-strain, not a multi-strain community.
- **Multi-species (Fig 9):** single log-normal depth regime; no depth gradient.

### Consistency — **NOT met currently**
Each figure simulates its **own** dataset with its **own** species/panel/seed:
C. acnes (Fig 3,5) · P. copri panels (Fig 4) · single C. acnes strain (Fig 3B) · 3 StrainScan species
(Fig 10) · 55-species community (Fig 9). There is **no single underlying simulated dataset** shared
across figures, and the reads are **error-free**.

---

## 2. StrainScan methodology (Liao et al., Microbiome 2023; PMID 37587527)

- **Read simulator:** `ART` — `art_illumina -p -l 250 -f <depth> -m 600 -s 150` (paired-end, 250 bp,
  insert 600±150, **with Illumina error model**).
- **Single-strain (their Fig 6):** per species, pick 1 strain, repeat 60× (50× for P. copri); depths
  **1×, 3×, 5×, 10×**. 6 species (C. acnes, E. coli, S. epidermidis, M. tuberculosis, A. muciniphila,
  P. copri). Min accepted depth 1×.
- **Co-existing strains (their Fig 7):** 2/3/5 strains, two difficulty strategies (strains from
  **different** clusters vs **same** cluster), fixed **uneven coverage profiles** —
  2→(100×,10×); 3→(100×,50×,10×); 5→(100×,70×,50×,20×,10×); 10 reps each → 60 sets/species.
- **Abundance metric:** Jensen–Shannon divergence (JSD) truth vs predicted.

---

## 3. Gap analysis — is re-simulation needed? **Yes.**

| Requirement | Current | Gap |
|---|---|---|
| One consistent dataset across figures | separate per figure | **re-simulate from one genome pool** |
| Single-species datasets for **all 15** Fig 2 species | only 6 species | **add 9 species** |
| Multi-strain, **uneven** abundance | ✅ (log-normal) | keep; adopt StrainScan fixed profiles for comparability |
| **Several depths per community** | only single-strain ladder | **add depth ladder per multi-strain community** |
| Multi-species across **3 depth gradients** | single regime | **add 3 depth gradients** |
| StrainScan-style reads (ART, errors, PE) | error-free 150 bp custom | **switch to ART** (adds realism; also closes the "error-free" critique) |

---

## 4. Proposed unified re-simulation (StrainScan-aligned, one consistent dataset)

**Master genome pool (pin once, reuse everywhere):** the **15 Fig-2 species**, complete/near-complete
genomes already downloaded (15–50 HQ genomes/species; `data/genome_qc_16s_panel.tsv`). This single pool
builds the Strain2bScan DBs, the Fig-2 distances, the single-species communities, and the multi-species
community — guaranteeing consistency.

**Read simulator:** `art_illumina -p -l 250 -f <depth> -m 600 -s 150` (paired-end, Illumina error model),
matching StrainScan. (Install: NIEHS ART macOS binary; fallback InSilicoSeq `iss` if ART won't build.)

**A. Single-species, multi-strain (15 datasets, one per Fig-2 species).** Per species, per cluster
(0.95) build the DB, then simulate communities:
- strain counts **2, 3, 5**; abundance = StrainScan uneven profiles (100:10 / 100:50:10 /
  100:70:50:20:10), normalized;
- **depth ladder 1×, 3×, 5×, 10×** (per StrainScan) applied to the *whole community* (scales all strains);
- 2 difficulty strategies (same-cluster / different-cluster) where the panel allows; **N=5 reps**.
- Truth: strain→relative abundance + cluster label, per sample.
- Deposit: `figure_raw_data/sim_single_species/<Species>/depth<X>/rep<r>.fq(.gz)` + `truth/…`.

**B. Multi-species, multi-strain community (1 dataset, 3 depth gradients).** Mix of the 15 species
(each contributing 1–2 strains at uneven abundance), simulated at **3 community depth gradients**
(e.g. low 1×, medium 5×, high 10× mean per-strain), **N samples** (e.g. 20). Truth: species/strain/RA.
Deposit: `figure_raw_data/sim_multi_species/depth_{low,med,high}/sample*.fq(.gz)` + truth.

**Figure re-mapping (all off the ONE dataset):**
- Fig 3A accuracy → single-species set at a fixed depth (e.g. 5×); Fig 3B depth → the depth ladder.
- Fig 4 robustness → subset panels of one species (P. copri) from the same pool.
- Fig 5 enzyme knob → one species (C. acnes) from the same pool.
- Fig 9 scale → the multi-species community + its depth gradients.
- Fig 10 head-to-head → the 3 StrainScan species (A. muciniphila, P. copri, M. tuberculosis), which are
  in the 15-species pool, at the ladder depths.

**Scale/compute estimate:** A = 15 species × 3 strain-counts × 4 depths × 5 reps ≈ 900 samples; at ART
250 bp PE, ~10× of a ~4 Mb community ≈ ~0.1–0.5 GB/sample → tens–~150 GB. B ≈ 60 samples. Feasible on the
~1 TB free; a few hours of ART. Reads git-ignored under `figure_raw_data/`, with pinned accessions +
truth committed.

**Decision points (please confirm):**
1. Simulator: **ART** (match StrainScan) vs InSilicoSeq (easier install) — recommend ART.
2. Depth ladder: **1/3/5/10×** (StrainScan) — add 0.5× to keep the low-depth story? (recommend add 0.5×).
3. Strain counts & abundance: adopt StrainScan **2/3/5 with fixed uneven profiles** (recommend yes).
4. Reps per condition: **5** (recommend), and same-cluster + different-cluster strategies (recommend yes).
5. Multi-species community: build from the **15-species pool** (recommend, for consistency) vs the
   existing 55-species panel.
