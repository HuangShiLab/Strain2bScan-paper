# Strain2bScan — detailed figure legends

Each legend has three parts: **(1) Data source** — data type (simulation / real metagenome / real 2bRAD),
production/collection method, and whether raw sequencing data is available locally
(`figure_raw_data/<Fig>/reads/`, git-ignored) or by accession; **(2) Key issue and conclusion**;
**(3) Results by subfigure**. Figure files are in `figures/` (numbered set in `figures/numbered/`);
source tables in `results/`; per-experiment docs in `docs/`.

Local raw-data index: `figure_raw_data/README.md`. Reference genomes are public assemblies (accession
manifests provided, not stored as reads).

---

## Figure 1 — Strain2bScan overview and the two input modes
**1. Data source.** Schematic; no sequencing data.
**2. Key issue & conclusion.** How one engine resolves strains from a sparse 2bRAD marker set and accepts
two input modes. Conclusion: reference genomes → 2bRAD tags → within-species clusters → cluster×marker DB;
a sample is digested once and profiled by a species gate (Layer-1) then within-species strain scoring
(Layer-2); the same tag space serves in-silico-digested shotgun and native BcgI 2bRAD, and is shared with
the Fast2bRAD-M species layer.
**3. Results by subfigure.** (A) Database construction pipeline. (B) Per-sample profiling pipeline and the
two data modes entering the shared tag space.

---

## Figure 2 — 2bRAD tags track genome-wide strain distance; 16S does not
**1. Data source.**
- *Data type:* reference genome assemblies (in-silico distance computation; no new sequencing).
- *Production/collection:* 15 pathogenic/commensal species, up to 50 genomes each pulled from NCBI
  (accession lists) via ENA FASTA, then **restricted to complete/near-complete assemblies** (CheckM
  completeness ≥ 97 %, contamination ≤ 5 %, assembly level Complete Genome/Chromosome), giving 14–50
  genomes/species. 16S genes extracted with barrnap 0.9 (HMMER `nhmmer`).
- *Raw data (local):* reference genomes are public — accession + QC manifest in
  `figure_raw_data/Fig02_16S_panel/genome_accessions_qc.tsv` (= `data/genome_qc_16s_panel.tsv`);
  per-pair distances in `results/pairdist/*.tsv`. No reads (assemblies only).
**2. Key issue & conclusion.** Does the reduced 2bRAD tag set, unlike the 16S gene, carry genome-wide
strain-level signal? Conclusion: **yes** — 2bRAD between-strain distances track whole-genome distance in
every species (median Spearman **0.94**), whereas 16S does not (median **0.36**), several 16S intervals
overlapping zero. 16S resolves species, 2bRAD tags resolve strains.
**3. Results by subfigure.**
- **(A)** Per-species Spearman correlation of 2bRAD (blue) and 16S (red) between-strain distance with
  whole-genome distance (bottom-3000 21-mer MinHash), sorted by 2bRAD advantage; 95 % CI error bars (500
  genome subsamples); genome count per species at right. Every 2bRAD CI is high and non-overlapping with
  its 16S CI; 16S CIs for *M. tuberculosis* (−0.04), *L. plantarum* (0.02) and *P. dorei* (−0.13) span
  zero (no strain signal).
- **(B)** 3×5 matrix of per-species **rank–rank** scatters: rank of each strain pair's whole-genome
  distance (x) vs its 2bRAD (blue) / 16S (red) distance (y), sorted worst→best 16S. 2bRAD hugs the
  diagonal in every species; 16S forms flat horizontal rank-bands (few discrete distances shared by many
  unrelated pairs), most extreme in *P. dorei* and *M. tuberculosis* (collapse to a single band).

---

## Figure 3 — Accurate strain profiling and depth sensitivity
**1. Data source.**
- *Data type:* simulation (short reads simulated from real reference genomes; closed-world).
- *Production/collection:* (A) real reference panels — *C. acnes* (64 genomes), *S. aureus* (60),
  *S. epidermidis* (60), NCBI accessions pinned — with simulated 2–5-strain mixtures at log-normal depth
  ≥1×; abundance evaluated at cluster resolution. (B) a single *C. acnes* strain simulated across a
  0.1–5× coverage ladder. All reads error-free 150 bp.
- *Raw data (local):* simulated reads are regenerated from `scripts/` + pinned accessions (not stored);
  accession lists in `data/accessions/`.
**2. Key issue & conclusion.** Is the sparse marker profiler accurate, and does the reduced tag set cost
low-depth sensitivity? Conclusion: **precision 1.0 across species** with high recall/low abundance error,
and **detection onset 0.5× coverage matching StrainScan** — no low-depth penalty.
**3. Results by subfigure.**
- **(A)** Precision (1.0 for all three), recall (0.75/0.79/1.0) and abundance accuracy (1−Bray–Curtis =
  0.76/0.67/0.98); cluster counts 16/17/10. Precision stays 1.0 even for low-diversity *S. epidermidis*.
- **(B)** Detection vs per-strain coverage for Strain2bScan (default and permissive floor) and StrainScan;
  all detect at ≥0.5×, all miss at 0.1× — curves coincide.

---

## Figure 4 — Reference-genome completeness limits strain ID under Jaccard; the `--containment` mode restores it
**1. Data source.**
- *Data type:* simulation — ART PE250 reads from the 15-species pool (same dataset as Figs 3/5/9/10);
  reference genomes are degraded assemblies. Two clustering modes compared.
- *Production/collection:* per species, a fixed multi-strain sample (2/3/5-strain uneven community at 5×,
  5 reps) is profiled while the truth strains' reference genomes are degraded across a completeness ladder
  **100/95/90/80/70/50 %** (contamination 0→10 %, fragmentation 1→400 contigs; `scripts/degrade.py`); the
  DB is rebuilt at each level under **Jaccard** and under **`--containment`** and the same reads re-profiled.
- *Raw data (local):* reads under `figure_raw_data/sim_single_species/`; genome pool by accession
  (`sim_pool_manifest.tsv`); results `results/refqual_15species.tsv` (Jaccard) and
  `results/refqual_15species_containment.tsv`.
**2. Key issue & conclusion.** Does reference incompleteness break strain ID, and does the tool address it?
Conclusion: under default **Jaccard**, precision and recall decline as references degrade (median precision
1.0→0.84→0.71, recall 0.96→0.80→0.74 at 100→90→70 %) because incomplete genomes split from complete
relatives. The **`--containment`** mode (max-containment) keeps them clustered, restoring median precision
to 0.98/0.92 and recall to 0.95/0.92 at 95/90 % completeness, converging with Jaccard only at ≤70 %
(genuinely low-quality). It also removes the near-clonal *M. tuberculosis* artifact. Containment is opt-in
(merges more aggressively); default stays Jaccard + the assembly-quality filter.
**3. Results by subfigure.**
- **(A)** Schematic of the mechanism: an incomplete genome is a subset of a complete relative → Jaccard
  splits them (shared tags demoted to non-discriminating *SharedPartial*) → max-containment merges them
  (one cluster, marker set = union of members → discriminating tags preserved).
- **(B)** Median precision vs completeness: Jaccard (grey dashed) vs `--containment` (solid), shaded gap;
  faint per-species containment lines show the spread.
- **(C)** Median recall vs completeness, same layout. **Inset:** *M. tuberculosis* recall — Jaccard collapses
  to ≈0.05 on any degradation (its single cluster shatters into spurious singletons), `--containment` holds
  1.0 to 90 % (the near-clonal cluster-fragmentation artifact is fixed).

---

## Figure 5 — The 2bRAD enzyme set is a resolution/cost knob enabling native BcgI operation
**1. Data source.**
- *Data type:* simulation (real *C. acnes* panel + simulated mixtures; closed-world).
- *Production/collection:* the *C. acnes* benchmark re-profiled while varying the type-IIB enzyme set from
  1 (BcgI) to 14 enzymes.
- *Raw data (local):* regenerated from scripts + pinned accessions; not stored as reads.
**2. Key issue & conclusion.** How many enzymes are needed, and can single-enzyme BcgI operation support
native 2bRAD libraries? Conclusion: **BcgI alone gives precision 1.0**; ~4 enzymes is the recall sweet
spot; more enzymes add markers/compute but not clusters — and single-enzyme BcgI is exactly what lets
native BcgI 2bRAD-M libraries be profiled directly.
**3. Results by subfigure.** Precision/recall, strain-specific marker count and build time vs number of
enzymes (recall 0.56→0.75 by 4 enzymes then plateau; markers 882→1 524→7 400; build 0.6→1.8→6.4 s; cluster
count constant at 16).

---

## Figure 6 — DNA mock, host-contamination series: native 2bRAD holds strain recovery where shotgun collapses
**1. Data source.**
- *Data type:* **real 2bRAD data + real shotgun metagenome** (known-truth mock community).
- *Production/collection:* ATCC **MSA-1002** 20-strain even mock genomic DNA, prepared as native BcgI
  **2bRAD-M** libraries and as **shotgun (WMS)** libraries, spiked with human DNA to **0/90/99 %**
  contamination (SRA BioProject **PRJNA1131785**; 2bRAD on Illumina, WMS on DNBSEQ). Profiled against a
  reconstructed 62-species BcgI panel (20 mock + 42 decoys) and scored vs the 20-species truth.
- *Raw data (local):* **available** — `figure_raw_data/Fig06_mock_hostcontam/reads/` (4 native-2bRAD +
  3 shotgun R1/R2 FASTQ) with `manifest.tsv` (SRR ↔ contamination level). Table `results/mock_hostcontam.tsv`.
**2. Key issue & conclusion.** Does native 2bRAD's wet-lab reduction preserve strain recovery under heavy
host contamination where shotgun cannot? Conclusion: **precision 1.0 for both data types at every host
level**, but native 2bRAD keeps **20/20 recall at 99 % host** while shotgun drops to **12/20** — because
2bRAD preserves ~10× more usable markers under host load.
**3. Results by subfigure.**
- **(A)** Species recall vs % human DNA (0/90/99) for native 2bRAD vs in-silico-digested shotgun;
  precision annotated 1.0 throughout; 2bRAD flat at 20/20, shotgun falls to 12/20 at 99 %.
- **(B)** Usable BcgI marker yield vs % host (log): 2bRAD ~340–370k regardless of host; shotgun
  96k→53k→33k as host rises — the mechanism of the recall gap.

---

## Figure 7 — Real saliva: strain profiles discriminate individuals and are temporally stable
**1. Data source.**
- *Data type:* **real native BcgI 2bRAD metagenome** (human saliva).
- *Production/collection:* saliva from **8 subjects sampled at 4 times of day** (9AM/11AM/1PM/5PM = 32
  native-2bRAD libraries), SRA **PRJNA1131785** (Illumina). Profiled against a 19-species oral-commensal
  reference panel (up to 25 genomes/species). Strain/species relative-abundance matrices → Bray–Curtis →
  PERMANOVA (subject/timepoint) and leave-one-timepoint-out 1-NN classification.
- *Raw data (local):* **available** — `figure_raw_data/Fig07_saliva_2bRAD/reads/` (32 R1 FASTQ) with
  `manifest.tsv` decoding `S<time>-<subject>` (prefix = time of day, suffix = subject 1–14). Tables
  `results/saliva_permanova.tsv`, `saliva_perspecies_subject.tsv`, `saliva_temporal_ml.tsv`,
  `saliva_strain_long.tsv`.
**2. Key issue & conclusion.** Does strain-level 2bRAD profiling distinguish individuals better than
species-level, and is the signature stable? Conclusion: **yes** — subject strain-level PERMANOVA
**R² 0.833 > species 0.822** (p 2e-4), host-ID **100 % (strain) vs 94 % (species)**, time-of-day
non-significant, and within-subject ≪ between-subject distance: an individual-specific, temporally stable
strain signature resolved in ~1 s/sample.
**3. Results by subfigure.**
- **(A)** Species-level PCoA coloured by subject (R² 0.82, 1-NN 88 %).
- **(B)** Strain-level PCoA coloured by subject (R² 0.83, 1-NN 91 %) — tighter within-subject clusters.
- **(C)** Per-species strain-level subject R² (bars); *Rothia mucilaginosa* 0.921 highest, 17/18 species
  significant.
- **(D)** Within- vs between-subject strain Bray–Curtis distance (0.19 vs 0.44, p 5e-25) — temporal stability.
- **(E)** Leave-one-timepoint-out host-ID accuracy, strain (100 %) vs species (94 %) vs chance (12.5 %).

---

## Figure 8 — Native 2bRAD confirms all shotgun strains and recovers the low-abundance ones shotgun misses
**1. Data source.**
- *Data type:* **real paired shotgun metagenome (WMS) + real native 2bRAD** (same saliva samples).
- *Production/collection:* paired shotgun WMS of the Fig 7 saliva samples (SRA PRJNA1131785, DNBSEQ). WMS
  R1 (~14 GB/sample) was truncated on download and used as a valid prefix subsample (in-silico BcgI
  digestion), compared per sample to the matched native-2bRAD strain calls; 3 subjects (7/1/14) with
  usable bacterial recovery.
- *Raw data (local):* **available (partial)** — `figure_raw_data/Fig08_saliva_shotgun/reads/` (WMS R1
  prefixes) with `manifest.tsv`; cached WMS profiles in `results/wms_preds/`; table
  `results/saliva_concordance.tsv`. (Full-depth WMS re-downloadable from PRJNA1131785.)
**2. Key issue & conclusion.** On high-host saliva, does native 2bRAD recover the strains shotgun finds
*and* the low-abundance strains shotgun cannot? Conclusion: **100 % of shotgun strain calls (65/65) are
confirmed by native 2bRAD** (the two modes agree; no in-silico false positives), and native 2bRAD
additionally recovers **128–163 strains/sample** that shotgun misses, which are **significantly
lower-abundance** — quantifying 2bRAD's sensitivity advantage on real clinical material.
**3. Results by subfigure.**
- **(A)** Per sample, strains detected: shotgun-confirmed (blue, = 100 % of shotgun calls) plus 2bRAD-only
  (orange, 128–163) — native 2bRAD is a strict superset.
- **(B)** Community relative abundance (log) of shared vs 2bRAD-only strains; 2bRAD-only significantly
  lower (median 0.0029 vs 0.0097; Mann–Whitney p 1.2e-23).

---

## Figure 9 — Fast, light, and scalable to whole communities
**1. Data source.**
- *Data type:* mixed — real *C. acnes* benchmark (A, B) + simulation (C, 55-species community).
- *Production/collection:* (A) per-sample time/memory on the real *C. acnes* mock vs StrainScan; (B) thread
  scaling of build and profile; (C) a simulated 55-species community (218 NCBI genomes, 30 samples mixing
  12 species at log-normal depth) profiled across increasing sample counts.
- *Raw data (local):* *C. acnes* mock reads from MockMetagenomes4Benchmark; simulated community regenerated
  from scripts + `data/accessions/multispecies_55x4.tsv`; not stored as reads.
**2. Key issue & conclusion.** How does the sparse-marker, digest-once design scale on conventional
metagenomes? Conclusion: **~8× faster / ~11× lighter per sample**, clean thread parallelism, and — because
per-sample cost is flat in #species — **~121–146× faster on a 55-species community** than a per-species
k-mer tool run once per species.
**3. Results by subfigure.**
- **(A)** Per-sample wall-time and peak memory, Strain2bScan (0.86 s / 78 MB) vs StrainScan (7.06 s / 828 MB).
- **(B)** Speedup of build and profile vs thread count (8.5× / 6.5× at 16 threads).
- **(C)** Community throughput: per-sample time flat in species count; fold-speedup vs a per-species tool
  (~132× at 100 samples, ~146× at ≥200).

---

## Figure 10 — Matches or exceeds StrainScan on its own databases
**1. Data source.**
- *Data type:* simulation on real reference panels (StrainScan's own reference sets; closed-world).
- *Production/collection:* *A. muciniphila*, *P. copri* and near-clonal *M. tuberculosis* (40-genome
  subsets), profiled head-to-head with StrainScan (v1.0) on the same panels; time/memory measured with
  `/usr/bin/time -l` on a 16-core Apple-silicon machine.
- *Raw data (local):* regenerated from scripts + pinned accessions; StrainScan DB build script provided
  (Linux-only dependency).
**2. Key issue & conclusion.** Head-to-head, does Strain2bScan match StrainScan's accuracy at far lower
cost? Conclusion: **matched precision 1.0**, **matched/exceeded recall** (0.93 vs 0.24; 0.94 vs 0.90),
**~17–23× faster / ~15–24× lighter**, and **completes near-clonal *M. tuberculosis* in ~1 s where
StrainScan does not** (>3.3 h, >25 GB).
**3. Results by subfigure.** Per species: precision, recall, profile time and memory for both tools;
*M. tuberculosis* marked StrainScan-DNF.

---

## Figure S1 — ATCC MSA-1002 DNA-input titration (native 2bRAD)
**1. Data source.** Real native BcgI 2bRAD of the ATCC MSA-1002 mock across a DNA-input titration
(0.001–100 ng); reads from **Figshare article 12272360** (2B-RAD-M unique-tags DB). Raw data: Figshare
(accession above); table `results/mock_msa1002_titration.tsv`.
**2. Key issue & conclusion.** Low-biomass limit of native 2bRAD. Conclusion: precision 1.0 with full
species recall down to **0.1 ng** input (19/20 at 0.01 ng).
**3. Results.** Species precision (1.0 throughout) and recall vs DNA input (log), annotated with read count
and detections per level.

## Figure S2 — Layer-1 gate calibration on the 55-species panel
**1. Data source.** Simulation (55-species community, normal and low ~0.62× depth). Raw data: regenerated
from scripts; tables `results/gate_calibration_*.tsv`.
**2. Key issue & conclusion.** Choosing the species-gate floor. Conclusion: the default floor gives species
precision 1.0 at both depths with leakage held in the middle tier; the breadth term is scale insurance.
**3. Results.** Species precision/recall vs marker floor at normal and low depth.

## Table S3 — Exploratory clinical oral cohort profiling
**1. Data source.** Real native BcgI 2bRAD of 4 clinical oral samples (SRA PRJNA1131785 `S_#` cohort).
Raw data: **available** — `figure_raw_data/TableS3_clinical_oral/reads/` with `manifest.tsv`; table
`results/clinical_exploratory.tsv`. No case/control labels in public metadata.
**2. Key issue & conclusion.** Does the tool run on the clinical cohort? Conclusion: yes — 15–17 species,
115–158 strain calls, ~2 s/sample; a differential tumour/normal test needs the study's sample labels.
**3. Results.** Per sample: marker count, species resolved, strain calls, runtime, top species.

## Table S4 — Genome quality-control for the 16S/2bRAD motivation panel
**1. Data source.** Reference-assembly metadata (NCBI). File: `data/genome_qc_16s_panel.tsv` /
`figure_raw_data/Fig02_16S_panel/genome_accessions_qc.tsv`.
**2. Key issue & conclusion.** Documents the completeness filter behind Fig 2. Conclusion: only
complete/near-complete genomes (CheckM ≥97 %/≤5 %, Complete/Chromosome) were used.
**3. Results.** Per genome: species, accession, assembly level, CheckM completeness/contamination, contig
count, length, high-quality flag.
