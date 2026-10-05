# Strain2bScan; detailed figure legends

Each legend has three parts: **(1) Data source**. data type (simulation / real metagenome / real 2bRAD),
production/collection method, and whether raw sequencing data is available locally
(`figure_raw_data/<Fig>/reads/`, git-ignored) or by accession; **(2) Key issue and conclusion**. **(3) Results by subfigure**. Figure files are in `figures/` (numbered set in `figures/numbered/`);
source tables in `results/`; per-experiment docs in `docs/`.

Local raw-data index: `figure_raw_data/README.md`. Reference genomes are public assemblies (accession
manifests provided, not stored as reads).

---

## Figure 1; Strain2bScan algorithm overview
**1. Data source.** Schematic; no sequencing data. Panel C plots the closed-form expectation
`1 - e^(-lambda)` and two measured points from the shadow-cluster experiment.
**2. Key issue & conclusion.** How one engine resolves strains from a ~1-2 % genomic subsample,
accepts both in-silico-digested shotgun and native BcgI 2bRAD-M through the same tag space, and
decides which calls to keep. Conclusion: reference genomes -> single-copy 2bRAD tags -> within-
species clusters -> cluster x marker database; a sample is digested **once** and matched against
every species database, gated at species level (Layer-1), restricted to panel-wide
species-specific markers, then scored within species (Layer-2); abundance is reported at three
scopes because per-species fractions cannot be concatenated into a community composition.
**3. Results by subfigure.** (A) Database construction, with the clustering threshold and the
occurrence-based definition of a unique marker. (B) Per-sample profiling: the two input modes, the
species gate with its depth-reachability scaling, the cross-species marker restriction, the three
Layer-2 detection criteria, the zero-inclusive depth estimator, and the three abundance scopes.
(C) The discriminator that a coverage floor cannot provide. Under Poisson sampling a genuinely
present cluster at depth lambda must show breadth `1 - e^(-lambda)` (blue curve). A genuine rare
strain at 0.44x with breadth 0.39 sits on the curve; a shadow -- a cluster called because the
strain in the sample carries a fraction of its distinguishing loci, so those markers appear at
full depth across only part of the panel -- sits far below it at 7.7x with breadth 0.35. The two
have near-identical breadth and differ 17x in depth, which is why the ratio, not the breadth, is
the test. Both points are measured.

---

## Figure 2; 2bRAD tags track genome-wide strain distance; 16S does not
**1. Data source.**
- *Data type:* reference genome assemblies (in-silico distance computation; no new sequencing).
- *Production/collection:* 15 pathogenic/commensal species, up to 50 genomes each pulled from NCBI
  (accession lists) via ENA FASTA, then **restricted to complete/near-complete assemblies** (CheckM
  completeness ≥ 97 %, contamination ≤ 5 %, assembly level Complete Genome/Chromosome), giving 14–50
  genomes/species. 16S genes extracted with barrnap 0.9 (HMMER `nhmmer`).
- *Raw data (local):* reference genomes are public; accession + QC manifest in
  `figure_raw_data/Fig02_16S_panel/genome_accessions_qc.tsv` (= `data/genome_qc_16s_panel.tsv`);
  per-pair distances in `results/pairdist/*.tsv`. No reads (assemblies only).
**2. Key issue & conclusion.** Does the reduced 2bRAD tag set, unlike the 16S gene, carry genome-wide
strain-level signal? Conclusion: **yes**. 2bRAD between-strain distances track whole-genome distance in
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

## Figure 3; Accurate strain profiling and depth sensitivity
**1. Data source.**
- *Data type:* simulation (short reads simulated from real reference genomes; closed-world).
- *Production/collection:* (A) real reference panels; *C. acnes* (64 genomes), *S. aureus* (60),
  *S. epidermidis* (60), NCBI accessions pinned; with simulated 2–5-strain mixtures at log-normal depth
  ≥1×; abundance evaluated at cluster resolution. (B) a single *C. acnes* strain simulated across a
  0.1–5× coverage ladder. Reads were ART error-modelled 250 bp paired-end reads.
- *Raw data (local):* simulated reads are regenerated from `scripts/` + pinned accessions (not stored);
  accession lists in `data/accessions/`.
**2. Key issue & conclusion.** Is the sparse marker profiler accurate, and does the reduced tag set cost
low-depth sensitivity? Conclusion: **precision 1.0 across species** with high recall/low abundance error,
and **detection onset 0.5× coverage matching StrainScan**. no low-depth penalty.
**3. Results by subfigure.**
- **(A)** Precision (1.0 for all three), recall (0.75/0.79/1.0) and abundance accuracy (1−Bray–Curtis =
  0.76/0.67/0.98); cluster counts 16/17/10. Precision stays 1.0 even for low-diversity *S. epidermidis*.
- **(B)** Detection vs per-strain coverage for Strain2bScan (default and permissive floor) and StrainScan;
  all detect at ≥0.5×, all miss at 0.1×; curves coincide.

---

## Figure 4; Reference-genome completeness limits strain ID under Jaccard; the `--containment` mode restores it
**1. Data source.**
- *Data type:* simulation; ART PE250 reads from the 15-species pool (same dataset as Figs 3/5/9/10);
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
- **(C)** Median recall vs completeness, same layout. **Inset:** *M. tuberculosis* recall; Jaccard collapses
  to ≈0.05 on any degradation (its single cluster shatters into spurious singletons), `--containment` holds
  1.0 to 90 % (the near-clonal cluster-fragmentation artifact is fixed).

---

## Figure 5; The 2bRAD enzyme set is a resolution/cost knob (14-species multi-species ladder)
**1. Data source.**
- *Data type:* simulation on the 15-species reference pool (closed-world, error-modelled ART 250 bp paired-end reads).
- *Production/collection:* 14 resolvable species of the simulation pool were re-profiled while varying the
  type-IIB enzyme set from 1 (BcgI) to 2/4/8/14 enzymes. Each species was sampled with 2/3/5 co-present
  strains at log-normal depth ≥1×; metrics are medians over species/replicates.
- *Raw data (local):* `results/enzyme_sweep_multi.tsv`; regenerated from `scripts/run_enzyme_sweep_multi.py`
  + pinned accessions; not stored as reads.
**2. Key issue & conclusion.** How many enzymes are needed across a broad species panel, and can
single-enzyme BcgI operation support native 2bRAD-M libraries? Conclusion: **precision stays 1.0 at every
step in every species**. Median recall rises monotonically with enzyme count; BcgI alone 0.50 → ~4 enzymes
0.72 → 14 enzymes 1.00; cluster count is invariant (more enzymes add marker density/recall, not
resolution). Single-enzyme BcgI operation is therefore sufficient for native BcgI 2bRAD-M libraries, while
~4 enzymes is the practical recall/cost sweet spot for in-silico shotgun profiling.
**3. Results by subfigure.**
- **(A)** Median precision (1.0 throughout) and recall vs enzyme number, with faint per-species recall
  trajectories. Recall is species-dependent: diverse species saturate early, near-clonal species need more
  enzymes.
- **(B)** Strain-specific marker yield vs enzyme number (median ~1 000 → ~12 000 markers from BcgI to
  14 enzymes); cluster count does not increase.

---

## Figure 6; Native 2bRAD strain-level identification and abundance across four ATCC DNA mocks
**1. Data source.**
- *Data type:* **real native BcgI 2bRAD-M metagenome** of four ATCC whole-cell mocks.
- *Production/collection:* MSA-1002 (20 strains, even), MSA-1003 (20 strains, staggered ~3 orders of
  magnitude), MSA-1005 and MSA-1007 (6 strains each). Native BcgI 2bRAD reads were profiled against a
  unified 28-species / 164-genome combined tree (each mock species = ATCC genome + up to 5 conspecific
  decoys, within-species ANI > 95 %), with `--enzyme BcgI --containment --similarity 0.95`,
  `--min-abundance 0 --min-coverage 0.2`. MSA-1002 was additionally run across a host-DNA ladder
  (90/95/99/99.9 %) and a low-biomass DNA-input ladder (1 → 0.001 ng); SRA PRJNA1131785.
- *Raw data (local):* **available**. `figure_raw_data/Fig06_mock_hostcontam/reads/` with `manifest.tsv`;
  table `data/fig6_fig12_metrics.tsv`; scorer `scripts/score_all.py`.
**2. Key issue & conclusion.** Does native BcgI 2bRAD-M resolve and quantify individual strains on
real DNA mocks, and how robust is it to host contamination and low input? Conclusion: **yes**. Strain2bScan
At 1e-4, F1 is 0.952 at 0.1 ng and 0.625 at 0.01 ng (precision 0.909 and 0.833; recall 1.0 and 0.5); F1 is
≈0.95–1.0 down to 99 % host DNA. Staggered mocks show a ~1× marker-depth noise floor that costs
single-threshold precision but is recovered by abundance-threshold AUPR.
**3. Results by subfigure.**
- **One row per sample.** Left: stacked per-genome relative abundance (truth vs Strain2bScan; false
  positives/decoys hatched grey). Right: precision, recall, F1, AUPR (abundance-threshold sweep),
  Bray–Curtis and L2 vs sequence-abundance ground truth.
- MSA-1002 rows: host-DNA ladder (90/95/99/99.9 %) and DNA-input ladder (1/0.1/0.01/0.001 ng).
- MSA-1003/1005/1007 rows: triplicates; MSA-1003 shows the ~1× marker-depth noise floor on the 28-species
  tree that the 20-species tree and AUPR sweep recover.

---

## Figure 12; Strain-level profiling on shotgun mocks: Strain2bScan vs StrainScan and inStrain,
including high host contamination
**1. Data source.**
- *Data type:* **whole-metagenome shotgun (WMS)** of the same four ATCC mocks as Fig 6.
- *Production/collection:* in-silico all-enzyme digestion of WMS reads against the 164-genome combined tree
  for Strain2bScan; StrainScan v1.0.14 with its per-species databases; inStrain 1.10.0 with a 98 %-ANI
  dereplicated reference (the non-dereplicated control is shown in Fig S4). Each tool scored in its own
  0.95-similarity cluster space against the mock ground truth. Host-contamination ladder 0/90/95/99 %
  human DNA for MSA-1002.
- *Raw data (local):* `figure_raw_data/Fig06_mock_hostcontam/reads/` WMS subsets; table
  `data/fig6_fig12_metrics.tsv`; scorer `scripts/score_all.py`; doc `docs/mock_hostcontam.md`.
**2. Key issue & conclusion.** On conventional shotgun input, does Strain2bScan match StrainScan/inStrain
accuracy while preserving detection *and* quantification under host contamination? Conclusion:
At the primary 1e-4 threshold, **Strain2bScan is the only tested tool that preserves both** (F1 = 1.0, Bray–Curtis dissimilarity 0.302, similarity
≥ 0.72). StrainScan keeps detection but its depth estimator diverges (Bray–Curtis similarity 0.03 at
99 % host); inStrain loses detection (recall 0.20). Clean samples are concordant across all three tools.
**3. Results by subfigure.**
- **One row per sample**, same layout as Fig 6: left = stacked per-genome abundance (truth + each tool),
  right = precision, recall, F1, AUPR, Bray–Curtis, L2.
- Core row block: MSA-1002 host-DNA ladder (0/90/95/99 %) showing the three-way separation.
- MSA-1003 rows: staggered mock on 20-species (`120`) and 28-species (`164`) trees.
- MSA-1005/1007 rows: 6-strain mocks at 0 % host, confirming clean-sample concordance.

---

## Figure 7; Real saliva: strain profiles discriminate individuals and are temporally stable
**1. Data source.**
- *Data type:* **real native BcgI 2bRAD-M metagenome** (human saliva).
- *Production/collection:* saliva from **8 subjects sampled at 4 times of day** (9AM/11AM/1PM/5PM = 32
  native-2bRAD libraries), SRA **PRJNA1131785** (Illumina). Profiled against a 19-species oral-commensal
  reference panel (up to 25 genomes/species). Strain/species relative-abundance matrices → Bray–Curtis →
  PERMANOVA (subject/timepoint) and leave-one-timepoint-out 1-NN classification.
- *Raw data (local):* **available**. `figure_raw_data/Fig07_saliva_2bRAD/reads/` (32 R1 FASTQ) with
  `manifest.tsv` decoding `S<time>-<subject>` (prefix = time of day, suffix = subject 1–14). Tables
  `results/saliva_permanova.tsv`, `saliva_perspecies_subject.tsv`, `saliva_temporal_ml.tsv`,
  `saliva_strain_long.tsv`.
**2. Key issue & conclusion.** Does strain-level 2bRAD profiling distinguish individuals better than
species-level, and is the signature stable? Conclusion: **yes**. subject strain-level PERMANOVA
**R² 0.757 versus species 0.755** (p 2e-4), host-ID **100 % (strain) vs 78.1 % (species)**, time-of-day
non-significant, and within-subject ≪ between-subject distance: an individual-specific, temporally stable
strain signature resolved in ~1 s/sample.
**3. Results by subfigure.**
- **(A)** Species-level PCoA coloured by subject (R² 0.755, 1-NN 78.1 %).
- **(B)** Strain-level PCoA coloured by subject (R² 0.757, 1-NN 96.9%).
- **(C)** Per-species strain-level subject R² (bars); *Neisseria subflava* 0.808 highest, 13/13 testable
  species significant.
- **(D)** Within- vs between-subject strain Bray–Curtis distance (0.327 vs 0.696, p 8.05e-19); temporal stability.
- **(E)** Leave-one-timepoint-out host-ID accuracy, strain (100 %) vs species (78.1 %) vs chance (12.5 %).

---

## Figure 8; Native 2bRAD confirms all shotgun strains and recovers the low-abundance ones shotgun misses
**1. Data source.**
- *Data type:* **real paired shotgun metagenome (WMS) + real native BcgI 2bRAD-M** (same saliva samples).
- *Production/collection:* paired shotgun WMS of the Fig 7 saliva samples (SRA PRJNA1131785, DNBSEQ). WMS
  R1 (~14 GB/sample) was truncated on download and used as a valid prefix subsample (in-silico BcgI
  digestion), compared per sample to the matched native-2bRAD strain calls; 3 subjects (7/1/14) with
  usable bacterial recovery.
- *Raw data (local):* **available (partial)**. `figure_raw_data/Fig08_saliva_shotgun/reads/` (WMS R1
  prefixes) with `manifest.tsv`; cached WMS profiles in `results/wms_preds/`; table
  `results/saliva_concordance.tsv`. (Full-depth WMS re-downloadable from PRJNA1131785.)
**2. Key issue & conclusion.** On high-host saliva, does native BcgI 2bRAD-M recover the strains shotgun finds
*and* the low-abundance strains shotgun cannot? Conclusion: All 65 shotgun strain-cluster calls (three usable paired samples, prefix-subsampled shotgun) are
present in native BcgI 2bRAD-M; this is directional concordance, not independent validation. Native BcgI 2bRAD-M
yields **128–163 additional candidate strain-cluster calls/sample**, which are **significantly
lower-abundance**. quantifying 2bRAD's sensitivity advantage on real clinical material.
**3. Results by subfigure.**
- **(A)** Per sample, strains detected: shared with shotgun (blue) plus 2bRAD-only candidate calls (orange)
  (orange, 128–163); this is a call-set concordance result and does not establish independent strain truth.
- **(B)** Community relative abundance (log) of shared vs 2bRAD-only strains; 2bRAD-only significantly
  lower (median 0.0029 vs 0.0097; Mann–Whitney p 1.2e-23).

---

## Figure 9; Fast, light, and scalable to whole communities
**1. Data source.**
- *Data type:* mixed; real *C. acnes* benchmark (A, B) + simulation (C, 55-species community).
- *Production/collection:* (A) per-sample time/memory on the real *C. acnes* mock vs StrainScan; (B) thread
  scaling of build and profile; (C) a simulated 55-species community (218 NCBI genomes, 30 samples mixing
  12 species at log-normal depth) profiled across increasing sample counts.
- *Raw data (local):* *C. acnes* mock reads from MockMetagenomes4Benchmark; simulated community regenerated
  from scripts + `data/accessions/multispecies_55x4.tsv`; not stored as reads.
**2. Key issue & conclusion.** How does the sparse-marker, digest-once design scale on conventional
metagenomes? Conclusion: **~8× faster / ~11× lighter per sample**, clean thread parallelism, and; because
per-sample cost is flat in #species; measured Strain2bScan cost is **121–146× lower than projected per-species querying** on a 55-species community.
**3. Results by subfigure.**
- **(A)** Per-sample wall-time and peak memory, Strain2bScan (0.86 s / 78 MB) vs StrainScan (7.06 s / 828 MB).
- **(B)** Speedup of build and profile vs thread count (4.6× / 5.8× at 16 threads; tracked in `results/parallel_and_build_scaling.tsv`).
- **(C)** Community throughput: per-sample time flat in species count; fold-speedup vs a per-species tool
  (~132× at 100 samples, ~146× at ≥200).

---

## Figure 10; Matches or exceeds StrainScan on its own databases
**1. Data source.**
- *Data type:* simulation on real reference panels (StrainScan's own reference sets; closed-world).
- *Production/collection:* *A. muciniphila*, *P. copri* and near-clonal *M. tuberculosis* (40-genome
  subsets), profiled head-to-head with StrainScan v1.0.14 on the same panels; time/memory measured with
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

## Figure 11; Systematic head-to-head on the 15-species simulated benchmark
**1. Data source.**
- *Data type:* simulation on a fixed 15-species genome pool (15–50 complete/near-complete NCBI genomes
  each; closed-world, error-modelled ART reads, 150 bp PE).
- *Production/collection:* single-species communities (2/3/5 co-present strains, same- or
  different-cluster, 0.5–10× per-strain ladder, 5 reps = 2 025 samples) and multi-species communities
  (~18 species, 3 depth gradients = 60 samples). Both tools build databases from the **same** pool and
  profile the **same** reads. Strain2bScan native (arm64); StrainScan v1.0.14 under Docker `linux/amd64`
  (QEMU) with a `linux/amd64` Strain2bScan run in the same container for the emulation-free speed ratio.
  Time/memory with `/usr/bin/time` (`-l` native, `-v` in container).
- *Scoring:* each tool in its own 0.95-cluster space (truth strains mapped to the tool's clusters);
  precision/recall/F1 per sample over cluster sets.
- *Raw data (local):* `figure_raw_data/sim_headtohead/*_persample.tsv`; aggregates
  `results/sim_headtohead_*.tsv`; tables `manuscript/tables.md`.
**2. Key issue & conclusion.** On a common, controlled benchmark, does Strain2bScan match StrainScan's
accuracy at far lower cost across many species and depths? Conclusion: **both hold precision 1.0**
(Strain2bScan in every species; StrainScan 0.80–0.83 in four), Strain2bScan **reaches full recall by 3×
vs StrainScan's 10×** (median R 0.80 vs 0.67, F1 0.89 vs 0.75), builds databases **249–614× faster /
43–138× lighter**, profiles **4–33× faster** (same env), and is **46–105× faster on multi-species
communities**. StrainScan **failed to build *K. pneumoniae*** entirely, which Strain2bScan built in 5.1 s.
**3. Results by subfigure.** (A–C) single-species precision/recall/F1 vs depth, both tools (14 species,
204 paired samples). (D) per-species DB build time (log), Strain2bScan vs StrainScan. (E) per-sample
profile time vs depth, same emulated container. (F) multi-species profiling time per community sample;
Strain2bScan one pass vs StrainScan Σ per-species runs, annotated with the fold-difference.

---

## Figure S1; ATCC MSA-1002 DNA-input titration (native BcgI 2bRAD-M)
**1. Data source.** Real native BcgI 2bRAD-M of the ATCC MSA-1002 mock across a DNA-input titration
(0.001–100 ng); reads from **Figshare article 12272360** (2B-RAD-M unique-tags DB). Raw data: Figshare
(accession above); table `results/mock_msa1002_titration.tsv`.
**2. Key issue & conclusion.** Low-biomass limit of native BcgI 2bRAD-M. Conclusion: precision 1.0 with full
species recall down to **0.1 ng** input in the separate 62-species titration panel (19/20 at 0.01 ng under that panel-specific gate).
**3. Results.** Species precision (1.0 throughout) and recall vs DNA input (log), annotated with read count
and detections per level.

## Figure S2; Layer-1 gate calibration on the 55-species panel
**1. Data source.** Simulation (55-species community, normal and low ~0.62× depth). Raw data: regenerated
from scripts; tables `results/gate_calibration_*.tsv`.
**2. Key issue & conclusion.** Choosing the species-gate floor. Conclusion: the default floor gives species
precision 1.0 at both depths with leakage held in the middle tier; the breadth term is scale insurance.
**3. Results.** Species precision/recall vs marker floor at normal and low depth.

## Table S3; Exploratory clinical oral cohort profiling
**1. Data source.** Real native BcgI 2bRAD-M of 4 clinical oral samples (SRA PRJNA1131785 `S_#` cohort).
Raw data: **available**. `figure_raw_data/TableS3_clinical_oral/reads/` with `manifest.tsv`; table
`results/clinical_exploratory.tsv`. No case/control labels in public metadata.
**2. Key issue & conclusion.** Does the tool run on the clinical cohort? Conclusion: yes; 15–17 species,
115–158 strain calls, ~2 s/sample; a differential tumour/normal test needs the study's sample labels.
**3. Results.** Per sample: marker count, species resolved, strain calls, runtime, top species.

## Table S4; Genome quality-control for the 16S/2bRAD motivation panel
**1. Data source.** Reference-assembly metadata (NCBI). File: `data/genome_qc_16s_panel.tsv` /
`figure_raw_data/Fig02_16S_panel/genome_accessions_qc.tsv`.
**2. Key issue & conclusion.** Documents the completeness filter behind Fig 2. Conclusion: only
complete/near-complete genomes (CheckM ≥97 %/≤5 %, Complete/Chromosome) were used.
**3. Results.** Per genome: species, accession, assembly level, CheckM completeness/contamination, contig
count, length, high-quality flag.


---

## Figure S3; Cost of unified-database expansion (20- vs 28-species combined tree)
- *Data type:* whole-metagenome shotgun (WMS) ATCC mock communities.
- *Comparison:* Strain2bScan on the 20-species (`120`) versus 28-species (`164`) combined tree.
- *Samples:* MSA-1002 at 0/90/95/99% host and MSA-1003 triplicates.
- *Result:* MSA-1002 is unchanged across database sizes; MSA-1003 loses single-threshold precision
  on the larger tree, while AUPR remains approximately 0.96.
- *File:* `figures/numbered/FigS3_tree_expansion.{png,pdf}`.

---

## Figure S4; inStrain requires a dereplicated reference (MSA-1002 shotgun)
- *Data type:* whole-metagenome shotgun ATCC MSA-1002 host-contamination ladder.
- *Comparison:* inStrain on the non-dereplicated 164-genome reference versus a 98%-ANI
  dereplicated reference.
- *Result:* the non-dereplicated reference inflates false positives; dereplication restores
  detection and precision. This is why Fig 12 uses the dereplicated inStrain reference.
- *File:* `figures/numbered/FigS4_instrain_derep.{png,pdf}`.
