# Strain2bScan — Methods (draft)

## Overview

Strain2bScan reimplements the two-layer StrainScan strategy — cluster near-identical strains,
then score samples on markers unique to a strain or cluster — but replaces the full k-mer set
with **2bRAD tags**, and is written in dependency-free Rust for speed and parallelism. The
pipeline is: (i) digest reference genomes and sample reads into 2bRAD-tag markers; (ii) build,
per species, a within-species cluster database annotated with unique markers; (iii) profile a
sample by detecting present clusters from their unique markers and estimating their abundance;
(iv) at the community level, digest each sample once and match it against every per-species
database, gated by a species-level Layer-1 check.

## 2bRAD tag extraction

Type-IIB restriction enzymes cut on both sides of their recognition site, releasing a
fixed-length fragment (the 2bRAD tag, 32–38 bp). Each of the 16 enzymes in the Fast2bRAD-M
table is modelled as a set of anchored sequence patterns; scanning every offset of a sequence
and testing the anchors reproduces the enzyme's digestion sites (forward and reverse patterns
allow scanning a single strand). Each tag is canonicalised (the lexicographically smaller of
the tag and its reverse complement) and hashed to a 64-bit integer marker (FNV-1a; genome and
sample tags use the same hash, so marker values are internally consistent). Two input modes
are supported: **BcgI only** for 2bRAD experimental libraries whose reads already are tags, or
a user-chosen enzyme set (`--enzyme all` for all 16) for in-silico digestion of conventional
shotgun reads, which enriches the marker set ~*n*-fold for *n* enzymes. For reference genomes
we retain only **single-copy** tags (occurring exactly once in the genome), following
StrainScan's and Fast2bRAD-M's use of single-copy markers for unbiased quantification.

## Reference database construction

**Within-species clustering.** For each species, genomes are grouped by single-linkage
hierarchical clustering at 0.95 marker-set similarity (0.05 distance), matching StrainScan's
`hclsMap_95`. Single-linkage at threshold τ is exactly the connected components of the graph
whose edges join genome pairs with Jaccard ≥ τ, computed with union-find. For panels of ≤96
genomes we use exact all-pairs Jaccard on the tag sets; above that we estimate Jaccard from
bottom-*k* MinHash sketches (*k* = 2000) of each genome's markers, which reduces the pairwise
cost from O(n²·m) to O(n²·k) with *k* ≪ *m* and yields partitions identical to exact on real
data (Results). Clusters are the finest reliable resolution unit: strains within one cluster
are too similar to separate from short reads.

**Containment clustering for uneven-completeness panels (`--containment`).** Jaccard penalises
incompleteness: an incomplete genome's markers are approximately a *subset* of a complete relative's,
so |A∩B|/|A∪B| falls below τ and the two spuriously split. The optional `--containment` mode instead
links on **max-containment**, |A∩B| / min(|A|,|B|), which stays ≈ 1 when one marker set is contained in
the other — the containment estimator used by Mash-screen and sourmash for uneven-completeness genomes.
It is exact for small panels; for large panels the intersection is estimated from the MinHash-sketch
Jaccard and the exact set sizes (|A∩B| = J·(|A|+|B|)/(1+J)), then divided by min(|A|,|B|). Because
max-containment ≥ Jaccard it merges at least as much, so it is opt-in (for reference sets of mixed
completeness) while the default stays Jaccard; the assembly-quality filter below is the
complementary first line of defence.

**Marker classification.** Within a species, each tag is labelled by its within-species
incidence — present in all clusters (*species-core*; detects the species, not strains), in one
cluster with ≥2 genomes (*cluster-specific*), in a single genome (*strain-specific*), or in
several but not all clusters (*shared-partial*). Cluster- and strain-specific tags are the
Layer-2 markers. Crucially these are derived from **all** tags of the species' genomes, not
from a pre-built species-unique database: species-unique markers (a genome compared against
genomes of *other* species) are computed for species detection and are orthogonal to
within-species strain structure. Each cluster's database is the union of its member genomes'
single-copy tags; a marker is *unique* to a cluster iff it occurs in exactly one cluster.

**Assembly-quality filtering.** Variable reference completeness biases Jaccard clustering
toward spurious splits: an incomplete genome's marker set is approximately a subset of its
complete twin's, so their Jaccard falls below 1 and they fail to cluster. Because CheckM is
not run in-line, two dependency-free proxies computed from data already at hand are used —
contig count (`--max-contigs`), and single-copy tag count relative to the conspecific median
(`--min-tag-fraction`, a completeness proxy). Genomes far below the median are always flagged;
they are removed only when a threshold is set.

## Layer-2 strain profiling

Sample reads are digested with the database's enzyme set (recorded in the database header) to
give per-marker counts. A cluster is called **present** iff at least *N* of its unique markers
are observed at count ≥ 2 (default *N* = 10, in tag units; the full-k-mer StrainScan floor
of ~1240 k-mers is inappropriate for the ~50–100× sparser tag set). Using only unique markers
makes detection immune to the shared-marker cross-talk that otherwise lets a greedy set-cover
over a large conspecific panel strip shared markers and starve true strains. Each present
cluster's **relative abundance** is estimated from the median sample count over its detected
unique markers — robust to repeat and contamination outliers, and less prone than a joint
regression over shared markers to mis-attributing signal between very similar co-present
strains (a non-negative Elastic Net solver is also provided). Calls are then filtered by a
minimum coverage fraction of their unique markers (`--min-coverage`, default 0.1; suppresses
spurious detection of large, similar clusters whose absolute unique-marker count clears the
floor at a tiny coverage fraction) and a minimum relative abundance (0.02), and renormalised.
When no cluster passes, Strain2bScan reports the species as detectable but not
strain-resolvable with the given enzyme set.

## Multi-species profiling and species selection

For community samples, the reads are digested **once** into a shared set of tag counts, matched
against every per-species cluster database in parallel; the per-species marginal cost is a
hash-set lookup rather than a re-count, so the total cost is independent of the number of species.

**Which species to strain-profile — the Layer-1 gate.** Strain markers are unique only *within* a
species, so a species absent from a sample can be spuriously hit by a present relative's shared
tags. Strain2bScan therefore decides per species from **absolute species-specific marker
evidence**, never relative abundance (which conflates community composition with sequencing
depth). Let *total* be the species-specific markers a species carries — tags unique to a single
species across the panel, the same tag space as the Fast2bRAD-M species layer — and *present* the
subset observed in the sample at count ≥ 2. The gate is

&nbsp;&nbsp;&nbsp;&nbsp;*resolve_gate* = max(*G*, ⌈*f* · *total*⌉), &nbsp; *detect_gate* = min(*d*, *resolve_gate*)

with an absolute floor *G* (`--min-species-markers`, default 200), a breadth fraction *f*
(`--min-species-marker-frac`, default 0) that scales the bar to each species' panel size, and a
low detection floor *d* (`--min-species-detect`, default 10). This yields three outcomes per
species: **strain-resolved** (*present* ≥ *resolve_gate*; Layer-2 runs), **detected but not
strain-resolvable** (*detect_gate* ≤ *present* < *resolve_gate*; reported at species level with its
observed marker breadth, no strain claim), or **absent**. The middle tier is the honest treatment
of a low-abundance species — present but too faint to support strain calls — rather than a binary
drop or an over-call. All inputs are computed by Strain2bScan from its own databases and a single
digest of the reads, so the gate needs no external abundance input; for open-world samples the
species presence call can instead be taken from an upstream Fast2bRAD-M step whose species
database is far broader than the strain panel.

**Gate calibration.** On the 55-species panel across normal and low (median 0.62×) depth, the
default floor gives species precision 1.0 at both depths, with leakage species correctly held in
the middle tier; at this panel the breadth term only trades recall, so *f* = 0 is optimal and is
the shipped default. The breadth term is scale insurance: when the floor is relaxed — or the panel
grows large enough for a fixed floor to be outrun by leakage — a small *f* (≈0.02) restores
precision to 1.0 at negligible recall cost, because it raises the bar in proportion to panel size,
where large-panel leakage concentrates (Results; `docs/gate_calibration.md`).

## Implementation

Strain2bScan is written in Rust with no third-party dependencies. Data-parallelism (genome
digestion, sketch construction, pairwise clustering, read digestion) uses scoped `std` threads
(`STRAIN2BSCAN_THREADS`; default = all cores). Databases are stored as sparse strain×marker
tables with an inverted unique-marker index and the enzyme set in the header.

## Benchmarking

**Datasets.** (i) A real *C. acnes* benchmark: a 64-genome reference panel (14 ground-truth
strains + 50 background, NCBI accessions pinned) and five paired-end mock samples from
MockMetagenomes4Benchmark (~100k read pairs each, ~12× total). (ii) A simulated multi-species
benchmark: 55 real species × ~4 strains (218 NCBI genomes) and 30 samples, each mixing strains
from twelve species at log-normal depth ≥1× (plus a low-depth variant, median 0.62×, used for
gate calibration). (iii) Cross-species mocks for *Staphylococcus
aureus* and *S. epidermidis* (60-genome panels each; 2–5 strains/sample, log-normal ≥1×,
matching the *C. acnes* design). (iv) A reference-degradation gradient in which the truth
strains' database genomes are degraded to completeness 100→50 % (with co-varying contamination
0→10 % and fragmentation), samples held fixed. All simulated reads are error-free 150 bp.

**Real-data and motivation datasets.** (v) *2bRAD-vs-16S motivation* (Fig 2): 15
pathogenic/commensal species, ~50 genomes each from NCBI accession lists (ENA FASTA), **restricted to
complete/near-complete assemblies** (CheckM completeness ≥ 97 %, contamination ≤ 5 %, assembly level
Complete Genome/Chromosome; `data/genome_qc_16s_panel.tsv`). Between-strain distance was computed in
three spaces — whole-genome (bottom-3000 canonical 21-mer MinHash), 2bRAD (Strain2bScan `build` BcgI
tags) and 16S (longest gene per genome via barrnap 0.9 + HMMER, 21-mer Jaccard) — all with the Mash
transform D(J) = −ln(2J/(1+J)); per species the 2bRAD and 16S pairwise vectors were correlated (Spearman)
against the whole-genome vector, with 95 % CIs from 500 genome subsamples. (vi) *DNA mock,
host-contamination series* (Fig 6, Fig S1): ATCC MSA-1002 20-strain even mock — native BcgI 2bRAD (SRA
PRJNA1131785, 0/90/99 % human DNA; Figshare 12272360 DNA-input titration) and shotgun WMS of the same
mock — profiled with `multi-profile --enzyme BcgI` against a reconstructed 62-species BcgI panel (20
mock + 42 decoys), scored against the 20-species truth. (vii) *Real saliva* (Fig 7, Fig 8): native BcgI
2bRAD (and paired shotgun WMS) saliva from PRJNA1131785, 8 subjects × 4 within-day timepoints, profiled
against a 19-species oral-commensal panel (up to 25 genomes/species). Strain- and species-level relative
abundances → Bray–Curtis → PERMANOVA (adonis, subject/timepoint factors) and leave-one-timepoint-out
1-NN host classification; shotgun R1 (in-silico BcgI) compared to native 2bRAD calls per sample. Full
per-dataset procedures and accessions are in `docs/` (`motivation_16s.md`, `mock_hostcontam.md`,
`saliva_individual_discrimination.md`, `saliva_temporal_ml.md`, `saliva_concordance.md`).

**Systematic head-to-head on a 15-species simulated benchmark (Fig 11, Table 1–3).** A common
benchmark was built from a fixed pool of 15 pathogenic/commensal species (15–50 complete/near-complete
NCBI genomes each; `figure_raw_data/sim_pool_manifest.tsv`). *Single-species* samples were generated for
every species as 2/3/5 co-present strains drawn either from the same or from different 0.95 clusters, at
per-strain coverages 0.5/1/3/5/10× with uneven abundance ratios (following StrainScan's simulation
design), 5 replicates per cell — 2 025 samples. *Multi-species* samples mixed ~18 co-present species
(one to a few strains each) across three community depth gradients — 60 samples. Reads were simulated
with ART (`art_illumina`, 150 bp paired-end) from the truth genomes; truth tables record each strain's
species, genome accession and 0.95-cluster assignment.

Both tools built their databases from the **same genome pool** and profiled the **same reads**.
Strain2bScan databases were built with `cluster --enzyme all --similarity 0.95` and profiled with
`profile` / `multi-profile --enzyme all` (reads decompressed, R1+R2 concatenated). StrainScan (v1.0.14,
bioconda) is Linux-x86-only — it ships `dashing_s128` and `jellyfish-linux` ELF binaries, a Python-3.7
`.so`, and an R reclustering step — so it was run inside a Docker `linux/amd64` container (QEMU emulation
on Apple Silicon; `strainscan_build`, then `strainscan -i R1 -j R2 -d DB`). Because the two tools cluster
genomes independently, **each tool was scored in its own cluster space**: predicted clusters were compared
against the truth strains mapped into that tool's clusters — for Strain2bScan via the truth `cluster`
column, for StrainScan via its `Cluster_Result/hclsMap_95_recls.txt` (report `Cluster_ID` = `C`+cluster
id) — and precision/recall/F1 computed over the cluster sets per sample. Strain2bScan profiled all 2 025 +
60 samples; StrainScan profiled a matched subset (different-cluster mixtures, k = 2/3/5, one replicate, all
depths; near-clonal *M. tuberculosis* via its same-cluster samples) — 204 depth-matched paired
single-species samples across 14 species, plus 4 multi-species samples per depth. StrainScan has no
multi-species mode, so each community sample was profiled once per species database and the per-sample
cost taken as the sum of wall-clock over species (peak RSS as the maximum).

*Timing.* Strain2bScan build and profile times are native (arm64). To compare profiling speed free of the
emulation confound, a `linux/amd64` Strain2bScan binary (zero-dependency `cargo build`) was run **inside
the same container** on the same subset, giving the same-environment ratio of Fig 11E/Table 2; StrainScan
build times are reported in the emulated environment (an upper bound). DB build for *K. pneumoniae*
(47 genomes × 5.5 Mb) did not complete under StrainScan (killed at a 100-min cap; still in the k-mer-matrix
step past 1 h 40 min on a longer retry), and that species is omitted from the paired accuracy set. The
container's VM memory was raised to 56 GB because StrainScan's build peaks at ~28 GB (vs ≤0.4 GB for
Strain2bScan). Scripts: `scripts/plot_sim_headtohead.py` and the drivers under `scratchpad/eval/`
(`run_s2b_{single,multi}.py`, `run_strainscan_{single,multi}.py`, `run_s2b_emulated_single.py`,
`analyze_headtohead.py`); raw per-sample tables in `figure_raw_data/sim_headtohead/`.

**Comparison to StrainScan (curated-DB and per-sample benchmarks).** In addition to the common benchmark
above, StrainScan (v1.0) was run on its **own** reference databases (Fig 10) and on the same *C. acnes*
per-sample profiling comparison (Fig 9A), using its low-depth modes for the depth series.

**Metrics.** Detection precision, recall and F1 at a 0.01 presence threshold; abundance error
by L1 distance and Bray–Curtis dissimilarity over the union of predicted and true labels,
evaluated at cluster resolution (ground-truth strains mapped to their clusters). Wall-clock
time and peak resident set size were measured with `/usr/bin/time -l` on a 16-core Apple
silicon machine. All scripts, pinned accession lists, result tables and figure code are in the
Strain2bScan-paper repository; every figure is regenerable with `make figures`.
