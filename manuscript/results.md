# Strain2bScan — Results (draft)

> Organized around Strain2bScan's **two input modes** (see `docs/manuscript_organization.md`): a shared
> foundation (concept, 2bRAD-vs-16S motivation, core accuracy/robustness), then **Part I — native
> 2bRAD-M for low-biomass/high-host microbiomes** and **Part II — conventional metagenomes at community
> scale**. Figures are `figures/*`; source tables `results/*.tsv`; per-experiment design in
> `docs/benchmark_datasets.md`. Simulated benchmarks are 150 bp, error-free, closed-world; the mock and
> saliva chapters use real (2bRAD and shotgun) reads. This caveat is restated in the Discussion.

## Part 0 — Concept, motivation, and core accuracy

### Overview of Strain2bScan and its two data modes (Fig 1)

Strain2bScan resolves within-species strains from the sparse marker set released by type-IIB restriction
(2bRAD) digestion rather than from a full k-mer profile. For each species, reference genomes are digested
in silico with up to 16 type-IIB enzymes into single-copy 32–38 bp tags — roughly 1–2 % of the genome and
50–100× sparser than a full k-mer index — clustered into within-species groups at 0.95 Jaccard similarity
(MinHash-accelerated), and reduced to a cluster × marker database of species-core, cluster-specific and
strain-specific tags (Fig 1A). A sample is profiled by digesting its reads **once**, counting canonical
markers, gating on species-specific markers (Layer-1), and detecting and quantifying strains within each
present species from unique markers (Layer-2), returning a strain-level profile (Fig 1B). Crucially, the
pipeline accepts **two input modes**: in-silico digestion of conventional shotgun metagenomes, or **native
BcgI 2bRAD experimental libraries whose reads *are* the tags**. Because the tag lengths and recognition
patterns are identical to Fast2bRAD-M / `2bRADExtraction.pl`, the Fast2bRAD-M species layer and the
Strain2bScan strain layer operate in one shared tag space (Fig 1, band).

### 2bRAD tags carry strain-level signal that 16S cannot (Fig 2)

The premise of strain-level 2bRAD is that its reduced tag set, unlike the 16S gene, tracks genome-wide
divergence. We tested this directly across the **15 pathogenic/commensal species** of the motivation
panel, using complete/near-complete genomes only (CheckM completeness ≥ 97 %, contamination ≤ 5 %; 14–50
genomes per species). For every strain pair we computed between-strain distance in three spaces —
whole-genome (bottom-3000 21-mer MinHash), 2bRAD tags, and the 16S rRNA gene — and correlated the 2bRAD
and 16S distance vectors against the whole-genome vector (Spearman; 95 % CIs from genome subsampling).
2bRAD tracks the genome in every species (**median Spearman 0.94**, range 0.59–0.99), whereas **16S is
weak to useless (median 0.36)**, with several confidence intervals overlapping zero — *M. tuberculosis*
−0.04, *L. plantarum* 0.02, *Phocaeicola dorei* −0.13 (Fig 2). Every 2bRAD CI is high and non-overlapping
with its species' 16S CI. The per-pair rank–rank scatter (Fig S1) shows the mechanism: 2bRAD distances
order strain pairs the same way the whole genome does, while 16S collapses to a few discrete distances
shared by many unrelated pairs. **16S resolves species, not strains; the 2bRAD tag set carries
genome-wide strain signal at ~1 % of the sequencing.**

### Accurate and robust strain profiling (Fig 3–4)

On real reference panels with simulated strain mixtures, Strain2bScan achieved **precision 1.0** for
every species tested — *Cutibacterium acnes* (64 genomes → 16 clusters), *Staphylococcus aureus* (60 →
17) and *S. epidermidis* (60 → 10) — with high recall (0.75, 0.79, 1.0) and low abundance error
(Bray–Curtis 0.24, 0.33, 0.02; Fig 3A). Two further species (*Akkermansia muciniphila*, *Prevotella
copri*) reach precision 1.0 in the head-to-head of Fig 10. Precision held even for
low-diversity *S. epidermidis*, where 60 genomes collapse to 10 clusters: the profiler did not
over-detect near-identical strains. Recall, not precision, is the species-dependent axis and tracks
*genuine* intra-species diversity. Sensitivity matched the k-mer state of the art: profiling a single
*C. acnes* strain across a coverage ladder, Strain2bScan detected it at **0.5× per-strain coverage**,
missing it only at 0.1× — the **same detection onset as StrainScan** (Fig 3B) — so the sparse marker set
carries no low-depth penalty at realistic coverage.

Strain2bScan was also robust to reference-database imperfection (Fig 4). On nested *P. copri* panels of
40/80/112 genomes it held **precision 1.0** at every size (recall 1.0/0.95/1.0; 23/43/51 clusters), and
the 112-genome partition matched the clustering of StrainScan's own pre-built *P. copri* database — an
external check on the MinHash clustering (Fig 4A). Robustness to reference *quality* has a defined floor
(Fig 4B): degrading the truth strains' genomes held precision 1.0 down to **70 % completeness** (recall
falling from 0.75 to 0.19), with detection collapsing only at 50 % completeness / 10 % contamination /
400 contigs — motivating the built-in assembly-quality filter that flags low-quality references before
clustering.

## Part I — Native 2bRAD-M for low-biomass, high-host microbiomes

### A tunable enzyme set enables native BcgI 2bRAD operation (Fig 5)

Because each added type-IIB enzyme releases more tags, the enzyme set is a resolution/cost knob. On
*C. acnes*, **BcgI alone** (one enzyme) already gave precision 1.0 with recall 0.56; recall rose to 0.75
by **four enzymes** then plateaued, while strain-specific marker count and build time kept rising roughly
linearly (882 → 1 524 → 7 400 markers; 0.6 → 1.8 → 6.4 s from 1 → 4 → 14 enzymes; Fig 5). Cluster count
was unchanged (16 throughout): more enzymes add markers and recall, not resolution. A **~4-enzyme set is
the sweet spot** — and, decisively for this pillar, **single-enzyme BcgI operation is what lets native
BcgI 2bRAD-M experimental libraries be profiled directly**, since their reads are BcgI tags.

### 2bRAD holds strain recovery under host contamination where shotgun collapses (Fig 6)

We tested the low-biomass claim on the ATCC MSA-1002 20-strain even mock with **real DNA and known
truth**, profiled against a 62-species BcgI panel (20 mock + 42 decoys). Across a **host-contamination
series (0/90/99 % human DNA)**, native BcgI 2bRAD and in-silico-digested shotgun (WMS) of the same mock
were profiled with the same database. **Precision was 1.0 for both data types at every contamination
level** — zero false strain calls even at 99 % human DNA, with out-of-panel relatives correctly held in
the detected-not-resolvable tier. But recall diverged sharply: native 2bRAD retained **full 20/20 species
recall at 99 % host**, whereas shotgun fell to **12/20** (Fig 6A). The mechanism is marker yield (Fig 6B):
native 2bRAD delivers ~340–370k usable BcgI markers regardless of host fraction — the reduction happens at
the wet-lab step, before host swamps the library — while digesting an already host-dominated shotgun
library in silico yields only 96k → 53k → 33k markers as host rises. In a separate DNA-input titration,
native 2bRAD held precision 1.0 with full recall down to **0.1 ng** input (Fig S2). This is the
low-biomass niche, quantified on real data: the wet-lab 2bRAD reduction preserves the strain-informative
markers that host contamination otherwise destroys.

### Real saliva: individual-specific, temporally stable strain signatures (Fig 7)

Moving to real, error-containing, open-world data, we profiled native BcgI 2bRAD saliva from **8 subjects
sampled at 4 times of day (32 libraries)** against a 19-species oral-commensal reference panel. (Panel
choice matters: a generic pathogen panel with few genomes per species gave no signal, because real
strains map uniformly across arbitrary clusters; a genome-rich oral panel is required — see Discussion and
`docs/saliva_individual_discrimination.md`.) Strain-level community structure discriminated individuals
**better than species-level structure**: subject PERMANOVA **R² = 0.833 at strain level vs 0.822 at
species level** (both p = 2×10⁻⁴), and a leave-one-timepoint-out nearest-neighbour classifier identified
the host with **100 % accuracy from strain features vs 94 % from species features** (Fig 7A–C). Per
species, ***Rothia mucilaginosa* reached R² 0.921** — the strongest individual marker — with **17 of 18
testable species significant** (Fig 7D). Time of day had no detectable effect (R² ≈ 0.05, ns), and
within-subject strain distance (0.19) was far below between-subject distance (0.44; p = 5×10⁻²⁵): an
individual's salivary strain profile is a **stable, person-specific signature** resolvable in ~1 s per
sample. This reproduces a prior heavy-pipeline result (strain2bfunc: ~40 min, ~30 GB per sample) at a
fraction of the compute, on real 2bRAD reads.

### Native 2bRAD recovers low-abundance strains that shotgun misses (Fig 8)

Paired shotgun (WMS) of the same saliva samples let us validate the two modes against each other and
quantify the 2bRAD sensitivity advantage on real high-host material. **Every strain the shotgun mode
called was confirmed by native 2bRAD (65/65, 100 %)** — the in-silico digestion produces no strain calls
that native 2bRAD contradicts, establishing that the two input modes agree at the call level. Native
2bRAD additionally recovered **128–163 strains per sample that host-limited shotgun missed**, and these
2bRAD-only strains were **significantly lower-abundance** (median community relative abundance 0.0029 vs
0.0097 for shared strains; Mann–Whitney p = 1.2×10⁻²³; Fig 8). On host-dominated saliva, shotgun reaches
only the most abundant, ubiquitous strains; native 2bRAD confirms all of those *and* resolves the
low-abundance tail — the same enrichment advantage seen on the controlled mock (Fig 6), now on real
clinical material. (Strain2bScan runs directly on an oral clinical cohort — 15–17 species, ~2 s per
sample — as an application demonstration; Fig S4.)

## Part II — Conventional metagenomes at community scale

### Fast, light, and scalable to whole communities (Fig 9)

On conventional-metagenome (shotgun) input, Strain2bScan's efficiency is its decisive advantage. Per
sample it profiled the *C. acnes* panel in **0.86 s using 78 MB** versus **7.06 s and 828 MB** for
StrainScan — **~8× faster and ~11× lighter** (Fig 9A); the memory gap is structural, since StrainScan
counts the full k-mer set whereas Strain2bScan streams reads and keeps only the sparse markers. Both
database build and profiling parallelise cleanly with dependency-free threads (8.5× and 6.5× from 1 → 16
threads; Fig 9B). The decisive gain is at **community scale**: because a sample is digested **once** and
each additional species is a cheap hash-set match rather than a re-count, per-sample cost is essentially
**independent of the number of species**. On a **55-species** community Strain2bScan profiled each sample
in ~2.6–3.1 s and ran **121–146× faster** than running a per-species tool once per species (~132× at 100
samples, ~146× at ≥200 samples; Fig 9C) — the regime where full-k-mer tools, which have no multi-species
mode, must re-count every sample once per species (a projected ~10 h vs ~5 min at 100 samples).

### Matches or exceeds StrainScan on its own databases (Fig 10)

Head-to-head against StrainScan on **StrainScan's own reference sets**, both tools reached **precision
1.0** on *A. muciniphila* (40 genomes → 19 clusters) and *P. copri* (40 → 23), but Strain2bScan **matched
or exceeded StrainScan's recall** (0.93 vs 0.24; 0.94 vs 0.90) while running **~17–23× faster** and
**~15–24× lighter** (0.7–0.9 s / 71–85 MB vs 12.6–15.6 s / 1.2–1.7 GB; Fig 10). On near-clonal
*M. tuberculosis* (40 genomes → only 5 clusters), **StrainScan did not complete** (>3.3 h, >25 GB)
whereas Strain2bScan finished in 0.89 s; its low recall there (0.25) reflects the intrinsic near-clonality
of the species — a resolution limit, not a tool failure — and is reported at precision 1.0.

The validity of this in-silico (shotgun) mode is established independently in Part I: it recovers 20/20
species from the mock at 0 % host (Fig 6) and its saliva strain calls are a confirmed subset of the
native-2bRAD calls (Fig 8) — so the fast shotgun mode and the sensitive 2bRAD mode are two faces of one
validated method.
