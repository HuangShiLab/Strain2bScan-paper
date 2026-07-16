# Strain2bScan: strain-level metagenomic profiling on 2bRAD-reduced markers for low-biomass microbiomes and community-scale cohorts

*Assembled manuscript. Authors, affiliations, funding and data-availability accessions to be completed.*

---

## Abstract

**Background.** Strain-level variation drives clinically important microbial phenotypes, but 16S
surveys cannot resolve strains and shotgun k-mer methods, while accurate, are held back from two
increasingly important settings: they are computationally heavy across the many-species, many-sample
cohorts of modern studies, and they fail on low-biomass, host-contaminated samples (saliva, tumour/FFPE,
skin) where most DNA is human. Reduced-representation 2bRAD sequencing samples a sparse, reproducible
subset of the genome, is robust to low input and host contamination, and has been used for species-level
profiling — but not for strain resolution.

**Results.** We present **Strain2bScan**, a dependency-free Rust strain profiler that ports the
StrainScan resolution framework (within-species clustering plus unique-marker scoring) onto **2bRAD
markers** — the 32–38 bp tags released by type-IIB restriction digestion — and uniquely accepts **two
input modes**. As a shared foundation it is accurate (**precision 1.0** across species, high recall, low
abundance error), detects strains to **0.5× coverage matching StrainScan**, achieves precision 1.0 with
complete references across all 15 species (declining as references degrade — complete/near-complete
genomes matter), and — unlike 16S, whose between-strain distances
are uncorrelated with genome divergence (median Spearman 0.36 across 15 species) — its tags track
whole-genome strain distance (median 0.94). **(1) Native 2bRAD-M for low-biomass microbiomes:** on the
ATCC MSA-1002 mock, native 2bRAD holds precision 1.0 and **20/20 strain recall at 99 % human DNA** where
in-silico-digested shotgun drops to 12/20; on real saliva (8 subjects × 4 timepoints) strain-level
profiles discriminate individuals better than species-level (PERMANOVA **R² 0.83 > 0.82**; leave-one-
timepoint-out host-ID **100 %**), are temporally stable, and recover **128–163 low-abundance strains per
sample that host-limited shotgun misses**. **(2) Conventional metagenomes at scale:** the sample is
digested once and matched against every species database, so per-sample cost is independent of the
number of species — **~8× faster and ~11× lighter** than StrainScan per sample and **~130–146× faster on
a 55-species community**. On a common 15-species simulated benchmark (both tools building databases from
the same genomes and profiling the same reads), Strain2bScan **matched StrainScan's precision (1.0) and
exceeded its recall** (median 0.80 vs 0.67, full recall by 3× vs 10×) while building databases
**249–614× faster and 43–138× lighter** and profiling **4–105× faster** — and it completed the one
species (*Klebsiella pneumoniae*) StrainScan could not build. The two modes agree, so one tool spans both
regimes.

**Conclusions.** Strain2bScan delivers accurate, genome-resolved strain profiling from a sparse marker
set: it uniquely enables strain-level analysis of low-biomass, high-host 2bRAD-M data, and scales
conventional-metagenome strain profiling to communities of many species across many samples.

**Availability.** Rust source: https://github.com/HuangShiLab/Strain2bScan · reproducible benchmarks and
figures: https://github.com/HuangShiLab/Strain2bScan-paper.

---

## Introduction

Microbial strains of the same species can differ sharply in phenotype — virulence, antibiotic
resistance, metabolic capacity and host interaction — so resolving *which* strains are present, and at
what abundance, is central to microbiome and clinical metagenomics. The 16S rRNA gene, the workhorse of
community surveys, cannot see this variation: its between-strain distances are essentially uncorrelated
with genome-wide divergence, so 16S resolves species but not strains (below, and Fig 2). Shotgun
metagenomes do carry strain information, and a family of tools now recovers it — reference k-mer methods
such as StrainScan and StrainGE, marker-gene methods such as StrainPhlAn, and sketch-based abundance
estimators such as sylph — making strain-level profiling routine for individual samples and species.

Two obstacles keep shotgun-based strain profiling from two increasingly important settings.

**Scale.** Modern studies span hundreds of samples and aim to resolve strains across the
dozens-to-hundreds of species in a community. Full k-mer methods index and query the entire k-mer
content of every reference genome; per-sample cost is dominated by counting the sample's k-mers, and is
paid again for every species database queried. For *S* species across *N* samples this scales as
*N × S × (k-mer count + search)* — hours-to-days of compute and hundreds of megabytes to gigabytes of
memory — a practical ceiling for population-scale, multi-species strain surveillance.

**Low biomass and host contamination.** Many of the most interesting clinical niches — saliva and other
oral sites, tumour and FFPE tissue, skin — yield little microbial DNA against a large human background
(often ≥90–99 % host). Shotgun sequencing spends most of its reads on the host, so the informative
bacterial fraction, and with it the rarer strains, is lost; strain resolution collapses exactly where it
would be most valuable.

Reduced-representation sequencing addresses both. Type-IIB restriction ("2bRAD") digestion releases a
sparse, reproducible set of fixed-length tags — roughly 1–2 % of the genome — and, applied as a wet-lab
protocol (2bRAD-M / Fast2bRAD-M), it enriches these informative markers *before* host DNA can swamp the
library, giving accurate profiles from picogram inputs, heavily host-contaminated samples and degraded
material. To date, however, 2bRAD has been used only for **species**-level profiling. Yet strain-resolution
methods do not need the whole genome — they score samples on *unique* markers that distinguish a strain or
a cluster of near-identical strains — and a 2bRAD tag set is exactly this kind of low-redundancy,
taxonomically informative marker set, 50–100× smaller than the full k-mer set.

We present **Strain2bScan**, which ports the StrainScan resolution framework — within-species clustering
into a search structure, followed by unique-marker scoring and abundance estimation — onto 2bRAD tags, in
dependency-free Rust. Its tag lengths and recognition patterns match Fast2bRAD-M / `2bRADExtraction.pl`
exactly, so the tags are interoperable with the Fast2bRAD-M species layer, and Strain2bScan uniquely
accepts **two input modes** that map onto the two obstacles above:

1. **Native BcgI 2bRAD experimental libraries** — enabling, for the first time, strain-level analysis of
   low-biomass and high-host microbiomes. Because the reduction happens at the bench, Strain2bScan holds
   precision 1.0 and full strain recall at 99 % host DNA where in-silico-digested shotgun loses most of
   its strains, and on real saliva it resolves individual-specific, temporally stable strain signatures —
   recovering low-abundance strains that host-limited shotgun cannot reach.

2. **In-silico digestion of conventional shotgun metagenomes** — enabling community-scale strain
   profiling. The sample is digested **once** and matched against every per-species database, so
   per-sample cost is independent of the number of species and linear in the number of samples: ~8×
   faster and ~11× lighter per sample than StrainScan, and ~130–146× faster on a 55-species community,
   while matching StrainScan's precision, its 0.5× detection onset and its recall on its own databases —
   and completing near-clonal *Mycobacterium tuberculosis* in ~1 s where StrainScan does not finish.

The two modes are shown to agree (in-silico and native digestion recover the same strains), so a single
tool, on one 2bRAD tag space, spans both the low-biomass clinical regime and the cohort-scale regime.

## Results

### Part 0 — Concept, motivation, and core accuracy

#### Overview of Strain2bScan and its two data modes (Fig 1)

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

#### 2bRAD tags carry strain-level signal that 16S cannot (Fig 2)

The premise of strain-level 2bRAD is that its reduced tag set, unlike the 16S gene, tracks genome-wide
divergence. We tested this directly across the **15 pathogenic/commensal species** of the motivation
panel, using complete/near-complete genomes only (CheckM completeness ≥ 97 %, contamination ≤ 5 %; 14–50
genomes per species). For every strain pair we computed between-strain distance in three spaces —
whole-genome (bottom-3000 21-mer MinHash), 2bRAD tags, and the 16S rRNA gene — and correlated the 2bRAD
and 16S distance vectors against the whole-genome vector (Spearman; 95 % CIs from genome subsampling).
2bRAD tracks the genome in every species (**median Spearman 0.94**, range 0.59–0.99), whereas **16S is
weak to useless (median 0.36)**, with several confidence intervals overlapping zero — *M. tuberculosis*
−0.04, *L. plantarum* 0.02, *Phocaeicola dorei* −0.13 (Fig 2). Every 2bRAD CI is high and non-overlapping
with its species' 16S CI. The per-pair rank–rank scatter (Fig 2B) shows the mechanism: 2bRAD distances
order strain pairs the same way the whole genome does, while 16S collapses to a few discrete distances
shared by many unrelated pairs. **16S resolves species, not strains; the 2bRAD tag set carries
genome-wide strain signal at ~1 % of the sequencing.**

#### Accurate and robust strain profiling (Fig 3–4)

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

Reference-genome **completeness** limits strain-identification accuracy — and Strain2bScan provides a
clustering mode that addresses it, tested across **all 15 species** (Fig 4). Holding the sample reads
fixed and progressively degrading the truth strains' reference genomes (completeness 100 → 50 %, with
co-varying contamination and fragmentation), default **Jaccard** clustering loses accuracy as
completeness falls — across the 14 **resolvable** species, median precision 1.0 → 0.84 (90 %) → 0.71
(70 %) and recall 0.96 → 0.80 (90 %) → 0.74 (70 %) — because an incomplete genome's markers are (approx.)
a **subset** of its complete relative's, dropping their Jaccard below the 0.95 cut and spuriously
splitting them into new singleton clusters. The **`--containment` clustering mode** dissolves this:
max-containment = |A∩B| / min(|A|,|B|) stays ≈ 1 when one genome is a subset of another (the estimator
Mash-screen/sourmash use), so incomplete genomes cluster with their complete relatives instead of
fragmenting. With `--containment`, median **precision rises to 0.98 / 0.92 and recall to 0.95 / 0.92 at
95 % / 90 % completeness** and stays ahead of Jaccard down to ~80 %, converging only at ≤ 70 % where the
references are genuinely low-quality (MIMAG-low; heavy contamination + hundreds of contigs). It also
**removes the near-clonal *M. tuberculosis* artifact**: its single 0.95 cluster stays intact, so recall
holds at 1.0 to 90 % under containment versus collapsing to ≈ 0.05 under Jaccard the moment references
degrade (Fig 4 inset) — confirming that collapse was cluster fragmentation, not a true completeness
effect. Because max-containment merges slightly more aggressively (it can coarsen clusters on
already-complete panels), it is an **opt-in mode for reference panels of uneven completeness**; the
default remains Jaccard, complemented by the built-in assembly-quality filter (`--min-tag-fraction` /
`--max-contigs`) that drops low-quality genomes before clustering and the complete-genome practice used
elsewhere (e.g. Fig 2).

### Part I — Native 2bRAD-M for low-biomass, high-host microbiomes

#### A tunable enzyme set enables native BcgI 2bRAD operation (Fig 5)

Because each added type-IIB enzyme releases more tags, the enzyme set is a resolution/cost knob. On
*C. acnes*, **BcgI alone** (one enzyme) already gave precision 1.0 with recall 0.56; recall rose to 0.75
by **four enzymes** then plateaued, while strain-specific marker count and build time kept rising roughly
linearly (882 → 1 524 → 7 400 markers; 0.6 → 1.8 → 6.4 s from 1 → 4 → 14 enzymes; Fig 5). Cluster count
was unchanged (16 throughout): more enzymes add markers and recall, not resolution. A **~4-enzyme set is
the sweet spot** — and, decisively for this pillar, **single-enzyme BcgI operation is what lets native
BcgI 2bRAD-M experimental libraries be profiled directly**, since their reads are BcgI tags.

#### 2bRAD holds strain recovery under host contamination where shotgun collapses (Fig 6)

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
native 2bRAD held precision 1.0 with full recall down to **0.1 ng** input (Fig S1). This is the
low-biomass niche, quantified on real data: the wet-lab 2bRAD reduction preserves the strain-informative
markers that host contamination otherwise destroys.

#### Real saliva: individual-specific, temporally stable strain signatures (Fig 7)

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

#### Native 2bRAD recovers low-abundance strains that shotgun misses (Fig 8)

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
sample — as an application demonstration; Table S3.)

### Part II — Conventional metagenomes at community scale

#### Fast, light, and scalable to whole communities (Fig 9)

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

#### Matches or exceeds StrainScan on its own databases (Fig 10)

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

#### Systematic head-to-head on a 15-species simulated benchmark (Fig 11, Table 1)

Where Fig 10 compares the tools on StrainScan's own curated databases, we also ran both tools **end to
end on a common, controlled benchmark**: the same 15-species reference pool, the same simulated reads,
and databases built by each tool from the same genomes. Single-species communities span 2/3/5 co-present
strains (from the same or different clusters) across a **0.5–10× per-strain coverage ladder** (2 025
samples); multi-species communities span **~18 co-present species** across three depth gradients (60
samples). Because the two tools cluster genomes independently, each is scored **in its own cluster space**
(the honest resolution unit for short reads; Methods). StrainScan is Linux-x86-only and ran under
`linux/amd64` emulation on Apple Silicon, so its wall-clock is emulation-inflated; to remove that
confound from the speed comparison we also ran a `linux/amd64` build of Strain2bScan **in the same
container**, giving a same-environment ratio (Fig 11E).

**Accuracy: equal precision, higher recall (Fig 11A–C, Table 1).** Across 14 resolvable species and 204
depth-matched samples, **both tools held median precision 1.0** — but only Strain2bScan held precision
1.0 in *every* species, whereas StrainScan dropped to 0.80–0.83 in four (*A. muciniphila*, *C. difficile*,
*P. dorei*, *S. enterica*), making occasional false cluster calls. Strain2bScan also **recovered
low-coverage strains at shallower depth**: median recall reached 1.0 by **3× coverage** versus **10×** for
StrainScan, and overall median recall/F1 were **0.80 / 0.89 (Strain2bScan) vs 0.67 / 0.75 (StrainScan)**.
Near-clonal *M. tuberculosis*, whose 29 genomes collapse to a single 0.95 cluster, was resolved correctly
by both (precision and recall 1.0 at the cluster level — the honest claim for a species short reads cannot
sub-resolve). On the multi-species communities the two tools were **comparable in accuracy** (StrainScan
marginally higher recall on the deepest communities, Strain2bScan higher precision at low depth; Table 1).

**Cost: two-to-three orders of magnitude cheaper (Fig 11D–F, Table 1).** The gap is largest in **database
construction**: Strain2bScan built each species' database in **0.7–5.1 s using 0.1–0.4 GB**, whereas
StrainScan required **5–43 min and 8–28 GB** — **249–614× faster and 43–138× lighter** per species
(e.g. *E. coli*: 4.3 s / 0.34 GB vs 43 min / 28 GB). The single heaviest species, *K. pneumoniae*
(47 genomes × 5.5 Mb), **StrainScan failed to build at all** (killed after >100 min, and stalled again
past 1 h 40 min in the k-mer-matrix step on a longer retry), while Strain2bScan built it in 5.1 s;
it is therefore excluded from the paired accuracy set. **Per-sample profiling** in the same emulated
environment was **4–33× faster and 5–39× lighter** for Strain2bScan (0.16–1.7 s / 20–160 MB vs
3.5–6.9 s / ~830 MB; Fig 11E). The advantage compounds at community scale (Fig 11F): because StrainScan
has no multi-species mode, each community must be profiled once per species, so its per-sample cost was
the **sum of 14 runs — 100–398 s — versus 1–9 s for Strain2bScan's single digest-once pass (46–105×
faster)**, consistent with the structural *S*-fold argument of Fig 9C. In short, on a common benchmark
Strain2bScan **matches or exceeds StrainScan's strain-identification accuracy while building databases
2–3 orders of magnitude faster and lighter and profiling 1–2 orders faster** — and completes the one
species StrainScan could not build.

## Discussion

Strain2bScan reframes strain-level profiling around a sparse, reproducible marker set, and does so as a
single engine serving **two input modes**: native 2bRAD-M libraries and in-silico digestion of
conventional shotgun. By porting StrainScan's clustering-plus-unique-marker framework onto 2bRAD tags in
Rust, it turns a ~50–100× smaller marker database into large, structural gains **without sacrificing
accuracy** — across every species tested it held precision 1.0, matched StrainScan's detection onset at
0.5× coverage (Fig 3), and matched or exceeded StrainScan's recall both on StrainScan's own databases
(Fig 10) and on a common 15-species simulated benchmark where the two tools built databases from the same
genomes and profiled the same reads (median recall 0.80 vs 0.67 at equal precision, while building
databases 249–614× faster and profiling 4–105× faster; Fig 11). The foundation for a strain-level *2bRAD* method is that its tags, unlike the 16S gene, carry
genome-wide strain signal: across 15 species 2bRAD between-strain distances tracked the whole genome
(median Spearman 0.94) while 16S did not (0.36), several 16S intervals overlapping zero (Fig 2A,B).
16S resolves species; 2bRAD tags resolve strains.

**Native 2bRAD for low-biomass, high-host microbiomes.** The distinctive advantage of the native-2bRAD
mode is that the reduction happens at the bench, before host DNA can swamp the library. On the ATCC
MSA-1002 mock at 99 % human DNA, native 2bRAD held precision 1.0 and full 20/20 strain recall where
in-silico-digested shotgun of the same material dropped to 12/20, because 2bRAD preserves ~10× more
usable markers under host load (Fig 6). On real saliva this translated into biology: strain-level
profiles discriminated individuals better than species-level profiles (PERMANOVA R² 0.83 vs 0.82;
leave-one-timepoint-out host-ID 100 % vs 94 %), were temporally stable within subject, and — validated
against paired shotgun of the same samples — recovered 128–163 low-abundance strains per sample that
host-limited shotgun could not reach, while calling nothing the shotgun mode contradicted (Fig 7, Fig 8).
This is the setting where strain resolution matters most (oral, tumour/FFPE, skin) and where shotgun is
weakest; a wet-lab reduction that concentrates sequencing on informative tags is the natural fit, and
Strain2bScan is, to our knowledge, the first tool to make native 2bRAD strain-resolved.

**Conventional metagenomes at community scale.** In the shotgun mode, the decisive gains are per-sample
speed and memory (~8× faster, ~11× lighter than StrainScan; Fig 9A) and, above all, scaling in the
number of species: because a sample is digested once and matched against every per-species database, the
marginal cost of an additional species is a hash lookup rather than a re-count. For *S* species over *N*
samples the cost is ≈ *N × (digest + S·ε)* versus ≈ *N × S ×* (k-mer count + search) for a full k-mer
tool run per species — an *S*-fold structural advantage that, measured on a 55-species community, reached
**~132× at 100 samples and ~146× beyond** (minutes versus projected hours) and grows with community
richness (Fig 9C). The two modes are not separate tools: the shotgun mode is validated by the mock
(20/20 at 0 % host, Fig 6) and by the saliva concordance (its calls are a confirmed subset of the
native-2bRAD calls, Fig 8), so the fast shotgun mode and the sensitive 2bRAD mode are two faces of one
method.

**Cluster resolution is the honest unit of strain analysis.** Short reads cannot separate strains that
share almost all of their sequence, so both StrainScan and Strain2bScan resolve to *clusters* of
near-identical strains and we evaluate at that resolution. The clusters-to-genomes ratio at the 0.95 cut
is a per-species property — for diverse species (*C. acnes*) clusters are essentially single strains,
whereas for near-clonal *M. tuberculosis* (5 clusters from 40 genomes) the cluster is the honest level of
claim. Resolving to clusters does not cost precision (it stayed 1.0 even for low-diversity
*S. epidermidis*); what varies with diversity is recall.

**Design choices.** Occurrence-based uniqueness — a marker is cluster-unique only if absent, at any copy
number, from every other cluster — eliminates false-unique markers from single-copy filtering and keeps
precision 1.0 in similar-strain species. MinHash-sketch clustering scales database construction while
producing partitions identical to exact Jaccard and to StrainScan's own pre-built *P. copri* clustering
(112 → 51 clusters; `results/panelsize_prevotella.tsv`). Reference incompleteness is the one factor that
genuinely degrades strain identification under Jaccard — an incomplete genome's markers are a subset of a
complete relative's, so the two fall below the 0.95 similarity cut and split (Fig 4). We address this two
ways: the built-in assembly-quality filter drops low-quality genomes before clustering
(`--min-tag-fraction`/`--max-contigs`), and the optional **`--containment` clustering mode**
(max-containment) keeps subset genomes with their complete relatives, restoring precision/recall down to
~80 % completeness and removing the near-clonal *M. tuberculosis* cluster-fragmentation artifact (Fig 4).
The same completeness concern motivated restricting the Fig 2 motivation panel to complete/near-complete
genomes, where the 2bRAD-over-16S advantage is unchanged (median 0.90→0.94), confirming it is not a
draft-assembly artifact. The enzyme count is an explicit resolution/cost control
(~4 enzymes captures most recall; Fig 5), and single-enzyme BcgI operation is what enables native
2bRAD-M libraries.

**Limitations and future work.** (i) **Near-clonality** caps recall on species where short reads cannot
distinguish strains (*M. tuberculosis*); this limit is intrinsic and shared by all short-read tools —
and where it bites hardest StrainScan failed to complete (>3.3 h, >25 GB) while Strain2bScan finished in
~1 s at precision 1.0 (Fig 10). Layering a within-cluster overlap/regression step on top of
occurrence-based uniqueness is the clearest algorithmic target for raising near-clonal recall.
(ii) **Reference panels must be niche-appropriate and genome-rich** for real communities: a generic
pathogen panel with few genomes per species gave no saliva signal because real strains map uniformly
across arbitrary clusters, whereas a genome-rich oral panel recovered the full individual-discrimination
signal — panel design is a real determinant of strain-level performance on open-world data.
(iii) **Host-limited shotgun** cannot reach the low-abundance strain tail on high-host samples; this is a
property of the input, not the tool, and is precisely the gap the native-2bRAD mode fills.
(iv) **Reference incompleteness is improved but not fully solved.** The `--containment` clustering mode
recovers accuracy when an incomplete genome has a complete relative in the panel (Fig 4), but it cannot
recover markers that are simply *absent* — a strain whose only reference is a partial assembly — nor undo
contamination that injects foreign tags, so it converges with Jaccard below ~70 % completeness; and
because it merges more aggressively it trades a little resolution on complete panels (hence opt-in, not
the default). Making strain identification *more resistant* to incomplete references is a clear direction:
a **completeness-aware detection gate** (scale the unique-marker floor by each genome's estimated
completeness, or gate on a *fraction* of a cluster's available markers rather than an absolute count —
analogous to the Layer-1 breadth term) so incomplete strains are not gated out; a
**best-quality-representative** marker set per cluster (define the cluster's markers from its most-complete
member); **upstream completeness/contamination estimation and decontamination** (CheckM2 / GUNC) feeding
the quality filter; and **pangenome-based imputation** of missing markers from complete conspecifics. The
irreducible case — a strain represented only by a low-completeness, contaminated genome — is a data limit
no clustering can overcome. (v) Aside
from the mock and saliva chapters, the accuracy benchmarks are simulated, error-free and closed-world.
Extensions include oral-cancer case/control analysis (needs the study's sample labels), FFPE and
degraded material, deeper multi-tool comparison (sylph, StrainGE), and completing the Fast2bRAD-M species
layer so species and strain calls come from one 2bRAD digest.

**Conclusion.** Reduced-representation 2bRAD markers, combined with a StrainScan-style resolution
framework and a fast Rust implementation, make accurate strain-level profiling practical at a fraction of
the compute and memory of full-k-mer methods. Uniquely, the same tool spans two regimes: native 2bRAD-M
for strain-level analysis of low-biomass, high-host microbiomes, and in-silico-digested shotgun for
strain profiling across communities of many species and many samples.

## Methods

### Overview

Strain2bScan reimplements the two-layer StrainScan strategy — cluster near-identical strains,
then score samples on markers unique to a strain or cluster — but replaces the full k-mer set
with **2bRAD tags**, and is written in dependency-free Rust for speed and parallelism. The
pipeline is: (i) digest reference genomes and sample reads into 2bRAD-tag markers; (ii) build,
per species, a within-species cluster database annotated with unique markers; (iii) profile a
sample by detecting present clusters from their unique markers and estimating their abundance;
(iv) at the community level, digest each sample once and match it against every per-species
database, gated by a species-level Layer-1 check.

### 2bRAD tag extraction

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

### Reference database construction

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

### Layer-2 strain profiling

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

### Multi-species profiling and species selection

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

### Implementation

Strain2bScan is written in Rust with no third-party dependencies. Data-parallelism (genome
digestion, sketch construction, pairwise clustering, read digestion) uses scoped `std` threads
(`STRAIN2BSCAN_THREADS`; default = all cores). Databases are stored as sparse strain×marker
tables with an inverted unique-marker index and the enzyme set in the header.

### Benchmarking

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

## Figure legends

### Main figures

**Figure 1. Strain2bScan overview and the two input modes.**
Pipeline schematic. (A) Reference construction: type-IIB (2bRAD) digestion of reference genomes into
single-copy 32–38 bp tags → within-species clustering at 0.95 Jaccard (MinHash-accelerated) → a
cluster × marker database of species-core, cluster-specific and strain-specific markers. (B) Profiling:
a sample is digested once into canonical markers, gated on species-specific markers (Layer-1), and
strains are detected and quantified within each present species from unique markers (Layer-2). The two
input modes — in-silico digestion of conventional shotgun, or native BcgI 2bRAD libraries whose reads
are the tags — enter the same tag space, shared with the Fast2bRAD-M species layer.

**Figure 2. 2bRAD tags track genome-wide strain distance; 16S does not.** **(A)** Per-species Spearman correlation between
whole-genome (bottom-3000 21-mer MinHash) between-strain distance and the corresponding 2bRAD-tag (blue)
or 16S rRNA (red) distance, across 15 species (complete/near-complete genomes only, CheckM ≥97 %/≤5 %;
14–50 genomes per species, shown at right); error bars are 95 % CIs over genome subsamples; genus and
species on separate label lines; 2bRAD median 0.94 vs 16S median 0.36. **(B)** 3×5 matrix of per-species
rank–rank scatters — rank of each strain pair's whole-genome distance (x) vs its 2bRAD (blue) or 16S (red)
distance (y): 2bRAD hugs the diagonal in every species, 16S forms flat rank-bands (most extreme in
*P. dorei*, *M. tuberculosis*).

**Figure 3. Accurate strain profiling and depth sensitivity.** (A)
— precision/recall/abundance error (Bray–Curtis) on real reference panels with simulated mixtures for
*C. acnes*, *S. aureus*, *S. epidermidis* (precision 1.0 throughout). (B)
— detection of a single *C. acnes* strain across a coverage ladder; onset at 0.5×, matching StrainScan.

**Figure 4. Reference-genome completeness controls strain-identification accuracy (all 15 species).** Sample reads held fixed while the truth strains' reference genomes are
degraded (completeness 100→50 %, with co-varying contamination + fragmentation), DB rebuilt and re-profiled.
Problem → mechanism → solution. **(A)** schematic (): an incomplete
genome's tags are a subset of a complete relative's, so Jaccard (6/10 = 0.60) drops below the 0.95 cut and
splits them — the shared tags land in two clusters and are demoted to non-discriminating *SharedPartial*;
max-containment (6/6 = 1.0) merges them into one cluster whose marker set is the **union** of members, so
the discriminating *ClusterSpecific* tags (and the complete member's markers) are preserved. **(B)** median
precision, **(C)** median recall vs reference completeness (14 resolvable
species; faint lines = per-species containment), comparing **default Jaccard** clustering (grey dashed)
with the **`--containment`** mode (solid, shaded gap). Degrading references collapses Jaccard accuracy
(precision 1.0→0.84→0.71, recall 0.96→0.80→0.74 at 100→90→70 %) because incomplete genomes split from
complete relatives; `--containment` (max-containment) keeps them together, restoring precision to
0.98/0.92 and recall to 0.95/0.92 at 95/90 %, converging with Jaccard only at ≤70 % (genuinely
low-quality). Inset (C): near-clonal ***M. tuberculosis*** recall — Jaccard collapses to ≈0.05 on any
degradation (single cluster shatters), containment holds 1.0 to 90 % (artifact fixed).

**Figure 5. The 2bRAD enzyme set is a resolution/cost knob enabling native BcgI operation.**
 — on *C. acnes*, precision/recall, strain-specific marker count and
build time vs number of enzymes (1→14). BcgI-alone gives precision 1.0; ~4 enzymes is the sweet spot;
cluster count is invariant (16). Single-enzyme BcgI is what enables native 2bRAD-M libraries.

**Figure 6. DNA mock, host-contamination series: native 2bRAD holds strain recovery where shotgun collapses.**
 — ATCC MSA-1002 mock across 0/90/99 % human DNA, native BcgI 2bRAD
vs in-silico-digested shotgun (WMS), same 62-species BcgI DB. (A) Species recall: precision 1.0 for both
at every level; 2bRAD 20/20 at 99 % host vs shotgun 12/20. (B) Usable BcgI marker yield: 2bRAD flat and
high; shotgun 96k→53k→33k as host rises.

**Figure 7. Real saliva: strain profiles discriminate individuals and are temporally stable.**
(A–C)  — PCoA (species vs strain) coloured by subject,
and per-species strain-level subject R² (8 subjects × 4 timepoints; strain R² 0.833 > species 0.822,
p 2e-4; *Rothia mucilaginosa* R² 0.921; 17/18 species significant). (D–E)
— within- vs between-subject strain distance (0.19 vs 0.44, p 5e-25) and leave-one-timepoint-out host-ID
accuracy (100 % strain vs 94 % species).

**Figure 8. Native 2bRAD confirms all shotgun strains and recovers the low-abundance ones shotgun misses.**
 — paired shotgun↔2bRAD on the same saliva samples. (A) Per sample,
strains detected: shotgun-confirmed (100 % of shotgun calls, 65/65, are in native 2bRAD) plus 2bRAD-only
(128–163/sample). (B) Community relative abundance of shared vs 2bRAD-only strains (2bRAD-only
significantly lower, p 1.2e-23).

**Figure 9. Fast, light, and scalable to whole communities.** (A)  —
per-sample time/memory vs StrainScan (~8×/~11×). (B)  — thread scaling of
build and profile (8.5×/6.5× to 16 threads). (C)  — 55-species
community: per-sample cost flat in #species, ~121–146× faster than a per-species tool run once per species.

**Figure 10. Matches or exceeds StrainScan on its own databases.**
— head-to-head on StrainScan's reference sets: *A. muciniphila*, *P. copri* (precision 1.0 both; recall
0.93 vs 0.24, 0.94 vs 0.90; ~17–23× faster, ~15–24× lighter) and near-clonal *M. tuberculosis*, where
StrainScan did not complete (>3.3 h, >25 GB) and Strain2bScan finished in 0.89 s.

**Figure 11. Systematic head-to-head on the 15-species simulated benchmark.** Both tools build databases from the same
genome pool and profile the same simulated reads, each scored in its own cluster space. **(A–C)** median
single-species precision/recall/F1 vs per-strain sequencing depth (14 species, 204 depth-matched paired
samples): both hold precision 1.0, Strain2bScan reaches full recall by 3× vs StrainScan's 10×.
**(D)** per-species database build time (log): Strain2bScan 0.7–5.1 s vs StrainScan 5–43 min (249–614×).
**(E)** per-sample profile time vs depth in the same emulated container (4–33× faster).
**(F)** multi-species profiling time per community sample — Strain2bScan's one digest-once pass vs
StrainScan's sum over 14 per-species runs (46–105×).

### Supplementary figures and tables

(The former Fig S1 — per-species rank–rank scatter — is now **Fig 2 panel B**.)

**Figure S1. ATCC MSA-1002 DNA-input titration (native 2bRAD).**
— precision 1.0 and full recall to 0.1 ng input.

**Figure S2. Layer-1 gate calibration on the 55-species panel.**  —
species precision/recall vs the marker floor at normal and low depth.

**Table S3. Exploratory clinical oral cohort profiling.**  — 4 clinical
oral 2bRAD samples, 15–17 species and 115–158 strain calls each, ~2 s/sample (no case/control labels
public). Doc: `docs/clinical_oral_exploratory.md`.

**Table S4. Genome quality-control for the 16S/2bRAD motivation panel.**  —
assembly level, CheckM completeness/contamination, contig count and length per genome; high-quality flag.

## Tables

**Table 1. Strain2bScan vs StrainScan on the 15-species simulated benchmark, per species.** Single-species accuracy is the median over depth-matched paired samples (2/3/5-strain mixtures across the 0.5–10× ladder), each tool scored in its own cluster space. Database build cost is per species (Strain2bScan native, arm64; StrainScan `linux/amd64` under emulation). n = genomes in the pool.

| Species | n | S2B P | S2B R | S2B F1 | SS P | SS R | SS F1 | S2B build | SS build |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| *A. muciniphila* | 50 | 1.000 | 0.800 | 0.889 | 0.833 | 0.667 | 0.667 | 2.8 s / 0.16 GB | 17 min / 15 GB |
| *C. difficile* | 47 | 1.000 | 0.800 | 0.889 | 0.800 | 0.667 | 0.667 | 3.2 s / 0.28 GB | 30 min / 14 GB |
| *C. acnes* | 43 | 1.000 | 0.667 | 0.800 | 1.000 | 0.667 | 0.800 | 1.5 s / 0.16 GB | 11 min / 8 GB |
| *E. coli* | 50 | 1.000 | 1.000 | 1.000 | 1.000 | 0.600 | 0.667 | 4.3 s / 0.33 GB | 43 min / 28 GB |
| *F. nucleatum* | 25 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 0.7 s / 0.10 GB | 5 min / 8 GB |
| *L. plantarum* | 50 | 1.000 | 1.000 | 1.000 | 1.000 | 0.600 | 0.750 | 5.0 s / 0.23 GB | 24 min / 17 GB |
| *M. tuberculosis* | 29 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.6 s / 0.25 GB | 16 min / 15 GB |
| *P. dorei* | 15 | 1.000 | 0.667 | 0.800 | 0.833 | 0.583 | 0.667 | 1.0 s / 0.25 GB | 5 min / 13 GB |
| *P. gingivalis* | 43 | 1.000 | 0.667 | 0.800 | 1.000 | 0.667 | 0.800 | 2.3 s / 0.15 GB | 19 min / 10 GB |
| *P. copri* | 19 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.8 s / 0.18 GB | 8 min / 25 GB |
| *S. enterica* | 45 | 1.000 | 1.000 | 1.000 | 0.833 | 0.800 | 0.800 | 4.2 s / 0.36 GB | 34 min / 16 GB |
| *S. aureus* | 48 | 1.000 | 0.667 | 0.800 | 1.000 | 0.500 | 0.667 | 1.7 s / 0.16 GB | 17 min / 12 GB |
| *S. epidermidis* | 50 | 1.000 | 0.667 | 0.800 | 1.000 | 0.600 | 0.750 | 1.6 s / 0.13 GB | 15 min / 8 GB |
| *S. pneumoniae* | 50 | 1.000 | 0.667 | 0.800 | 1.000 | 0.600 | 0.750 | 1.5 s / 0.15 GB | 11 min / 11 GB |

**Median (14 resolvable species):** Strain2bScan P 1.00 / R 0.80 / F1 0.89; StrainScan P 1.00 / R 0.67 / F1 0.75. Build speed-up 249–614×; build memory 43–138× lighter.

**Table 2. Accuracy and per-sample cost vs sequencing depth (single-species, 14 species, 204 paired samples; medians).** Profile time/memory in the same emulated container (Strain2bScan `linux/amd64` vs StrainScan).

| Depth (×) | n | S2B P/R/F1 | SS P/R/F1 | S2B time/mem | SS time/mem |
|--:|--:|:--:|:--:|:--:|:--:|
| 0.5 | 40 | 1.000/0.500/0.667 | 1.000/0.667/0.800 | 0.16 s / 21 MB | 5.3 s / 831 MB |
| 1 | 41 | 1.000/0.667/0.800 | 1.000/0.333/0.500 | 0.24 s / 30 MB | 3.5 s / 831 MB |
| 3 | 41 | 1.000/1.000/1.000 | 1.000/0.600/0.750 | 0.56 s / 60 MB | 4.2 s / 831 MB |
| 5 | 41 | 1.000/1.000/1.000 | 1.000/0.667/0.800 | 0.89 s / 91 MB | 4.8 s / 831 MB |
| 10 | 41 | 1.000/1.000/1.000 | 1.000/1.000/0.889 | 1.66 s / 161 MB | 6.9 s / 832 MB |

**Table 3. Multi-species community profiling (4 samples/depth, matched to the 14 species with StrainScan databases; medians).** Strain2bScan profiles each community in one digest-once pass; StrainScan (no multi-species mode) profiles once per species, so its cost is the sum over 14 databases.

| Community depth | S2B P/R/F1 | SS P/R/F1 | S2B time/mem | SS time/mem |
|---|:--:|:--:|:--:|:--:|
| low | 0.926/0.678/0.782 | 0.898/0.767/0.827 | 1.0 s / 311 MB | 100 s / 1112 MB |
| med | 0.863/0.853/0.869 | 0.911/0.856/0.895 | 4.3 s / 670 MB | 228 s / 1696 MB |
| high | 0.773/0.872/0.819 | 0.814/0.972/0.895 | 8.7 s / 1119 MB | 398 s / 2028 MB |

## References

*Working bibliography compiled from the tools and methods cited in the text. Author-facing note:
verify exact volume/issue/page numbers and publication years against the primary sources before
submission; DOIs/accession numbers should be added.*

1. Liao H, Ji Y, Sun Y. **StrainScan: a highly accurate and sensitive computational tool for strain-level
   detection of the microbiome.** *Microbiome* 2023;11:186.
2. Sun Z, Huang S, Zhu P, *et al.* **Species-resolved sequencing of low-biomass or degraded microbiomes
   using 2bRAD-M.** *Genome Biology* 2022;23:36.
3. Wang S, Meyer E, McKay JK, Matz MV. **2b-RAD: a simple and flexible method for genome-wide
   genotyping.** *Nature Methods* 2012;9(8):808–810.
4. Truong DT, Tett A, Pasolli E, Huttenhower C, Segata N. **Microbial strain-level population structure
   and genetic diversity from metagenomes.** *Genome Research* 2017;27(4):626–638.
5. van Dijk LR, Walker BJ, Straub TJ, *et al.* **StrainGE: a toolkit to track and characterize
   low-abundance strains in complex microbial communities.** *Genome Biology* 2022;23:74.
6. Shaw J, Yu YW. **Rapid species-level metagenome profiling and containment estimation with sylph.**
   *Nature Biotechnology* 2024 (advance online).
7. Ondov BD, Treangen TJ, Melsted P, *et al.* **Mash: fast genome and metagenome distance estimation using
   MinHash.** *Genome Biology* 2016;17:132.
8. Ondov BD, Starrett GJ, Sappington A, *et al.* **Mash Screen: high-throughput sequence containment
   estimation for genome discovery.** *Genome Biology* 2019;20:232.
9. Brown CT, Irber L. **sourmash: a library for MinHash sketching of DNA.** *Journal of Open Source
   Software* 2016;1(5):27.
10. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. **CheckM: assessing the quality of
    microbial genomes recovered from isolates, single cells, and metagenomes.** *Genome Research*
    2015;25(7):1043–1055.
11. Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. **CheckM2: a rapid, scalable and accurate tool for
    assessing microbial genome quality using machine learning.** *Nature Methods* 2023;20(8):1203–1212.
12. Orakov A, Fullam A, Coelho LP, *et al.* **GUNC: detection of chimerism and contamination in
    prokaryotic genomes.** *Genome Biology* 2021;22:178.
13. Bowers RM, Kyrpides NC, Stepanauskas R, *et al.* **Minimum information about a single amplified genome
    (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea.** *Nature Biotechnology*
    2017;35(8):725–731.
14. Huang W, Li L, Myers JR, Marth GT. **ART: a next-generation sequencing read simulator.**
    *Bioinformatics* 2012;28(4):593–594.
15. Marçais G, Kingsford C. **A fast, lock-free approach for efficient parallel counting of occurrences of
    k-mers.** *Bioinformatics* 2011;27(6):764–770.
16. Baker DN, Langmead B. **Dashing: fast and accurate genomic distances with HyperLogLog.** *Genome
    Biology* 2019;20:265.
17. Minkin I, Medvedev P. **Scalable multiple whole-genome alignment and locally collinear block
    construction with SibeliaZ.** *Nature Communications* 2020;11:6327.
18. Eddy SR. **Accelerated profile HMM searches.** *PLoS Computational Biology* 2011;7(10):e1002195.
19. Seemann T. **Barrnap: basic rapid ribosomal RNA predictor.** Software, https://github.com/tseemann/barrnap.
20. Anderson MJ. **A new method for non-parametric multivariate analysis of variance.** *Austral Ecology*
    2001;26(1):32–46. (PERMANOVA)
