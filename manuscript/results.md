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
with its species' 16S CI. The per-pair rank–rank scatter (Fig 2B) shows the mechanism: 2bRAD distances
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

## Part I — Native 2bRAD-M for low-biomass, high-host microbiomes

### A tunable enzyme set enables native BcgI 2bRAD operation (Fig 5)

Because each added type-IIB enzyme releases more tags, the enzyme set is a resolution/cost knob — and
this holds **across all 14 resolvable species** of the simulation pool, not just one species (Fig 5).
Sweeping the enzyme set from 1 → 14 enzymes, **precision stayed 1.0 at every step in every species**,
while median **recall rose monotonically — 0.50 (BcgI alone) → 0.61 → 0.72 → 0.89 → 1.00** at 1/2/4/8/14
enzymes — as the strain-specific marker yield grew ~13-fold (median 976 → 12 647 markers; Fig 5A,B).
Cluster count was invariant with enzyme number (more enzymes add markers and recall, not resolution),
confirming the single-species observation generalises. The recall/enzyme trade-off is species-dependent
(BcgI alone suffices for diverse species but under-resolves near-clonal ones), so **~4 enzymes is a
practical sweet spot** that recovers most recall at a fraction of the markers, with full recall available
at higher enzyme counts. Decisively for this pillar, **single-enzyme BcgI operation is what lets native
BcgI 2bRAD-M experimental libraries be profiled directly**, since their reads are BcgI tags.

### Strain-level identification and quantification on a real 20-strain DNA mock (Fig 6)

Fig S3 established species-level recovery under host load; here we test whether native 2bRAD **resolves
and quantifies the individual strains** of the ATCC MSA-1002 (even) and MSA-1003 (staggered) mocks. We
built a single **combined 120-genome database** in which each of the 20 mock species carries its true
ATCC genome plus **5 additional conspecific strains** (high-quality, within-species ANI > 95 %), so a
correct call must pick the right genome out of six same-species candidates; the mocks contain only the
ATCC strain, so any other cluster is a false call. Against this database Strain2bScan resolved **19/20
strains** to their correct ATCC cluster at 90–99 % human DNA (100 ng) with near-zero false strain calls
(Fig 6A). The one exception, *L. gasseri*, has fewer than 10 unique BcgI tags in the combined tree —
a single-enzyme resolution limit (Fig 5) — and is still detected at species level. Detection scaled with
DNA input, reaching 19/20 by 0.1 ng (Fig 6B). Species abundance, read from each cluster's marker
**depth** (∝ genome copy number) normalised across species, recovered the mock design: roughly uniform
for the even MSA-1002 mix and spanning three orders of magnitude for the staggered MSA-1003 mix
(Fig 6C). Because the > 95 %-ANI decoys are near-clonal, they fragment the sparse single-enzyme BcgI
tree and strip unique markers; **`--containment` clustering (Fig 4) merges them and restores detection** —
the same mechanism, now on real 2bRAD reads. Crucially, this native-2bRAD accuracy is **on par with
conventional shotgun of the same mock**: in-silico-digested shotgun resolves 20/20 strains and the
shotgun k-mer/mapping tools StrainScan and inStrain likewise resolve 20/20 (Fig 12) — so 2bRAD's ~50–100×
reduced marker set costs essentially nothing in strain-identification accuracy, while the wet-lab
reduction is what preserves that accuracy under the host contamination and low input where shotgun fails
(Fig 6A,B; Fig S3).

### 2bRAD holds strain recovery under host contamination where shotgun collapses (Fig S3)

We tested the low-biomass claim on the ATCC MSA-1002 20-strain even mock with **real DNA and known
truth**, profiled against a 62-species BcgI panel (20 mock + 42 decoys). Across a **host-contamination
series (0/90/99 % human DNA)**, native BcgI 2bRAD and in-silico-digested shotgun (WMS) of the same mock
were profiled with the same database. **Precision was 1.0 for both data types at every contamination
level** — zero false strain calls even at 99 % human DNA, with out-of-panel relatives correctly held in
the detected-not-resolvable tier. But recall diverged sharply: native 2bRAD retained **full 20/20 species
recall at 99 % host**, whereas shotgun fell to **12/20** (Fig S3A). The mechanism is marker yield (Fig S3B):
native 2bRAD delivers ~340–370k usable BcgI markers regardless of host fraction — the reduction happens at
the wet-lab step, before host swamps the library — while digesting an already host-dominated shotgun
library in silico yields only 96k → 53k → 33k markers as host rises. In a separate DNA-input titration,
native 2bRAD held precision 1.0 with full recall down to **0.1 ng** input (Fig S1). This is the
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
low-abundance tail — the same enrichment advantage seen on the controlled mock (Fig S3), now on real
clinical material. (Strain2bScan runs directly on an oral clinical cohort — 15–17 species, ~2 s per
sample — as an application demonstration; Table S3.)

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
species from the mock at 0 % host (Fig S3) and its saliva strain calls are a confirmed subset of the
native-2bRAD calls (Fig 8) — so the fast shotgun mode and the sensitive 2bRAD mode are two faces of one
validated method.

### Systematic head-to-head on a 15-species simulated benchmark (Fig 11, Table 1)

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

### The same strain-level identification via in-silico shotgun (Fig 12)

Digesting conventional shotgun of the **same** MSA-1002/1003 mocks in silico (all 16 enzymes) and matching
against the all-enzyme combined 120-genome tree, Strain2bScan resolved **20/20 strains** to their correct
ATCC cluster with **zero false positives** on both the even and staggered mixes (Fig 12A) — the richer
all-enzyme marker set avoids the sparse-tree fragmentation seen with single-enzyme BcgI. Depth-based
species abundances again recovered the community structure: uniform for MSA-1002 and spanning
*E. coli* 25 % → *Bifidobacterium* 0.02 % (by mass) for the staggered MSA-1003 (Fig 12B). With
`--min-abundance 0` (no low-abundance floor, so no true strain is dropped) precision is governed by
marker **coverage**: every true cluster retains ≥ 28 % breadth while the occasional near-sibling
cross-call sits at ~10 %, so `--min-coverage 0.2` removes all false positives while keeping every strain.
The native-2bRAD (Fig 6) and shotgun (Fig 12) modes thus deliver the same strain-level identification and
abundance on a real DNA mock from one combined database — two entry points to one method.

We benchmarked the shotgun mode head-to-head against the two standard shotgun strain tools,
**StrainScan** (reference k-mer) and **inStrain** (read-mapping + microdiversity), on the same
MSA-1002-even shotgun sample and 120-genome reference (Fig 12C,D). **All three resolved 20/20 strains**
— accuracy is comparable — but the compute cost differs by orders of magnitude: Strain2bScan profiled the
sample in **~66 s** (one digest-once pass over the combined tree) versus StrainScan **~35 min** (20
per-species k-mer runs; Docker `linux/amd64` emulation) and inStrain **~52 min** (bowtie2 mapping of
~19 M read pairs + microdiversity profiling). On the more diverse staggered MSA-1003 sample **inStrain did
not finish** (35 % after 1 h 08, ~5–6 h projected), whereas Strain2bScan completed it in ~66 s. So on
conventional shotgun Strain2bScan matches the accuracy of both established tools at a fraction of the wall
time (`results/mock_tool_comparison.tsv`) — the compute advantage that makes it practical for
host-contaminated and cohort-scale shotgun where per-sample cost dominates.
