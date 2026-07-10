# Strain2bScan — Results (draft)

> Numbers below are from the corrected-enzyme binary (Fast2bRAD-M tag lengths, single-strand scan
> + canonical hash). Figures are `figures/*`; source tables are `results/*.tsv`; the experimental
> design and per-experiment conclusions are in `docs/benchmark_datasets.md`. Section order follows
> the figure plan in `docs/results_figure_plan.md`. All benchmark reads except the saliva chapter
> are simulated (150 bp, error-free, closed-world); this caveat is stated in the Discussion.

## Overview of Strain2bScan (Fig 1)

Strain2bScan resolves within-species strains from the sparse marker set released by type-IIB
restriction (2bRAD) digestion rather than from a full k-mer profile. For each species, reference
genomes are digested in silico with up to 16 type-IIB enzymes into single-copy 32–38 bp tags —
roughly 1–2 % of the genome and 50–100× sparser than a full k-mer index — clustered into
within-species groups at 0.95 Jaccard similarity (MinHash-accelerated), and reduced to a
cluster × marker database of species-core, cluster-specific and strain-specific tags (Fig 1A).
A sample is profiled by digesting its reads once, counting canonical markers, gating on
species-specific markers (Layer-1), and detecting and quantifying strains within each present
species from unique markers (Layer-2), returning a strain-level profile (Fig 1B). The pipeline
accepts two data modes — in-silico digestion of conventional shotgun metagenomes, or native BcgI
2bRAD experimental libraries whose reads *are* the tags. Because the tag lengths and recognition
patterns are identical to Fast2bRAD-M / `2bRADExtraction.pl`, the Fast2bRAD-M species layer and
the Strain2bScan strain layer operate in one shared tag space (Fig 1, band).

## Accurate strain profiling across species (Fig 2)

On real reference panels with simulated strain mixtures, Strain2bScan achieved **precision 1.0**
for every species tested — *Cutibacterium acnes* (64 genomes → 16 clusters), *Staphylococcus
aureus* (60 → 17) and *S. epidermidis* (60 → 10) — with recall 0.75, 0.79 and 1.0 and abundance
error (Bray–Curtis) 0.24, 0.33 and 0.02, respectively (Fig 2). Precision remained 1.0 even for
low-diversity *S. epidermidis*, where 60 genomes collapse into only 10 clusters: the profiler did
not over-detect closely related strains. Recall, not precision, is the species-dependent axis and
tracks *genuine* intra-species diversity — high where strains are distinct, lower where the
reference panel is dominated by near-identical genomes.

## Detection down to 0.5× coverage, matching StrainScan (Fig 3)

Sensitivity was benchmarked by profiling a single *C. acnes* strain present in both tools'
databases across a coverage ladder. Strain2bScan detected the strain at **0.5× per-strain
coverage** and missed it only at 0.1× — the **same detection onset as StrainScan** (Fig 3). The
sparse 2bRAD marker set therefore carries no low-depth penalty at realistic coverages; an earlier
report of reduced low-depth sensitivity was an artifact of the forward-only digestion bug, now
fixed.

## The 2bRAD enzyme set is a resolution/cost knob (Fig 4)

Because each added enzyme yields more tags, the enzyme set is a tunable knob. On *C. acnes*,
**BcgI alone** (one enzyme) already gave precision 1.0 with recall 0.56; recall rose to 0.75 by
**four enzymes** and then plateaued, while strain-specific marker count and build time kept rising
roughly linearly (882 → 1 524 → 7 400 markers and 0.6 → 1.8 → 6.4 s from 1 → 4 → 14 enzymes;
Fig 4). Cluster count was unchanged by enzyme number (16 throughout): more enzymes add markers and
recall, not resolution. A **~4-enzyme set is the sweet spot** — most of the recall at a fraction
of the markers and compute of the full 16 — and, critically, single-enzyme BcgI operation is what
lets native BcgI 2bRAD-M libraries be profiled directly.

## Robust to reference panel size and moderate reference degradation (Fig 5)

Strain2bScan was insensitive to reference database size: on nested *Prevotella copri* panels of
40, 80 and 112 genomes it held **precision 1.0** at every size, with recall 1.0 / 0.95 / 1.0 and
cluster counts 23 / 43 / 51 (Fig 5A). The 112-genome partition (51 clusters) matches the
clustering of StrainScan's own pre-built *P. copri* database, an external check that the
MinHash-accelerated clustering is correct. Robustness to reference *quality* has a floor
(Fig 5B): degrading the truth strains' reference genomes held precision at 1.0 down to **70 %
completeness** (recall declining from 0.75 to 0.19 as completeness fell), but detection collapsed
entirely at 50 % completeness / 10 % contamination / 400 contigs (precision 0). This motivates the
built-in assembly-quality filter (contig count and single-copy-tag count relative to the
conspecific median) that flags low-quality references before clustering.

## Fast, light, and scalable to whole communities (Fig 6)

Per sample, Strain2bScan profiled the *C. acnes* panel in **0.86 s using 78 MB** versus **7.06 s
and 828 MB** for StrainScan — **~8× faster and ~11× lighter** (Fig 6A). The memory gap is
structural: StrainScan counts the full k-mer set, whereas Strain2bScan streams reads and keeps
only the sparse 2bRAD markers. Both database build and profiling parallelise cleanly with
dependency-free threads (8.5× and 6.5× from 1 → 16 threads; Fig 6B). The decisive advantage is at
community scale: because a sample is digested **once** and each additional species is a cheap
hash-set match rather than a re-count, per-sample cost is essentially **independent of the number
of species**. On a **55-species** community Strain2bScan profiled each sample in ~2.6–3.1 s and
ran **121–146× faster** than running a per-species tool once per species (~132× at 100 samples,
~146× at ≥200 samples; Fig 6C) — the regime where full-k-mer tools, which have no multi-species
mode, must re-count every sample once per species.

## Matches or exceeds StrainScan on its own databases (Fig 7)

Head-to-head against StrainScan on StrainScan's own reference sets, both tools reached
**precision 1.0** on *Akkermansia muciniphila* (40 genomes → 19 clusters) and *P. copri* (40 →
23), but Strain2bScan **matched or exceeded StrainScan's recall** (0.93 vs 0.24; 0.94 vs 0.90)
while running **~17–23× faster** and **~15–24× lighter** (0.7–0.9 s / 71–85 MB vs 12.6–15.6 s /
1.2–1.7 GB; Fig 7). On near-clonal *Mycobacterium tuberculosis* (40 genomes → only 5 clusters),
**StrainScan did not complete** (>3.3 h, >25 GB) whereas Strain2bScan finished in 0.89 s; its low
recall there (0.25) reflects the intrinsic near-clonality of the species — a resolution limit, not
a tool failure — and is reported at precision 1.0.

## Real-data validation on paired saliva metagenomes (Fig 8 — planned)

The benchmarks above are simulated and closed-world. The planned real-data chapter uses a paired
**shotgun + BcgI-2bRAD saliva** dataset (8 subjects × 4 within-day timepoints) to test three
things simulation cannot: (i) **individual discrimination** — whether strain-level community
structure separates subjects better than species-level structure (PERMANOVA R²), reproducing a
prior lab result at a fraction of the compute; (ii) **paired shotgun↔2bRAD concordance**, which
validates both real-read performance and the in-silico-digestion ↔ native-2bRAD equivalence
without external ground truth; and (iii) **within-subject temporal stability**. Design and
validation logic are detailed in `docs/real_metagenome_plan.md`.
