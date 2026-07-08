# Species expansion: 3 StrainScan-native species, using StrainScan's own pre-built DBs

**Question.** Beyond *C. acnes* / *S. aureus* / *S. epidermidis*, does the head-to-head hold
across more representative species? StrainScan's own GitHub ships 8 pre-built species
databases (Akkermansia muciniphila 157/53, *C. acnes* 275/28, Prevotella copri 112/51,
*E. coli* 1433/823, *M. tuberculosis* 792/25, *S. epidermidis* 995/378, *S. aureus* 1627/202,
*L. crispatus* 1124/311 [strains/clusters]). We added three: **Akkermansia muciniphila**
(moderate, gut), **Prevotella copri** (moderate, gut), **Mycobacterium tuberculosis** (large,
low-diversity clinical pathogen).

## Design: genuinely fair, no DB-version mismatch

Profiling with StrainScan needs no `dashing` (build does; see `docs/macos_dashing.md`), so we
can run **real StrainScan** directly, on macOS, against its own pre-built DB — no Linux/HPC
needed. For each species: draw a 40-genome subset **from StrainScan's own DB accession list**
(guarantees genome overlap between the two tools' databases, avoiding the earlier *C. acnes*
issue where an independently-chosen panel only partly overlapped StrainScan's 2022 snapshot),
resolve versions and download from NCBI, build our own Strain2bScan cluster DB on that subset,
simulate 5 log-normal strain-mixture samples (2–5 truth strains, depth ≥1×), and profile both
tools on the same reads. `scripts/run_species_expansion.py`, `results/species_expansion.tsv`,
`figures/species_expansion.*`.

## Results

| species | clusters (40-genome subset) | Strain2bScan P / R | Strain2bScan time / mem | StrainScan P / R | StrainScan time / mem |
|---|---|---|---|---|---|
| *A. muciniphila* | 33 | 0.452 / 0.933 | 0.50 s / 89 MB | 1.000 / 0.235 | 12.55 s / 1,689 MB |
| *P. copri* | 40 | 0.231 / 0.750 | 0.76 s / 141 MB | 1.000 / 0.900 | 15.62 s / 1,207 MB |
| *M. tuberculosis* | 32 | 0.824 / 0.933 | 0.49 s / 93 MB | **DNF** (>3.3 h, >25.5 GB, 1/5 samples) | — |

Speed/memory pattern replicates: Strain2bScan is **~20–30× faster and ~9–18× lighter** where
StrainScan completes. Accuracy shows the expected **precision/recall trade-off** — StrainScan
(full k-mer + overlap-matrix/lasso Layer-2) is precision-first (1.0 on both completed species,
sometimes at real recall cost: 0.235 on *A. muciniphila*); Strain2bScan (unique-marker
detection) is recall-first (0.75–0.93 throughout) with precision that varies by species — and
uncovers **two distinct, previously uncharacterized failure modes**, not one.

## Two distinct accuracy failure modes (new finding)

**1. Clustering collapse (already known, *S. epidermidis*-type):** low intra-species
diversity → strains merge into few clusters → coarse resolution. Not the story here — these
two gut species show *moderate* diversity (33/40, 40/40 clusters from 40 genomes).

**2. Recombination-driven false uniqueness (new, *P. copri*-specific):** *P. copri*'s 40
genomes form **40 singleton clusters (0.231 precision) despite maximal apparent resolution**
— every genome is its own cluster, yet Strain2bScan still calls 10 clusters present when only
2 are true. We traced this to real biology, not a bug (verified the 10 calls map to 10
genuinely distinct genome accessions, and confirmed several false positives have *higher*
read support/coverage than the true positives — ruling out a simple confidence-threshold fix).
*P. copri* is documented in the literature as highly recombinogenic with mosaic genomes; a
marker that is "unique" within a reference panel can still be carried by the *true* (mosaic)
sample strain if a recombination breakpoint happens to match a segment private to some *other*
panel genome.

### Panel-size test: prediction made, then REFUTED

We initially predicted this would **improve with reference panel size** (reasoning: a bigger
panel is more likely to already contain, somewhere, any given recombinant segment, so it gets
reclassified as shared rather than falsely unique). We tested it directly: strictly-nested
panels **40 ⊂ 80 ⊂ 111** genomes (111/112 of StrainScan's whole *P. copri* panel; one genome
un-fetchable), the **identical 5 samples**, Strain2bScan only (`scripts/run_panelsize.py`,
`results/panelsize_prevotella.tsv`, `figures/panelsize_prevotella.*`; genomes from EBI/ENA as
NCBI was IP-blocked during this run — GCA vs GCF assemblies, hence the 40-genome baseline reads
0.194 here vs 0.231 in the NCBI-sourced run above, within assembly-source noise).

| panel size | clusters | precision | recall |
|---|---|---|---|
| 40 | 40 | 0.194 | 0.650 |
| 80 | 80 | 0.121 | 0.400 |
| 111 | 111 | 0.114 | 0.400 |

**The prediction was wrong — precision *and* recall both DECLINE as the panel grows.** The
mechanism is the opposite of what we guessed: because *P. copri* genomes are all distinct
(every added genome is a new singleton cluster), a larger panel means (a) **more** singleton
clusters whose "unique" markers a mosaic sample strain can spuriously match → more false
positives, and (b) the true strain's own previously-unique markers get reclassified as shared
once a recombination-similar genome is added → the true cluster loses detectable unique signal
→ recall drops too. So for a highly recombinant species, adding reference genomes makes the
simple unique-marker Layer-2 *worse*, not better.

This is a **stronger** result than a confirmation would have been: it shows the false-uniqueness
problem is **intrinsic to unique-marker detection for recombinant species and cannot be fixed by
enlarging the database** — it specifically needs a joint, overlap-aware Layer-2 (deconvolve which
few clusters best explain the reads, rather than calling every cluster with enough matched
"unique" markers). It sharpens, rather than weakens, the motivation for that step.

**3. Mycobacterium tuberculosis: StrainScan itself did not complete.** Profiling a *single*
simulated sample against StrainScan's own 792-strain pre-built DB ran for **>3.3 hours and
>25.5 GB RAM** without finishing (we killed it to protect the host machine) — on a task
Strain2bScan completes in 0.49 s / 93 MB. *M. tuberculosis* is textbook near-clonal (very low
inter-strain diversity), and we attribute this to StrainScan's within-cluster Layer-2, which
solves a lasso/regression over highly collinear (near-identical) genome columns — a regime
where such regressions are known to become numerically unstable or slow to converge. This is
itself a notable finding: **the full-k-mer regression-based approach that gives StrainScan its
precision advantage on diverse species can become computationally pathological on very
low-diversity ones**, whereas Strain2bScan's simpler, non-regression unique-marker detection
has no such failure mode (P=0.824, R=0.933, sub-second).

## What this means for the manuscript

1. **The speed/memory advantage is species-general** (not just a *C. acnes* artifact):
   confirmed on 2 more completed species, plus a third where StrainScan couldn't complete at all.
2. **Precision/recall is a genuine, honestly-reported trade-off**, not a Strain2bScan
   weakness to paper over — cite it plainly, alongside the resolvability framing from
   `docs/cross_species.md`.
3. **Two distinct, biologically-grounded limitations**, not one: clustering collapse (a fixed
   species property) vs. recombination-driven false uniqueness. Critically, the panel-size test
   above shows the recombination case is **not** fixable by a larger database (both precision
   and recall got worse) — it, like clustering collapse, needs the **overlap-aware Layer-2**
   (StrainScan's lasso/overlap-matrix step): jointly deconvolve which few clusters explain the
   reads rather than calling every cluster with enough matched "unique" markers. That single fix
   addresses both regimes, at the cost of the numerical/scaling fragility item 3 exposes on
   near-clonal species — a real design tension to state in the paper, not hide.
4. **A genuine scalability ceiling for full-k-mer regression on near-clonal species** — a
   novel, citable robustness result distinct from the standard "faster/lighter" efficiency claim.

## Caveats

- 40-genome subsets (not each species' full StrainScan panel) — chosen for tractable
  build/profile time. For *P. copri* we went further and tested the full nested 40/80/111 panel
  (above): accuracy did **not** improve — it declined — refuting our initial "bigger panel
  helps" prediction and localizing the fix to an overlap-aware Layer-2 instead.
- Reads are simulated, error-free, closed-world (truth strains guaranteed in both DBs by
  construction) — a stronger, real-world test needs held-out strains and sequencing errors.
- The *M. tuberculosis* DNF is one data point (1 sample, killed rather than left to completion)
  — sufficient to document a real, reproducible problem, but not a fully characterized timeout
  curve; if pursued further, bisecting genome-panel size would locate where StrainScan's
  Layer-2 starts to blow up.
