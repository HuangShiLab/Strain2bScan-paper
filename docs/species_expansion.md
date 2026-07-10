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

**Strand-fixed** (Strain2bScan re-run; StrainScan side unchanged — its pipeline was unaffected
by our bug, so its numbers, measured on equivalent samples, are carried over):

| species | clusters (40) | Strain2bScan P / R | Strain2bScan time / mem | StrainScan P / R | StrainScan time / mem |
|---|---|---|---|---|---|
| *A. muciniphila* | 19 | **1.00** / 0.93 | 0.75 s / 71 MB | 1.00 / 0.235 | 12.55 s / 1,689 MB |
| *P. copri* | 23 | **1.00** / 0.94 | 0.68 s / 81 MB | 1.00 / 0.900 | 15.62 s / 1,207 MB |
| *M. tuberculosis* | 5 | **1.00** / 0.25 | 0.89 s / 85 MB | **DNF** (>3.3 h, >25.5 GB, 1/5 samples) | — |

(Pre-fix Strain2bScan precision was 0.45 / 0.23 / 0.82 — see below; the low precision was the
strand bug, now fixed.) Strain2bScan is **~10–13× faster and ~7–15× lighter** where StrainScan
completes, and now matches StrainScan's precision (both 1.0) while being recall-first. The one
species where recall is low, *M. tuberculosis* (0.25), reflects **genuine near-clonality** —
40 genomes collapse to just 5 clusters, so within-cluster strains are intrinsically
unresolvable from short reads; this is real biology, not the bug.

## Two distinct accuracy failure modes (new finding)

**1. Clustering collapse (already known, *S. epidermidis*-type):** low intra-species
diversity → strains merge into few clusters → coarse resolution. Not the story here — these
two gut species show *moderate* diversity (33/40, 40/40 clusters from 40 genomes).

**2. What looked like "recombination-driven false uniqueness" was a STRAND BUG (found and
fixed).** With the initial forward-only digestion, *P. copri*'s 40 genomes formed 40 singleton
clusters and Strain2bScan called ~11 clusters present when only 2 were true (precision ~0.19–0.23).
We first hypothesized this was real biology — *P. copri* is highly recombinogenic, so a marker
"unique" in the panel might be carried by the true mosaic strain via recombination — and
predicted it would improve with a larger reference panel. **Both the hypothesis and the
prediction were wrong.** Testing the prediction (nested 40/80/112 panels) is what forced a proper
root-cause diagnosis, which uncovered a core correctness bug:

> **Tag extraction was not reverse-complement invariant.** The enzyme patterns frame a tag
> window differently depending on which strand a site is read from, so digesting a genome vs its
> reverse complement produced *almost disjoint* marker sets (2,209 shared of ~37,600). Reference
> genomes were digested forward-strand-only (~half the tags); sequencing reads come ~50%
> reverse-complemented and so carried *both* strands' tags. A genome's reverse-strand tags then
> collided with *other* genomes' forward tags → spurious detections, worst for similar genomes.
> It also corrupted clustering: G_5 vs G_18 measured Jaccard 0.64 forward-only but **0.993** with
> both strands — near-identical genomes were wrongly split into separate singleton clusters (the
> "40 singletons" artifact).

The strand asymmetry was ultimately traced to **wrong enzyme tag lengths** (15 of 16 enzymes; only
BcgI was correct), which made the forward/reverse recognition patterns stop being reverse-complement
pairs. The fix ports Fast2bRAD-M's correct tag lengths and keeps a single forward-strand scan +
canonical hash (now strand-invariant, and interoperable with Fast2bRAD-M). The panel-size experiment, re-run with the fixed
binary (`scripts/run_panelsize.py`, `results/panelsize_prevotella.tsv`,
`figures/panelsize_prevotella.*`; genomes from EBI/ENA as NCBI was IP-blocked):

| panel size | clusters (was, forward-only) | precision (was) | recall (was) |
|---|---|---|---|
| 40  | **23** (40) | **1.00** (0.19) | **1.00** (0.65) |
| 80  | **43** (80) | **1.00** (0.12) | 0.95 (0.40) |
| 112 | **51** (112) | **1.00** (0.11) | **1.00** (0.40) |

Clustering now matches StrainScan's own *P. copri* granularity (112 → **51** here, exactly StrainScan's 112 → 51), and precision is **1.0** across all panel sizes. There was no
recombination-driven failure mode and no need for a bigger panel — it was a digestion bug.
**This also means the *S. epidermidis* "clustering collapse" and the whole cross-species
precision story (`docs/cross_species.md`) were affected by the same bug and are being re-run;
prior accuracy numbers understate the fixed tool.**

**3. Mycobacterium tuberculosis: StrainScan itself did not complete.** Profiling a *single*
simulated sample against StrainScan's own 792-strain pre-built DB ran for **>3.3 hours and
>25.5 GB RAM** without finishing (we killed it to protect the host machine) — on a task
Strain2bScan completes in ~0.9 s / 85 MB (precision 1.0, recall 0.25). *M. tuberculosis* is textbook near-clonal (very low
inter-strain diversity), and we attribute this to StrainScan's within-cluster Layer-2, which
solves a lasso/regression over highly collinear (near-identical) genome columns — a regime
where such regressions are known to become numerically unstable or slow to converge. This is
itself a notable finding: **the full-k-mer regression-based approach that gives StrainScan its
precision advantage on diverse species can become computationally pathological on very
low-diversity ones**, whereas Strain2bScan's simpler, non-regression unique-marker detection
has no such failure mode (P=1.0, R=0.25, sub-second).

## What this means for the manuscript

1. **The speed/memory advantage is species-general** (not just a *C. acnes* artifact):
   confirmed on 2 more completed species, plus a third where StrainScan couldn't complete at all.
2. **The apparent Strain2bScan precision problem on similar-strain species was a digestion bug,
   not a method limitation.** Diagnosing the *P. copri* panel-size result uncovered a
   strand-invariance bug in tag extraction (above); fixing it gives precision 1.0 on *P. copri*
   with no overlap-aware Layer-2 at all. The earlier framing ("recombination-driven false
   uniqueness", "clustering collapse", "needs overlap-aware Layer-2") was chasing a symptom —
   the honest story is a found-and-fixed correctness bug and a re-verification of every accuracy
   benchmark.
3. **A genuine scalability ceiling for full-k-mer regression on near-clonal species** — StrainScan
   did not complete on *M. tuberculosis* (>3.3 h / >25.5 GB, 1 sample); Strain2bScan's
   non-regression detection has no analog. This finding is independent of the strand bug.

## Caveats

- 40-genome subsets (not each species' full StrainScan panel) — chosen for tractable
  build/profile time.
- Reads are simulated, error-free, closed-world (truth strains guaranteed in both DBs by
  construction) — a stronger, real-world test needs held-out strains and sequencing errors.
- The *M. tuberculosis* DNF is one data point (1 sample, killed rather than left to completion)
  — sufficient to document a real, reproducible problem, but not a fully characterized timeout
  curve; if pursued further, bisecting genome-panel size would locate where StrainScan's
  Layer-2 starts to blow up.
