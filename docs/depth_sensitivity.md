# Depth-sensitivity: Strain2bScan vs StrainScan (the <1× question)

**Question (user):** at <1× depth for a strain, is Strain2bScan *more* sensitive than StrainScan?

**Answer: No — StrainScan is more sensitive at low depth**, for digitally-digested shotgun
metagenomes. Confirmed empirically below and expected from first principles: markers observed
≈ depth × (#marker loci); 2b tags are ~50–100× fewer loci than the full k-mer set, so at low
depth Strain2bScan runs out of observed markers first.

## Experiment
Single C. acnes strain **GCF_009737125** (present in *both* tools' DBs), simulated shotgun
reads (150 bp, error-free) at 0.1× / 0.5× / 1× / 5× genome coverage. Strain2bScan profiled
against the 64-genome panel (detection floor 10 and a permissive 3); StrainScan against its
275-strain DB with its low-depth modes (`-l 2` <1×, `-l 1` at 1×).

| depth | Strain2bScan (floor 10) | Strain2bScan (floor 3) | StrainScan |
|---|---|---|---|
| 0.1× | miss | miss | miss |
| **0.5×** | miss | miss | **DETECT** |
| 1× | miss | DETECT | DETECT |
| 5× | DETECT | DETECT | DETECT |

**Detection-onset depth: StrainScan ≈ 0.5×, Strain2bScan ≈ 1× (permissive) to ≈ 5× (default).**
Even with a permissive floor of 3, Strain2bScan cannot detect at 0.5× — fewer than 3 of the
strain's unique tags are observed. This is a fundamental marker-sparsity limit, not a
threshold artifact.

## Implication for the manuscript (important reframe)
- Do **not** claim low-depth sensitivity superiority for digital 2bRAD — the data refute it.
- The correct, defensible claim: *Strain2bScan matches StrainScan at **sufficient depth
  (≈≥5×)** while being ~2× faster and ~8× lighter; below ~1× its sparse markers lose
  sensitivity faster than full k-mers.* The depth sweep **defines "sufficient depth."**
- The genuine low-input advantage of 2bRAD belongs to **real wet-lab BcgI 2bRAD libraries**,
  where sequencing is concentrated on tags (deep per-tag coverage at low total cost) — a
  separate experiment (claim 3 / real 2bRAD data), not digital digestion of shotgun reads.

## To strengthen for publication
- Replicate across strains/species (S. aureus, S. epidermidis) and add ≥3 replicates per
  depth; report precision/recall + abundance error (not just detect/miss).
- Add a realistic error model (errors hurt 2b tags slightly more — longer markers).
- For a fully fair head-to-head, build both DBs on the identical panel (StrainScan build
  needs Linux/dashing).
- Add the matched **real-2bRAD-vs-shotgun at equal sequencing budget** experiment — that is
  where 2bRAD can actually win, and it is the stronger story.
