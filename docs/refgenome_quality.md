# Reference-genome quality vs strain profiling

**Question.** Public reference assemblies vary in completeness, fragmentation and
contamination. Does that variation degrade strain-level profiling? It should: Jaccard-based
clustering (StrainScan's Dashing step, and Strain2bScan's tag-set clustering) is biased by
completeness — an incomplete genome's marker set is ~a subset of its complete twin's, so their
Jaccard (|A∩B|/|A∪B|) falls well below 1 and the same strain fails to cluster (a **spurious
split**), while contamination injects foreign markers.

## Design (15-species reference-degradation gradient)

We vary the **reference DB** quality and keep the **samples fixed**. Degrading the genomes used to
*simulate* the samples too would confound the result: a region missing from both the DB and the
sample cannot mismatch, so the missing region partially cancels. So:

- **Degrade**: the truth strains in a 15-species reference pool across a MIMAG-like gradient
  (completeness 100→95→90→80→70→50%, contamination 0→1→2→5→8→10%, contigs 1→20→50→100→200→400
  co-varying).
- **Keep fixed**: background genomes and the simulated mock samples (reads from the original
  complete genomes at the original abundances — the real metagenome has complete organisms).
- **Measure**: rebuild the cluster DB at each quality level, profile the fixed samples, score
  precision / recall / Bray–Curtis vs reference completeness. Run both **default Jaccard**
  clustering and the optional **`--containment`** mode.

Degradation is simulated by `scripts/degrade.py` (drop random 5-kb windows to the target
completeness, split into the target contig count, append a contaminant fraction taken from an
*E. coli* genome). Because we control it, **completeness/contamination are ground truth**
(no CheckM needed for the x-axis).

## Results (15 species; `results/refqual_15species.tsv`, `results/refqual_15species_containment.tsv`)

### Default Jaccard clustering

Across the 14 resolvable species, median precision and recall decline as references degrade:

| completeness | median precision | median recall | mechanism |
|---|---|---|---|
| 100% | 1.00 | 0.96 | complete references cluster correctly |
| 95% | 0.98 | 0.90 | slight fragmentation begins |
| 90% | 0.84 | 0.80 | incomplete genomes split from complete relatives |
| 80% | 0.71 | 0.74 | many spurious singleton clusters |
| 70% | 0.71 | 0.74 | continued degradation |
| 50% | <0.60 | <0.60 | catastrophic failure for most species |

The near-clonal species *Mycobacterium tuberculosis* is the most extreme case: its single 0.95
cluster shatters into ~18 singletons the moment references degrade, so recall collapses from 1.0
to 0.04–0.12 even at 90–95% completeness. This is a **cluster-fragmentation artifact**, not a
true completeness effect.

### `--containment` mode

The optional `--containment` mode links on **max-containment** (|A∩B| / min(|A|,|B|)) instead of
Jaccard. Because an incomplete genome's marker set is approximately a subset of its complete
relative's, max-containment stays ≈1 and the genomes remain clustered. Results:

| completeness | median precision | median recall |
|---|---|---|
| 100% | 1.00 | 0.96 |
| 95% | 0.98 | 0.95 |
| 90% | 0.92 | 0.92 |
| 80% | 0.86 | 0.86 |
| 70% | 0.75 | 0.75 |

Containment removes the *M. tuberculosis* artifact: the single cluster stays intact and recall
remains 1.0 down to 90% completeness. It converges with Jaccard only at ≤70% completeness, where
references are genuinely low-quality (MIMAG-low; heavy contamination + hundreds of contigs).

Because max-containment merges more aggressively than Jaccard, it is offered as an **opt-in mode**
for reference panels of uneven completeness; the default remains Jaccard, complemented by the
built-in assembly-quality filter (`--min-tag-fraction` / `--max-contigs`).

## Historical single-species *C. acnes* table

An earlier version of this analysis degraded only the 14 *C. acnes* truth strains in the 64-genome
panel. The qualitative pattern was identical — precision robust to moderate degradation, recall
and abundance error degrading progressively, catastrophic failure at 50% completeness / 10%
contamination / 400 contigs — and is archived in the git history. The 15-species gradient
reported here is the current manuscript version.

## Caveats
- Reads are error-free; an error model would add noise on top of the quality effect.
- Completeness/contamination are designed values; report CheckM2 estimates of the simulated
  assemblies for the final figure (Linux).
- The two modes score at the same cluster resolution; near-clonal species are compared at the
  honest cluster level.

## Running it (one `make` target per arm)

**Strain2bScan arm + figure** (any platform; needs the built binary):
```bash
export STRAIN2BSCAN_BIN=/path/to/strain2bscan      # or have `strain2bscan` on PATH
make refqual        # downloads the 15-species pool + runs the gradient
# -> results/refqual_15species.tsv, figures/refqual_figure.{png,pdf}
```

**Containment arm:**
```bash
make refqual CONTAINMENT=1
# -> results/refqual_15species_containment.tsv
```

**CheckM2 validation** (Linux HPC; for the published x-axis = measured scores):
```bash
conda create -n checkm2 -c bioconda -c conda-forge checkm2 && conda activate checkm2
checkm2 database --download        # ~3 GB diamond DB, once
make refqual-checkm2               # -> results/refqual_checkm2.tsv (designed vs measured)
```

Then copy `work/refqual/*.tsv` and `work/figures/*` back into `results/` and `figures/` to
refresh the committed snapshots.
