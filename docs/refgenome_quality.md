# Reference-genome quality vs strain profiling (supplementary)

**Question.** Public reference assemblies vary in completeness, fragmentation and
contamination. Does that variation degrade strain-level profiling? It should: Jaccard-based
clustering (StrainScan's Dashing step, and Strain2bScan's tag-set clustering) is biased by
completeness — an incomplete genome's marker set is ~a subset of its complete twin's, so their
Jaccard (|A∩B|/|A∪B|) falls well below 1 and the same strain fails to cluster (a **spurious
split**), while contamination injects foreign markers.

## Design (clean isolation of the reference-quality effect)

We vary the **reference DB** quality and keep the **samples fixed**. Degrading the genomes
used to *simulate* the samples too would confound the result: a region missing from both the
DB and the sample cannot mismatch, so the missing region partially cancels. So:

- **Degrade**: the 14 *C. acnes* truth strains in the 64-genome panel, across a MIMAG-like
  gradient (completeness 100→50%, contamination 0→10%, contigs 1→400 co-varying).
- **Keep fixed**: the 50 background genomes, and the 5 mock samples (reads from the original
  complete genomes at the original abundances — the real metagenome has complete organisms).
- **Measure**: rebuild the cluster DB at each quality level, profile the 5 fixed samples,
  score precision / recall / Bray–Curtis vs reference completeness.

Degradation is simulated by `scripts/degrade.py` (drop random 5-kb windows to the target
completeness, split into the target contig count, append a contaminant fraction taken from an
*E. coli* genome). Because we control it, **completeness/contamination are ground truth**
(no CheckM needed for the x-axis); `CheckM2` validation of the simulated assemblies is a
Linux follow-up (it has no osx-arm64 build, like the StrainScan build deps).

The simulator is faithful: an 80%-complete *C. acnes* genome yields **28,397 single-copy tags
vs 34,387 for the original (0.83×)** — tags drop ~proportionally to completeness, which also
validates Strain2bScan's tag-count completeness proxy (`quality.rs`).

## Results (Strain2bScan arm; `results/refqual_degradation.tsv`, `figures/refqual_figure.*`)

| completeness | contamination | contigs | clusters | precision | recall | Bray–Curtis |
|---|---|---|---|---|---|---|
| 100% | 0% | 1 | 60 | **0.90** | 0.56 | **0.25** |
| 95% | 1% | 20 | 61 | 0.80 | 0.50 | 0.36 |
| 90% | 2% | 50 | 61 | 0.80 | 0.50 | 0.39 |
| 80% | 5% | 100 | 61 | 0.80 | 0.50 | 0.37 |
| 70% | 8% | 200 | 61 | 0.80 | 0.50 | 0.36 |
| 50% | 10% | 400 | 61 | **0.64** | 0.44 | **0.46** |

- **Precision and abundance error are the sensitive readouts.** Even mild degradation
  (95% / 1%) drops precision 0.90→0.80 and nearly doubles abundance error (Bray–Curtis
  0.25→0.36) — contamination adds spurious markers that trigger false-positive clusters.
- **A clear threshold appears at low quality**: at 50% completeness / 10% contamination,
  precision falls to 0.64, recall to 0.44, Bray–Curtis to 0.46.
- Cluster count ticks 60→61 (a degraded truth strain splitting off — the spurious-split
  mechanism), modest here because the panel's clusters are mostly singletons already.

This supports the hypothesis: **below a quality threshold (~MIMAG medium quality), reference
assembly quality has a non-negligible effect on strain profiling.** Recall is dominated by the
sparse-marker low-abundance-minor limit (≈0.5 even at 100%), so precision and Bray–Curtis are
the cleaner indicators of the reference-quality effect.

## StrainScan arm (run on Linux)

`scripts/run_refqual_strainscan.py` builds a StrainScan DB at each quality level on the same
degraded panels and profiles the same samples. It is **Linux-only** (StrainScan_build needs
`dashing`, no osx-arm64 build). `scripts/plot_refqual.py` overlays both tools once that TSV
exists. StrainScan uses the same Jaccard/Dashing clustering, so the same completeness bias is
expected; the figure quantifies whether full-k-mer density buffers it relative to sparse tags.

## Caveats
- Reads are error-free; an error model would add noise on top of the quality effect.
- Completeness/contamination are designed values; report CheckM2 estimates of the simulated
  assemblies for the final figure (Linux).
- The two arms score at slightly different resolutions (Strain2bScan at cluster level,
  StrainScan at strain level) — noted on the figure.

## Running it (one `make` target per arm)

**Strain2bScan arm + figure** (any platform; needs the built binary):
```bash
export STRAIN2BSCAN_BIN=/path/to/strain2bscan      # or have `strain2bscan` on PATH
make refqual        # downloads the C. acnes 64-panel + 5 mock reads + E. coli, runs the gradient
# -> work/refqual/refqual_degradation.tsv, work/figures/refqual_figure.{png,pdf}
```

**StrainScan arm** (Linux HPC; StrainScan_build needs `dashing`):
```bash
conda env create -f env/strainscan-linux.yml && conda activate strainscan
export STRAINSCAN=/path/to/StrainScan  STRAINSCAN_PY=$(which python)
make refqual              # ensure data + Strain2bScan arm are in place
make refqual-strainscan   # StrainScan DB per quality level; plot_refqual.py overlays both tools
```

**CheckM2 validation** (Linux HPC; for the published x-axis = measured scores):
```bash
conda create -n checkm2 -c bioconda -c conda-forge checkm2 && conda activate checkm2
checkm2 database --download        # ~3 GB diamond DB, once
make refqual-checkm2               # -> work/refqual/refqual_checkm2.tsv (designed vs measured)
```

Then copy `work/refqual/*.tsv` and `work/figures/*` back into `results/` and `figures/` to
refresh the committed snapshots.
