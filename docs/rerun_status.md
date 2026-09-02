# What needs re-running, and what does not

Two engine changes on `strainscan-port` landed after the committed results were produced. This
records which figures they touch, measured rather than assumed, so nothing is re-run on
suspicion and nothing stale is quoted.

## 1. `6329f9e` — a digestion locus shared by several enzymes is counted once

Ten of the sixteen type-IIB enzymes emit 27 bp tags and their recognition sites genuinely
coincide (every AloI site with a C at offset 16 is also a BsaXI site). Scanning per enzyme gave
such a locus a copy number of 2 or 3, so the single-copy filter dropped it — a systematic hole
in multi-enzyme panels.

**Measured effect** (19-cluster *E. coli* panel, 14-enzyme set, synthetic 150 bp read mixtures):

| | old | new | truth |
|---|---|---|---|
| database markers | 166,224 | 176,734 | — |
| unique markers | 104,488 (62.9%) | 111,046 (62.8%) | — |
| 70/30 mixture, abundance | 72.13 / 27.87 % | 72.22 / 27.78 % | 70 / 30 |
| 90/10 mixture, abundance | 90.95 / 9.05 % | 90.95 / 9.05 % | 90 / 10 |
| depth | 7.320 / 2.828 x | 7.330 / 2.820 x | — |
| coverage | 98.46 / 93.74 % | 98.47 / 93.74 % | — |
| `support` | 4,982 / 8,375 | 5,327 / 8,907 | — |

So: the database gains ~6.3% of markers that were being wrongly discarded, and `support` rises
with it, but **abundance moves by at most 0.09 percentage points and is unchanged at 90/10**.
Coverage and depth are stable to the third digit. The commit's own summary says the same thing
— "panel recovery, not an abundance correction" — and this is an independent check of it.

**Single-enzyme databases are byte-identical.** Verified: BcgI gives 25,425 markers before and
after.

### Consequence

| result | enzyme set | action |
|---|---|---|
| Fig 9 mock host contamination | BcgI | none — byte-identical |
| MSA-1002 titration | BcgI (native 2bRAD) | none |
| saliva (Fig 7/10) | BcgI | none *from this change* — but see §2 |
| clinical exploratory | BcgI | none — verified `profile_clinical.py:19` |
| simulated multi-species, panel size, species expansion, cross-species, reference quality | 14-enzyme | rebuild DBs; abundance-based numbers will barely move, any `support`/marker count will rise ~6% |

Nothing here changes a recall, precision, Bray–Curtis or abundance figure by a visible amount.
Rebuild the multi-enzyme databases when convenient; do not re-run the abundance figures on
account of this change alone.

## 2. Saliva community matrices were built on `support`

Independent of the engine, and more serious. See `saliva_individual_discrimination.md` for the
full note. In short: `profile_saliva.py` parsed multi-profile's output positionally, stopped at
index 4, and so never captured `sample_fraction`; every downstream matrix used index 7, which is
`support` — a marker count, not an abundance.

`scripts/profile_saliva.py` must be re-run on the saliva reads (not in this repo). The scripts
now refuse to run against the stale table rather than silently substituting a column.

**Any number in this repo or in the manuscript that came from `saliva_permanova.tsv`,
`saliva_strain_long.tsv` or the concordance/temporal analyses is stale until that happens** —
including the "strain R² 0.833 > species 0.822" claim, which reverses on the only abundance the
stale table contains (0.6886 vs 0.7101 on `within_abund`). Note this is the one place where §1
also bites: `support` is precisely the column that moves ~6% with the engine change.

## 3. `--marker-source kmer` is undocumented in the manuscript

The engine gained a FracMinHash-style k-mer marker source (`63c4f50`). Methods describes only
enzyme digestion. See `methods_algorithm.md`.
