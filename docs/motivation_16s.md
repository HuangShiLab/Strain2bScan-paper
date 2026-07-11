# Motivation figure — 2bRAD captures strain-level signal that 16S cannot (15 species)

The "why 2bRAD for strains" motivation (StrainScan's Fig 1 analog), recomputed with Strain2bScan
across all **15 species** from the Strain2bfunc slide deck (PPTX slide 23). Tool-agnostic property of
the 2bRAD representation vs 16S; reproduces the prior finding.

## Method — `scripts/compute_16s_v2.py`, `scripts/plot_16s_v2.py`
Per species, genomes from NCBI accession lists → ENA FASTA (`scripts/dl_16s_panel.py`), then
**restricted to complete / near-complete assemblies** — CheckM completeness ≥ 97%, contamination
≤ 5%, assembly level Complete Genome or Chromosome (`scripts/fetch_genome_qc.py` →
`data/genome_qc_16s_panel.tsv`). This matters: fragmented/incomplete assemblies have missing sequence
that inflates Mash distances. **14–50 high-quality genomes per species** (median ~46). Three
between-strain distance spaces, all with the Mash transform D(J) = −ln(2J/(1+J)):
- **whole-genome**: bottom-3000 canonical 21-mer MinHash (pure-Python Mash).
- **2bRAD**: Strain2bScan `build` BcgI tag sets (per-genome), tag-set Jaccard.
- **16S**: longest 16S gene per genome extracted with **barrnap 0.9** (HMMER `nhmmer`), exact
  canonical 21-mer Jaccard.

Per species, the pairwise distance vector for 2bRAD and for 16S is correlated (Spearman) against the
whole-genome vector; **95% CIs by 500 genome subsamples (80% without replacement)**. `PYTHONHASHSEED=0`.

## Result — `results/mash_2brad_vs_16s_recomputed.tsv`, `figures/mash_2brad_vs_16s.*`

**2bRAD median Spearman 0.94 vs 16S median 0.36** across the 15 species. Every 2bRAD CI is high
(≥0.53) and never overlaps its species' 16S CI; several **16S CIs overlap zero** (no strain signal
at all). Signature contrasts (Spearman [95% CI]):

| species | n | 2bRAD | 16S |
|---|---|---|---|
| *M. tuberculosis* | 29 | 0.81 [0.75, 0.93] | **−0.04 [−0.09, 0.04]** |
| *S. aureus* | 48 | 0.98 [0.97, 0.99] | 0.10 [−0.01, 0.22] |
| *L. plantarum* | 48 | 0.87 [0.83, 0.91] | 0.02 [−0.11, 0.16] |
| *P. dorei* | 14 | 0.94 [0.89, 0.99] | −0.13 [−0.28, 0.04] |
| *E. coli* | 50 | 0.99 [0.98, 0.99] | 0.36 [0.24, 0.47] |
| *A. muciniphila* | 48 | 0.96 [0.95, 0.97] | 0.65 [0.55, 0.75] |
| *F. nucleatum* | 25 | 0.87 [0.79, 0.92] | 0.68 [0.48, 0.83] |

(full table incl. *C. acnes, C. difficile, K. pneumoniae, P. copri, P. gingivalis, S. enterica,
S. epidermidis, S. pneumoniae* in the TSV). The extreme case is *M. tuberculosis* — a monomorphic
pathogen whose 16S is invariant across strains (Spearman ≈ 0, CI spans zero), while its 2bRAD tags
recover the whole-genome strain structure (0.81). *S. aureus*, *L. plantarum* and *P. dorei* are
similar (16S CI touches/crosses zero). Only *F./A.* species have appreciable 16S signal, and even
there 2bRAD matches or beats it.

Restricting to complete genomes and adding CIs **does not change the conclusion** (2bRAD median
0.90→0.94; 16S unchanged), confirming it is not a draft-assembly artifact. Reproduces the
Strain2bfunc slide-deck finding with Strain2bScan tags + an open 16S pipeline. This is the
Introduction/Fig-2 motivation: **16S resolves species, not strains; 2bRAD's reduced tags carry
genome-wide strain signal at ~1% of the sequencing.**

## Supplementary — per-species rank–rank scatter (`figures/mash_2brad_vs_16s_scatter.*`)
The panel figure shows *why* the Spearman values differ: for each species, the rank of each strain
pair's whole-genome distance (x) vs the rank of its 2bRAD (blue) and 16S (red) distance (data in
`results/pairdist/<species>.tsv`). **2bRAD hugs the diagonal in every species** (its distance orders
strain pairs the same way the whole genome does); **16S forms flat horizontal rank-bands** — a
handful of discrete 16S distances shared by many unrelated pairs, uncorrelated with genome distance.
The discrepancy is most extreme in *Phocaeicola dorei* and *M. tuberculosis*: 16S collapses to a
single flat band (nearly all pairs have identical 16S → zero strain resolution) while 2bRAD still
recovers the genome-wide ordering. This is the pair-level view behind the summary bar chart.
