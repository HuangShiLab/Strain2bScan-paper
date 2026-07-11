# Motivation figure — 2bRAD captures strain-level signal that 16S cannot (15 species)

The "why 2bRAD for strains" motivation (StrainScan's Fig 1 analog), recomputed with Strain2bScan
across all **15 species** from the Strain2bfunc slide deck (PPTX slide 23). Tool-agnostic property of
the 2bRAD representation vs 16S; reproduces the prior finding.

## Method — `scripts/compute_16s_concordance.py`, `scripts/plot_16s_motivation.py`
For each species, ~50 genomes (NCBI accession lists → ENA FASTA; `scripts/dl_16s_panel.py`). Three
between-strain distance spaces, all with the Mash transform D(J) = −ln(2J/(1+J)):
- **whole-genome**: bottom-3000 canonical 21-mer MinHash (pure-Python Mash).
- **2bRAD**: Strain2bScan `build` BcgI tag sets (per-genome), tag-set Jaccard.
- **16S**: longest 16S gene per genome extracted with **barrnap 0.9** (HMMER/nhmmer), exact
  canonical 21-mer Jaccard.

Per species, the pairwise distance vector for 2bRAD and for 16S is correlated (Spearman) against the
whole-genome vector. Run with `PYTHONHASHSEED=0`.

## Result — `results/mash_2brad_vs_16s_recomputed.tsv`, `figures/mash_2brad_vs_16s.*`

**2bRAD median Spearman 0.90 vs 16S median 0.36** across the 15 species. 2bRAD between-strain
distances track whole-genome distance in every species (0.59–0.99); 16S is weak to useless
(−0.14 to 0.78). Signature contrasts (2bRAD vs 16S):

| species | 2bRAD | 16S | | species | 2bRAD | 16S |
|---|---|---|---|---|---|---|
| *M. tuberculosis* | 0.895 | **−0.135** | | *E. coli* | 0.989 | 0.360 |
| *S. aureus* | 0.983 | 0.086 | | *C. acnes* | 0.977 | 0.470 |
| *L. plantarum* | 0.867 | 0.023 | | *S. enterica* | 0.976 | 0.608 |
| *K. pneumoniae* | 0.734 | 0.239 | | *P. copri* | 0.976 | 0.370 |
| *S. epidermidis* | 0.900 | 0.205 | | *A. muciniphila* | 0.958 | 0.648 |
| *S. pneumoniae* | 0.589 | 0.203 | | *C. difficile* | 0.942 | 0.356 |
| *P. dorei* | 0.896 | 0.179 | | *P. gingivalis* | 0.865 | 0.419 |
| | | | | *F. nucleatum* | 0.813 | 0.782 |

The extreme case is *M. tuberculosis* — a monomorphic pathogen whose 16S is essentially invariant
across strains (Spearman ≈ 0, even slightly negative), while its 2bRAD tags still recover the
whole-genome strain structure (0.90). *S. aureus* (0.98 vs 0.09) and *L. plantarum* (0.87 vs 0.02)
are similar. Only *F. nucleatum* has appreciable 16S signal (0.78), and even there 2bRAD matches it.

Reproduces the Strain2bfunc slide-deck values (2bRAD median ~0.90; 16S median low) with Strain2bScan
tags + an open 16S pipeline. This is the Introduction/Fig-2 motivation: **16S resolves species, not
strains; 2bRAD's reduced tags carry genome-wide strain signal at ~1% of the sequencing.**
