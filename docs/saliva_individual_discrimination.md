# Saliva case study (Fig 7) — strain-level individual discrimination

Real-data case study: does **strain-level** native BcgI 2bRAD-M profiling distinguish individuals
better than species-level? The numbers below are the corrected `sample_fraction` values used in
the manuscript.

## Data
SRA **PRJNA1131785** (npj Biofilms Microbiomes 2025), native BcgI 2bRAD-M saliva, the diurnal design:
**8 subjects × 4 timepoints = 32 samples** (R1 used). Sample aliases decode as `S<time>-<subject>`:
prefix `S9/S11/S13/S17` = time of day (9 AM / 11 AM / 1 PM / 5 PM), suffix `1,2,6,7,9,11,12,14` = the
8 subjects. (Confirmed empirically: PERMANOVA subject R² 0.757 vs timepoint R² 0.038.)

## Reference panel
A **saliva/oral-commensal panel** — 19 abundant oral species, up to 25 genomes each (359 total,
NCBI accession lists → ENA FASTA; `scripts/dl_oral_panel.py`), clustered with
`cluster --enzyme BcgI --similarity 0.95`. Many genomes/species is essential: it lets clustering
form real strain clusters that individuals differ on. (A first attempt with the generic 62-species
*pathogen* panel gave a null result — 4 genomes/species and no oral commensals, so real saliva
strains map ~uniformly across arbitrary clusters. Panel choice is the key methodological point.)

Species: *Rothia mucilaginosa, R. dentocariosa, Streptococcus salivarius/mitis/oralis/sanguinis/
parasanguinis/gordonii, Haemophilus parainfluenzae, Veillonella parvula/atypica/dispar, Prevotella
melaninogenica, Neisseria subflava, Actinomyces odontolyticus, Gemella haemolysans, Granulicatella
adiacens, Porphyromonas gingivalis, Fusobacterium nucleatum.*

## Analysis
`multi-profile --enzyme BcgI --min-species-markers 50 --min-species-detect 5` on all 32 samples
(each ~4–18 s; up to 18 species resolved, 40–188 strain calls/sample). Strain- and species-level
relative-abundance matrices (`sample_fraction`, per-sample normalized) → Bray–Curtis → PERMANOVA
(adonis, 4999 perms) + leave-one-out 1-NN subject classification. `scripts/profile_saliva.py`,
`scripts/analyze_saliva.py`.

## Result — `results/saliva_permanova.tsv`, `results/saliva_perspecies_subject.tsv`, `figures/saliva_individual_discrimination.*`

| factor | level | R² | p | LOO 1-NN acc |
|---|---|---|---|---|
| **subject** | species | 0.755 | 2e-4 | 78.1% |
| **subject** | **strain** | **0.757** | 2e-4 | **100%** |
| timepoint | species | 0.038 | 0.99 | — |
| timepoint | strain | 0.038 | 0.99 | — |

- **Individuals are strongly separable, and strain-level identification is perfect**: leave-one-
  timepoint-out nearest-neighbour classification identified the host with **100% accuracy from strain
  features versus 78.1% from species features** (both PERMANOVA R² ≈ 0.76, p = 2 × 10⁻⁴). Strain and
  species PERMANOVA R² are comparable (0.757 versus 0.755), so the advantage is strongest for subject
  identification and low-abundance strain recovery.
- **Time of day has no detectable effect** (R² = 0.038, p = 0.99): an individual's oral strain profile
  is stable across the day, so the separation is truly by person, not batch/collection time.
- **Per-species strain-level subject R²**: the strongest single-species marker is ***Neisseria
  subflava*** (R² = 0.808); all 13 testable species were significant (p < 0.05).

This is the biological headline: on real, error-containing native BcgI 2bRAD-M saliva, Strain2bScan
resolves person-specific strain signatures that identify the host perfectly, validating the tool on
real open-world data at a fraction of the compute of full-k-mer methods.

## Open extensions
- Paired shotgun↔2bRAD concordance (WMS saliva exists in PRJNA1131785; large — grab a few).
- Temporal-stability quantification (within-subject across the 4 timepoints; ICC).
- ML host-ID accuracy on a held-out timepoint (Strain2bfunc reported 100%).
