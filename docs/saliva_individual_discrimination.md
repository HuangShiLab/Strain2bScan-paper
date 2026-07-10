# Saliva case study (Fig 10) — strain-level individual discrimination

Real-data case study closing the "all-simulated" gap: does **strain-level** 2bRAD profiling
distinguish individuals better than species-level? Reproduces the prior Strain2bfunc result with the
fast/light Strain2bScan.

## Data
SRA **PRJNA1131785** (npj Biofilms Microbiomes 2025), native BcgI 2bRAD saliva, the diurnal design:
**8 subjects × 4 timepoints = 32 samples** (R1 used). Sample aliases decode as `S<time>-<subject>`:
prefix `S9/S11/S13/S17` = time of day (9 AM / 11 AM / 1 PM / 5 PM), suffix `1,2,6,7,9,11,12,14` = the
8 subjects. (Confirmed empirically: PERMANOVA subject R²≈0.83 vs timepoint R²≈0.05.)

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
relative-abundance matrices (cluster support, per-sample normalized) → Bray–Curtis → PERMANOVA
(adonis, 4999 perms) + leave-one-out 1-NN subject classification. `scripts/profile_saliva.py`,
`scripts/analyze_saliva.py`.

## Result — `results/saliva_permanova.tsv`, `results/saliva_perspecies_subject.tsv`, `figures/saliva_individual_discrimination.*`

| factor | level | R² | p | LOO 1-NN acc |
|---|---|---|---|---|
| **subject** | species | 0.822 | 0.0002 | 87.5% |
| **subject** | **strain** | **0.833** | 0.0002 | **90.6%** |
| timepoint | species | 0.071 | 0.67 | 9% |
| timepoint | strain | 0.053 | 0.89 | 9% |

- **Individuals are strongly separable, and strain-level beats species-level** on both R²
  (0.833 > 0.822) and 1-NN classification (90.6% > 87.5%) — reproducing Strain2bfunc's finding
  (strain > species for host discrimination) with a tool that runs in ~1 s / ~80 MB instead of
  ~40 min / ~30 GB.
- **Time of day has no detectable effect** (R² ≈ 0.05, ns): an individual's oral strain profile is
  stable across the day, so the separation is truly by person, not batch/collection time.
- **Per-species strain-level subject R²**: *Rothia mucilaginosa* **0.921** — the exact species
  Strain2bfunc highlighted (they reported 0.87) — then *Neisseria subflava* 0.885,
  *Streptococcus salivarius* 0.850, *Gemella haemolysans* 0.780, *Prevotella melaninogenica* 0.770…
  **17 of 18 testable species significant** (p < 0.05); only *S. mitis* not.

This is the biological headline: on real, error-containing native-2bRAD saliva, Strain2bScan resolves
person-specific strain signatures, with *Rothia mucilaginosa* the strongest individual marker —
validating the tool on real open-world data and reproducing the prior result at a fraction of the compute.

## Open extensions
- Paired shotgun↔2bRAD concordance (WMS saliva exists in PRJNA1131785; large — grab a few).
- Temporal-stability quantification (within-subject across the 4 timepoints; ICC).
- ML host-ID accuracy on a held-out timepoint (Strain2bfunc reported 100%).
