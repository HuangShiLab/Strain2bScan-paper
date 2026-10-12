# Figure layout

- `numbered/` is the publication figure set. It contains final Figures 1–11, Supplementary Figures
  S1–S5, and `MANIFEST.tsv` with PNG/PDF SHA-256 checksums.
- `source_panels/` contains only the current panel images used to build `numbered/`. Regeneration
  scripts write here; `scripts/number_figures.py` writes the final set to `numbered/`.

Superseded figures are in `../archive/legacy_figures/`.
