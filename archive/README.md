# Archive

This directory separates historical and bulk material from the active manuscript.

## Contents

- `legacy_docs/`: superseded planning, session-handoff, and earlier review records.
- `legacy_figures/`: superseded exploratory figures retained for provenance.
- `legacy_results/`: diagnostic and intermediate tables that are not current manuscript evidence.
- `large_intermediates/`: local bulk reads and regenerated simulation inputs (ignored by git).
- `work_legacy/`: earlier working directories and HPC monitoring artifacts (ignored by git).

The active manuscript sources, final numbered figures, source panels, current result tables, and
frozen provenance remain under `manuscript/`, `figures/numbered/`, `figures/source_panels/`,
`results/`, `data/`, and `scripts/`.

Bulk simulation reads are reachable through compatibility symlinks in `figure_raw_data/`, so the
committed manifests and rerun scripts continue to resolve the same paths without duplicating data.
