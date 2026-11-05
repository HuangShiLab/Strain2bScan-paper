# Strain-sharing analysis

Inputs are already-gated Strain2bScan 0.95 strain-cluster calls.

For each sample pair and each species detected in both samples:
- `any_cluster`: at least one identical cluster label was observed in both samples.
- `dominant_cluster`: the highest within-species abundance cluster was identical in both samples.
- The pair-level primary sharing rate is `n_shared_species_dominant_cluster / n_analyzable_species`.
- A coexistence-aware secondary rate is `n_shared_species_any_cluster / n_analyzable_species`.
- Per-species cluster-set Jaccard is also reported for multi-strain species.

Saliva comparisons are within-subject timepoint contrasts and between-subject controls.
Monkey IDs beginning with O and Y denote old and young cohorts; GI region and compartment are used for
corresponding-site contrasts.

Reproduce with:

```bash
python3 scripts/compute_strain_sharing.py
python3 scripts/bootstrap_strain_sharing.py
/usr/bin/python3 scripts/plot_strain_sharing.py
```
