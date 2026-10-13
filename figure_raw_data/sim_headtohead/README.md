# Strain2bScan vs StrainScan on simulated metagenomes

Both tools use the same 15-species reference pool and simulated reads under `../sim_single_species/`
and `../sim_multi_species/`. Accuracy is scored separately in each tool's own 0.95-cluster space.
Strain2bScan ran natively on Apple Silicon. StrainScan is Linux-x86-only and ran under Docker
`linux/amd64` QEMU emulation, so its wall-clock is inflated; a same-environment Strain2bScan subset is
provided for profiling-speed comparison.

## Current accuracy evidence

The current Strain2bScan regression used source commit `6ae760c2f3bd2f39efc0e6aa06c53a84366cb16e`
(tag `benchmark-v0.1.1-6ae760c`) with current depth-adaptive defaults. It reran the 225 matched
single-species entries and the 12 matched multi-species communities. Aggregate accuracy is in:

- `../../results/strainscan_rerun/current_single_by_species.tsv`
- `../../results/strainscan_rerun/current_single_by_depth.tsv`
- `../../results/strainscan_rerun/current_multi_by_depth.tsv`
- `../../results/strainscan_rerun/current_uncertainty_summary.tsv`

The paired point estimates and caveats are reported in Table 1–3 of the manuscript. The important
trade-off is lower current Strain2bScan precision/recall on *K. pneumoniae* but higher recall on most
other species.

## Files

- Historical Strain2bScan timing and 2025 source tables remain in this directory.
- `../../results/strainscan_rerun/current_single_persample.tsv` and
  `../../results/strainscan_rerun/current_multi_persample.tsv` supersede the corresponding 2025
  Strain2bScan accuracy columns.
- `sim_headtohead_build_cost.tsv` and `profile_timing_sameenv.tsv` remain the archived timing evidence.
- Archived timing aggregates with a `superseded for accuracy` role are indexed in
  `../../results/sim_headtohead_provenance.tsv`.

## Reproduce

From the paper repository:

```bash
python3 scripts/regression_strain2bscan_current.py
python3 scripts/aggregate_strain2bscan_regression.py
```

The regression script defaults to the sibling Strain2bScan repository, rebuilds species databases,
reruns the matched reads, scores cluster sets, and compares against the reproducible StrainScan table.
