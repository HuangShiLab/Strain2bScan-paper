# Strain2bScan mock diagnostic stdout (2026-09-26)

Generated with local Strain2bScan binary built from `strainscan-port` branch
(includes commit `bb16eb4 Recalibrate the abundance floor to the unbiased depth estimator`).

## 99.9 % host sample (MSA-1002)

Sample: `BcgI_MSA1002_99.9_100ng_1` (R1+R2 merged).

| File | Command gist | Outcome |
|---|---|---|
| `99.9host_combined_default.stdout` | `profile --db MSA_combined164_bcgi_cont.tsv --reads ...` | Only **C42 (E. coli)** called |
| `99.9host_combined_threshold.stdout` | `+ --min-abundance 0.02` | Same as default |
| `99.9host_combined_oldlike.stdout` | `+ --min-abundance 0.02 --fixed-gate` | Same as default |
| `99.9host_combined_minconsistency0.stdout` | `+ --min-consistency 0` | **14 clusters** called |
| `99.9host_Escherichia_coli_default.stdout` | per-species E. coli DB | C0 called (depth ~50x) |
| `99.9host_Escherichia_coli_oldlike.stdout` | per-species E. coli DB, old-like | Same |

Interpretation: the 99.9 % host "recall collapse" is driven by the **consistency gate**
(`--min-consistency`, default 0.5), not by `--min-abundance`. Relaxing consistency to 0
recovers 13 additional clusters.

## MSA-1003 staggered mock

Sample: `BcgI_MSA1003_0_100ng_1` (R1+R2 merged).

| File | Outcome |
|---|---|
| `MSA1003_default.stdout` | `multi-profile` default: 18/20 species strain-resolved, **14 strain calls** |

MSA1003 is robust to the current default gating.
