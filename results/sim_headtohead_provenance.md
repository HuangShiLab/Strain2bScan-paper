# 15-species head-to-head provenance status

Date: 2026-10-07

The original `scratchpad/eval` drivers and raw StrainScan outputs were not retained. That limitation is now
addressed by a compact reproducible StrainScan v1.0.14 rerun. `results/strainscan_rerun/` contains the
run manifest, database-build logs, StrainScan final reports, failed/no-call list, run-level scores, median
aggregates, multi-species coverage audit, and SHA-256 checksums. The rerun paired 204 completed
single-species samples across 14 species from 225 intended matched entries; 21 intended entries ended
without calls, including all *Salmonella enterica* entries. Twelve multi-species samples were rerun, with
158 of 180 possible species-by-sample reports; absent reports were treated as no detection.

Accuracy and cost are deliberately kept in separate evidence classes. Accuracy in Figure 11A–C and Tables
1–2 uses the rerun. Database-build, same-environment profiling, and multi-species timing panels in Figure
11D–F use the archived timing/build tables because those runs included controlled cost measurements. These
classes are not pooled. The archived timing run did not complete *Klebsiella pneumoniae* under StrainScan,
whereas the reproducible accuracy rerun did complete it; both facts are reported rather than replaced.

`results/sim_headtohead_provenance.tsv` records current evidence files, SHA-256 values, row counts, roles,
and status. `results/strainscan_rerun/rerun_provenance.json` records the rerun scope, container metadata,
bootstrap intervals, and analysis entry points. `results/benchmark_configuration.json` records the frozen
primary ATCC mock configuration and the distinction between primary mock and comparator-rerun evidence.
