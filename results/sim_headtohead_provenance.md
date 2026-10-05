# 15-species head-to-head provenance status

Date: 2026-10-06

The paper repository retains complete Strain2bScan native and same-container per-sample tables, aggregate
head-to-head tables, build-cost tables, and the timing subset. The original `scratchpad/eval` drivers and the
complete raw StrainScan per-sample outputs were not retained. Consequently, the 204-sample StrainScan
comparison is reproducible from the committed aggregate evidence only at the reporting level; the underlying
run-level StrainScan artifacts are incomplete.

`results/sim_headtohead_provenance.tsv` records the available files, SHA-256 values, row counts, expected
role, and provenance status. `results/benchmark_configuration.json` records the frozen primary ATCC mock
configuration. The manuscript states this limitation rather than claiming complete run-level reproducibility.

Restoration requires one of the following:

1. archive the original StrainScan container, commands, cluster mappings, and complete per-sample outputs;
2. rerun the 204 paired single-species samples and 12 multi-species samples with a pinned StrainScan image;
3. replace the current head-to-head with a newly archived, deterministic benchmark at one immutable commit.
