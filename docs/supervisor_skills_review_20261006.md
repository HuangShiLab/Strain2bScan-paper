# Supervisor-Skills second review, Strain2bScan manuscript

Review date: 2026-10-06
Commit reviewed: `be06b2d6d217c7c207e5ea8bbf3043e03e916780`
Primary manuscript: `manuscript/full_manuscript.md`
Skills applied: `benchmark-paper-template`, `pre-submission-reviewer`, and `paper-polish` principles.
Automated gates: `scripts/check_manuscript_claims.py` passed; all 16 numbered figures pass `figures/numbered/MANIFEST.tsv` checksum checks.

## Executive verdict

**Major revision, but substantially closer to submission.** The primary ATCC mock configuration is now frozen, Figures 6 and 12 use primary rows only, projected runtime is labelled, and the saliva and cohort claims are much more cautiously framed. The remaining blockers are no longer mainly claim tone. They are structural assembly quality, a mismatch between the benchmark-frozen implementation and the implementation described in the Methods/Introduction, incomplete run-level comparator artifacts, and missing inferential uncertainty for the central benchmark and real-data claims.

## Scores

- Benchmark-paper five-pillar readiness: 6.5 / 10
- Pre-submission readiness: 5 / 10
- Overall recommendation: **needs major revision before submission**

## Five-pillar audit

| Pillar | Assessment | Main remaining issue |
|---|---|---|
| Research gap | Strong. The scale barrier of per-species k-mer profilers and the low-biomass/high-host barrier of shotgun are clear and well paired with native 2bRAD-M. | The Introduction should explicitly state the evaluation question as a benchmark contribution rather than only a method advantage. |
| Construction pipeline | Much stronger after freezing `f26f234`, archiving primary mock predictions, database hashes, and configuration. | The 15-species StrainScan comparison remains aggregate-level because original run drivers and complete raw StrainScan outputs are unavailable. |
| Evaluation framework | Broad and multidimensional: reference completeness, enzyme number, host fraction, DNA input, mocks, saliva, public WGS subsets, runtime, and memory. | Primary metrics are now clearer, but the manuscript still needs prespecified paired tests, confidence intervals, and multiplicity handling. |
| Empirical findings | More credible after revising clean-mock precision and cohort language. The host-contamination result is now bounded and primary. | The projected 55-species runtime result should be either measured natively or moved from headline framing to clearly secondary evidence. |
| Companion method | Appropriate: Strain2bScan is the implementation being evaluated. | The benchmark-frozen software needs a public tagged release or archived source snapshot separate from the moving development branch. |

## Critical issues

### C1. Assembled manuscript has a duplicated Results heading

The assembled manuscript contains both:

- line 97: `## Results`
- line 99: `### Results`

This comes from `manuscript/results.md` beginning with `## Results`, while `scripts/assemble_manuscript.py` also prepends `## Results`. The result is a structurally malformed main manuscript.

**Fix.** Remove the top-level `## Results` heading from `manuscript/results.md`, or make `results()` strip it before assembly. Add an assembly freshness and heading-collision check to `scripts/check_manuscript_claims.py`.

### C2. The described software version and the benchmark-frozen version are not yet reconciled

The primary benchmark is frozen at Strain2bScan commit `f26f234b817ba7772a3f1df59ce720751e9b45b9`, with binary SHA-256 `a4cf7a4043a06c55a99d6abf1e9fd312f6243aa5a5ca021bd41915636fb56b18`. That is good. However, the current `strainscan-port` branch has moved to `967a8f5`, which uses Clap as a dependency. The Introduction and Methods still describe the tool as “dependency-free Rust” and “no third-party dependencies”.

This creates two different software identities:

1. the frozen benchmark implementation used for the paper;
2. the current moving implementation described in the text.

**Fix.** Either:

- tag and archive `f26f234` as the benchmark release and state explicitly that the benchmark binary was dependency-free at that commit, while the current development branch uses Clap; or
- rerun the primary benchmark with the new tagged release and update all hashes and commands.

For submission, create a Zenodo or GitHub release snapshot containing source, binary build recipe, tests, and the binary hash.

### C3. The central 15-species StrainScan comparison still lacks complete run-level evidence

The new provenance files honestly state that only 3 raw StrainScan single-species rows and 3 multi-species rows are available, versus 204 and 12 required. The original `scratchpad/eval` drivers were not retained. Therefore Tables 1–3, Figure 11, and the abstract recall comparison cannot yet be independently reproduced at run level.

**Fix.** Choose one:

1. recover the original StrainScan container, commands, cluster mappings, outputs, and analysis drivers;
2. rerun the 204 paired single-species samples and 12 multi-species samples with a pinned StrainScan image;
3. downgrade the 15-species comparison to an aggregate exploratory result and make the ATCC primary benchmark the central accuracy claim.

The third option is weakest for a benchmark paper but is currently the only fully supported option without new runs.

## Major issues

### M1. The strongest scalability number remains projected

The Abstract and Introduction report a projected 121–146× advantage on 55-species communities. The manuscript now labels this correctly, which is a real improvement, but it remains a projected comparator rather than a measured native benchmark.

**Fix.** Rerun StrainScan natively for at least 10, 100, and 200 samples, or move the projected 55-species result after the measured single-panel and same-container comparisons and avoid using it as the main scalability claim.

### M2. Benchmark accuracy still lacks uncertainty

Tables 1–3 report medians without confidence intervals. There are no paired bootstrap or permutation intervals for recall differences, no classification confidence interval for the 100% saliva result, no sample-aware model for repeated strain calls, and no multiplicity adjustment for species-level saliva tests.

**Fix.** Add paired bootstrap or permutation intervals for primary benchmark differences, binomial or cross-validated intervals for host-ID accuracy, sample-aware inference for call-level abundance comparisons, and explicit multiplicity handling.

### M3. Introduction still contains overclaims that Results now avoids

Examples:

- “for the first time, strain-level analysis” at Introduction line 42;
- “The two modes are shown to agree ... recover the same strains” at Introduction lines 55–56;
- “dependency-free Rust” at line 77;
- “uniquely accepts two input modes” at line 78.

These are stronger than the current qualified Results and Discussion.

**Fix.** Replace “for the first time” with “to our knowledge” plus the precise setting, or remove it. Replace “recover the same strains” with “65/65 shotgun calls in three usable paired samples were also present in native calls.” Specify dependency-free only for the frozen benchmark binary.

### M4. Assembled figure legends contain broken file-reference artifacts

After assembly, several legends have empty references or stray punctuation:

- Figure 3: `(A) ;` and `(B) ;`
- Figure 9: `**(B)** : thread scaling`
- Figure 10: `**Figure 10 ...** :`
- Figure S2: `**Figure S2 ...** ;`

The assembler removes file paths but leaves orphan punctuation. Figure 9 also still needs visible A/B/C panel labels in the rendered figure.

**Fix.** Make `figure_legends()` preserve source-file names or remove the trailing punctuation cleanly. Add visible panel letters to Figure 9 and rewrite its caption so every rendered panel is specified.

### M5. References are not integrated into the text

The bibliography is still called a “Working bibliography”, and the author-facing note says DOIs and publication details must be checked. Tools named in the text, including StrainScan, StrainGE, sylph, inStrain, SPAdes, ART, and PERMANOVA, are not consistently cited at their use sites.

**Fix.** Convert to a citation manager or BibTeX workflow, add in-text citations at first mention, verify every DOI and year, remove the author-facing note, and add a Related Work or benchmark-comparison paragraph/table if the target journal expects it.

### M6. Public cohort reporting still lacks complete run-level metadata

The cohort results are now correctly described as exploratory, but Table 4 does not include SRR accessions, read counts, estimated depth, body site, or selection rationale for every library. Human-data provenance and data-use statements are also not complete.

**Fix.** Add a supplementary cohort manifest with sample/run accession, body site, timepoint, read pairs, host fraction if available, panel, command, and output checksum. Add explicit data-use and ethics statements.

### M7. Runtime reporting lacks a full environment lock

The frozen configuration includes software commit, binary hash, Rust version, and database hashes, which is strong. But the manuscript does not yet report CPU model, OS/kernel, memory configuration, I/O state, thread pinning, container image digest for StrainScan, number of timing repetitions, and variance.

**Fix.** Add a machine-readable environment manifest for every timing benchmark. For StrainScan, include container image ID/digest, command, versions, and host resource limits.

## Minor issues

1. The abstract title is “strain-level profiling from 2bRAD-reduced markers”, while the assembled title is “strain-level metagenomic profiling on 2bRAD-reduced markers”. Choose one.
2. The Introduction overuses semicolons where commas, colons, or separate sentences are clearer. For example, lines 3–4 separate phenotype examples with semicolons.
3. Supplementary material lists S3 before S1 and S2. Reorder or explain the ordering.
4. Draft wording remains in modular headings such as “Introduction (draft)” and “Discussion (draft)”.
5. “Part 0” and “Part I/II” are unconventional for a journal manuscript. Consider “Overview”, “Native 2bRAD-M”, and “Community-scale shotgun profiling”.
6. Several Results paragraphs still mix thresholded and threshold-free metrics in one sentence. Keep `AUPR`, `P@1e-4`, `R@1e-4`, and `F1@1e-4` visually separate in text and figures.
7. Some figure legends still say “strain” where the measured unit is a 0.95-similarity cluster or merged panel unit.

## What improved since the previous review

1. The primary ATCC benchmark is now frozen at source, binary, database, metric, and per-output levels.
2. Figures 6 and 12 no longer mix primary and sensitivity variants.
3. Clean-mock precision trade-offs are now reported rather than hidden.
4. Runtime projections are labelled as projections.
5. Saliva concordance is separated from independent validation.
6. Cohort isolate recovery is called a panel/read compatibility check.
7. The main text is much shorter and easier to follow.
8. Automated claim and figure checksum gates now exist.

## Immediate next goal

The next goal should be a submission-ready benchmark release, not further algorithm changes.

Recommended order:

1. Fix assembly and create a release candidate: `v0.2.0-benchmark.1`.
2. Tag/archive the benchmark software commit and current paper commit.
3. Recover or rerun the 204-sample StrainScan comparison.
4. Add uncertainty and multiplicity reporting.
5. Complete author, data-use, accession, and reference metadata.
6. Regenerate figures with panel labels and clean legends.
7. Run the claim checker, assembly freshness check, figure checksum check, and human figure review.
8. Freeze a DOCX/PDF submission package and inspect every page.

## Revised target after fixes

If C1 through C3 and M1 through M3 are resolved, the paper should be materially stronger. The likely submission story is then:

> Strain2bScan makes strain-resolved profiling practical in native low-biomass/high-host 2bRAD-M data and scales digest-once shotgun profiling to many samples, with a frozen benchmark and transparent operating limits.

That is a stronger and safer claim than asserting universal clean-mock superiority or cohort-level biological discovery.
