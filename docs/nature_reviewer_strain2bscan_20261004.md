# Nature-reviewer assessment, Strain2bScan manuscript

Review date: 2026-10-04
Repository: `/Users/macstudio/Downloads/Strain2bScan-paper`
Commit under review: `b3161a1`
Primary manuscript: `manuscript/full_manuscript.md`
Method: three isolated reviewer passes followed by one post-review synthesis. No reviewer saw another report before freezing its own assessment.

## Review setup

### Input scope

The review covered the assembled manuscript at commit `b3161a1`, the modular manuscript sections, the figure map and detailed legends, the numbered figure manifest, tracked result tables, and result files directly cited by the manuscript. The primary manuscript is 1223 lines and approximately 14,610 words. All 16 numbered figure entries are present as PNG and PDF and pass the SHA-256 checks in `figures/numbered/MANIFEST.tsv`. The text references Figures 1 through 12 and S1 through S4.

### Assessment boundary

This is a pre-submission review, not an editorial decision. No external literature search was performed. Claims were assessed against the repository evidence available at the stated commit. Scratch and untracked files were treated as outside the immutable evidence base unless directly cited by the manuscript. Missing material is labelled as not assessable rather than treated as an automatic flaw.

### Shared manuscript claim summary

The manuscript claims that Strain2bScan enables strain-resolved profiling on native BcgI 2bRAD-M libraries and in-silico-digested shotgun metagenomes. It reports controlled simulations, reference-completeness experiments, ATCC mock communities, real saliva analyses, public WGS cohort applications, isolate-derived panel compatibility checks, and runtime and memory comparisons.

### Missing materials affecting confidence

The following directly cited or central benchmark materials were not complete at the reviewed commit:

- A pinned Strain2bScan source release, binary checksum, compiler version, and per-benchmark binary identity.
- Complete raw per-sample StrainScan outputs for the 204 single-species and 4-per-depth multi-species comparisons.
- The cited `scratchpad/eval` analysis drivers.
- Complete machine-independent scripts and environment locks for every figure and benchmark.
- Complete run-level and panel-accession manifests for all public and reference-panel analyses.

---

## Reviewer 1

### Overall assessment

The manuscript presents a potentially useful dual-mode strain profiler for native 2bRAD and in-silico-digested shotgun data. The concept is clearly motivated, the evaluation spans simulation, genome completeness, mock communities, saliva, public cohorts, and runtime benchmarks, and the repository contains many aggregated result tables and checksummed figures.

However, the central benchmark and accuracy case is not established in the current commit. There are direct contradictions between manuscript claims and tracked result tables, missing per-sample evidence for key comparator runs, projected performance claims that are not consistently labelled, weak treatment of dependence and uncertainty, and reproducibility gaps. These issues affect headline claims in the Abstract and Results.

### Who would be interested in the results, and why

Microbiome researchers, microbial population geneticists, clinical metagenomics groups, and bioinformatics method developers would be interested if the claims were supported. Native 2bRAD profiling could reduce host-dominated sequencing waste while preserving strain-level resolution. Digest-once shotgun profiling could reduce the cost of multi-species strain searches.

### Major strengths

1. The dual-input design is well motivated and clearly distinguishes native BcgI libraries from in-silico shotgun digestion.
2. The 2bRAD versus 16S motivation analysis is supported by the tracked per-species correlation table.
3. The benchmark scope is broad and spans completeness, enzyme number, host fraction, DNA input, mocks, saliva, public cohorts, and scaling.
4. The manuscript candidly acknowledges cluster-level resolution, closed panels, application-only cohort subsets, and nonindependent isolate self-recovery.
5. The figure provenance system is unusually useful.

### Major Concerns

#### R1-M1 [experimental-design]

**Severity** Major
**Blocking** Yes
**Claim pointer** The Methods describe the systematic benchmark reads as 150 bp paired-end ART reads and elsewhere as error-free 150 bp reads.
**Evidence pointer** `scripts/simulate_single_species.py`, `scripts/simulate_multi_species.py`, `figure_raw_data/sim_README.md`, and `docs/simulation_audit.md` describe 250 bp ART reads with an error model.
**Concern** The stated simulator specification does not match the tracked scripts and documentation.
**Why it matters** Read length and error model affect marker recovery, singleton policy, detection thresholds, and comparator performance.
**Resolution test** State the exact simulator, read length, insert size, quality model, seeds, versions, and binary provenance. Audit or regenerate all reads against one fixed specification.

#### R1-M2 [reproducibility]

**Severity** Major
**Blocking** Yes
**Claim pointer** The manuscript claims a 204-sample paired StrainScan comparison and says drivers and raw per-sample tables are available.
**Evidence pointer** `figure_raw_data/sim_headtohead/strainscan_single_persample.tsv` has four data rows, `strainscan_multi_persample.tsv` has three, and the cited `scratchpad/eval` drivers are absent.
**Concern** The primary comparator evidence is incomplete and the cited analysis code is unavailable.
**Why it matters** Tables 1 through 3, Figure 11, and abstract recall claims depend on these runs.
**Resolution test** Add complete per-sample outputs, cluster mappings, command lines, versions, timings, seeds, and executable analysis code. Demonstrate regeneration of every published median and range.

#### R1-M3 [statistical-rigor]

**Severity** Major
**Blocking** Yes
**Claim pointer** The manuscript claims improved recall, comparable multi-species accuracy, perfect host identification, and significantly lower abundance for 2bRAD-only calls.
**Evidence pointer** Tables 1 through 3 report point estimates without intervals; saliva includes 8 subjects for classification and only 3 usable paired samples for concordance; multiple species-level tests are reported without multiplicity adjustment.
**Concern** Comparative and real-data claims lack uncertainty, dependence-aware inference, and multiplicity control.
**Why it matters** Small, filtered, or repeated observations can produce unstable point estimates and inflated significance.
**Resolution test** Predefine primary paired comparisons, report bootstrap or permutation intervals, give classification confidence intervals, model within-sample dependence, and adjust or clearly label nominal p values.

#### R1-M4 [claim-moderation]

**Severity** Major
**Blocking** Yes
**Claim pointer** The manuscript states that Strain2bScan and StrainScan had F1 = 1.0 on clean MSA-1003, MSA-1005, and MSA-1007 shotgun samples.
**Evidence pointer** `data/fig6_fig12_metrics.tsv` reports lower raw F1 values for default WMS Strain2bScan on MSA-1005 and MSA-1007. Figure 12 displays F1 at 1e-4 rather than raw F1.
**Concern** Clean-mock accuracy claims conflict with the tracked metrics and rendered figure.
**Why it matters** The mock validation supports the shotgun-mode central claim.
**Resolution test** Define one primary threshold and one primary swept metric. Reconcile every text, figure, legend, and table number. Report clean-sample false positives explicitly.

#### R1-M5 [experimental-design]

**Severity** Major
**Blocking** Yes
**Claim pointer** The Abstract reports 130 to 146 times faster community profiling, and the text reports 8.5 times database and 6.5 times profiling parallel speedups.
**Evidence pointer** `results/community_throughput.tsv` labels StrainScan throughput as projected. `results/parallel_and_build_scaling.tsv` reports 4.6 and 5.8 times.
**Concern** Some runtime claims conflate measured and projected values and one thread-scaling claim conflicts with the source table.
**Why it matters** Scalability is a central claim.
**Resolution test** Rerun natively where possible or label every projected and emulated number consistently. Reconcile speedups with the tracked table and figure.

#### R1-M6 [statistical-rigor]

**Severity** Major
**Blocking** Yes
**Claim pointer** The manuscript reports 100 percent leave-one-timepoint-out host identification and a highly significant lower-abundance 2bRAD-only set.
**Evidence pointer** `results/saliva_permanova.tsv` reports 0.969 generic leave-one-out accuracy; `results/saliva_temporal_ml.tsv` reports 1.000 leave-one-timepoint-out accuracy; `results/saliva_concordance.tsv` has only three usable paired samples.
**Concern** Different validation protocols are conflated, intervals are absent, and the abundance test does not model dependence among calls.
**Why it matters** These are biological headline claims.
**Resolution test** Separate generic leave-one-out from leave-one-timepoint-out, report intervals, use sample-aware statistics, and avoid calling concordance validation.

#### R1-M7 [claim-moderation]

**Severity** Major
**Blocking** Yes
**Claim pointer** The Abstract presents 6/6 isolate self-calls and persistent strain tracking as cohort validation.
**Evidence pointer** Methods state that the same isolate reads generated the assemblies and were then used for self-recovery.
**Concern** Circular self-recovery cannot establish strain-calling accuracy.
**Why it matters** This affects the real-data accuracy narrative.
**Resolution test** Reword as panel and read compatibility plus hypothesis generation. Add independent isolates, spike-ins, or independently sequenced samples for accuracy validation.

#### R1-M8 [reproducibility]

**Severity** Major
**Blocking** Yes
**Claim pointer** The manuscript claims all scripts, pinned accessions, result tables, and figure code are in the repository and every figure is regenerable.
**Evidence pointer** The cited `results/real_data_strain_benchmark.md` and related result directories are untracked, `scratchpad/eval` is absent, and `make figures` does not regenerate every main figure. Some scripts contain machine-specific paths.
**Concern** Reproducibility is overstated.
**Why it matters** Auditable benchmark artifacts are central to a benchmark paper.
**Resolution test** Commit or archive all cited inputs, outputs, drivers, environment locks, and complete figure recipes at one immutable revision. Demonstrate a clean-room figure and analysis run.

### Minor Comments

- R1-m1 [writing-clarity]. Remove the duplicated public-cohort Methods section at lines 893 through 951.
- R1-m2 [writing-clarity]. Rename and clearly separate generic leave-one-out and leave-one-timepoint-out saliva metrics.
- R1-m3 [figures-and-tables]. Define raw precision, threshold precision, AUPR, and abundance metrics once and label every figure panel.
- R1-m4 [writing-clarity]. Replace broad `strain` terminology with `strain-resolved cluster` or `panel unit` where appropriate.
- R1-m5 [reproducibility]. Reconcile stale thread-scaling values in modular files with tracked tables and figures.
- R1-m6 [reproducibility]. Integrate numbered citations into the text and verify all bibliographic metadata.
- R1-m7 [data-resource-quality]. Provide complete Supplementary Tables S3 and S4, not only short legends.
- R1-m8 [ethical-governance]. Complete authors, affiliations, funding, ethics or data-use statements, accession mappings, and code archival identifiers.

### Assessment against Nature-style criteria

- Originality. Moderate to high. The native 2bRAD application is promising, although the algorithm is presented as a port of the StrainScan resolution framework.
- Scientific importance. Potentially high for low-biomass, high-host profiling and scalable multi-species screening.
- Interdisciplinary readership. Broad in principle, but currently dense and specialist.
- Technical soundness. Not established in the current version because of benchmark inconsistency and reproducibility gaps.
- Readability for nonspecialists. Mixed. The motivation is clear, but terminology and metric definitions need simplification.

### Recommendation posture

The current case is not established. Major revision with mandatory reanalysis, complete comparator artifacts, corrected claims, and independent reproducibility checks is required.

---

## Reviewer 2

### Overall assessment

The manuscript addresses an important gap by coupling a sparse 2bRAD marker space to a StrainScan-like clustering and unique-marker framework. The repository contains unusually detailed result tables, figure manifests, and analysis scripts. However, the central algorithmic and benchmarking case is not yet established because headline claims are not consistently supported by cited tables, the primary manuscript is stale relative to modular methods, the profiler source is not version-pinned, and key drivers are absent or machine-specific.

### Who would be interested in the results, and why

Microbiome researchers, metagenomic software developers, clinical microbiologists, and genome informaticians working on low-biomass, host-rich, or cohort-scale samples would be interested. Native 2bRAD-M strain profiling and digest-once multi-species profiling are both useful if reproducible.

### Major strengths

1. The biological motivation is clear.
2. The 15-species controlled benchmark is ambitious and includes single-species and multi-species communities.
3. Methods describe clustering, unique markers, species gates, adaptive singleton policy, consistency tests, abundance scopes, and cross-species restriction.
4. The figure package has checksums and a manifest.
5. The manuscript candidly discusses nonindependent cohort validation and closed-panel limitations.

### Major Concerns

#### R2-M1 [claim-moderation]

**Severity** Major
**Blocking** Yes
**Claim pointer** The Abstract and Results claim precision 1.0 at low DNA input, accurate abundance at 99 percent host, F1 = 1.0 on clean shotgun mocks, and a unique host-contamination advantage.
**Evidence pointer** `data/fig6_fig12_metrics.tsv` and Figures 6 and 12 show lower default-variant metrics at several operating points. The source table also contains `Strain2bScan-port`, `164-tracegap`, and `164-layers` variants not clearly defined in the primary manuscript.
**Concern** Text reports the most favourable interpretation across variants and metric families while figures display default or thresholded values.
**Why it matters** Readers cannot identify the primary algorithm configuration or determine whether the benchmark is selective.
**Resolution test** Predefine one primary variant and metric family. Define optional trace-gap and layer variants separately. Regenerate all abstract, figure, and table claims from the same specification.

#### R2-M2 [reproducibility]

**Severity** Major
**Blocking** Yes
**Claim pointer** The manuscript links to Rust source and makes detailed implementation and regression-test claims.
**Evidence pointer** The paper repository does not contain the profiler source, release tag, source commit, build hash, compiler version, or binary checksum.
**Concern** The implementation is not assessable at a reproducible version.
**Why it matters** This is primarily a software and algorithm paper.
**Resolution test** Archive the exact source for each benchmark in Zenodo or a tagged release. State source commit, release tag, Rust compiler, build command, binary SHA-256, and test provenance.

#### R2-M3 [writing-clarity]

**Severity** Major
**Blocking** Yes
**Claim pointer** The assembled manuscript is intended to be generated from modular sources, but contains contradictory and duplicated material.
**Evidence pointer** `manuscript/methods.md` documents `--marker-source kmer`; `full_manuscript.md` describes k-mer marker sources as future work. The full manuscript has duplicate public-cohort Methods sections at lines 893 through 951.
**Concern** The primary manuscript is not a fresh, consistent assembly of the latest modular sections.
**Why it matters** Readers cannot tell whether the k-mer mode is implemented, benchmarked, proposed, or unsupported.
**Resolution test** Regenerate the full manuscript from modular sources, add an assembly freshness check, remove duplicate sections, and state whether k-mer mode was benchmarked.

#### R2-M4 [reproducibility]

**Severity** Major
**Blocking** Yes
**Claim pointer** The manuscript claims complete scripts, pinned accessions, result tables, and figure code.
**Evidence pointer** Cited `scratchpad/eval` drivers are absent; several scripts hard-code private paths; some random selections use Python process-salted `hash`; a cited depth table appears to predate the corrected binary; and truth-cluster and profiling enzyme sets may differ.
**Concern** The benchmark pipeline is not demonstrably versioned or portable.
**Why it matters** Exact genomes, reads, cluster mappings, software, and scoring rules determine reproducibility.
**Resolution test** Restore missing drivers, parameterize paths, use stable seeds, lock environments, and reconcile truth clustering with the profiling marker space.

#### R2-M5 [statistical-rigor]

**Severity** Major
**Blocking** Yes
**Claim pointer** The manuscript reports 130 to 146 times faster 55-species profiling, 8.5 times database scaling, 6.5 times profiling scaling, and 204 same-environment timing samples.
**Evidence pointer** `results/community_throughput.tsv` labels StrainScan as projected; `results/parallel_and_build_scaling.tsv` reports 4.6 and 5.8 times; timing rows total 189 rather than 204.
**Concern** Runtime claims mix measured, projected, and different subsets without consistent labels.
**Why it matters** Scalability is one of the two central claims.
**Resolution test** Separate measured, projected, and emulated comparisons. Match accuracy and timing sets or report both denominators. Provide hardware, threads, repetitions, variability, and environment details.

#### R2-M6 [clinical-validity]

**Severity** Major
**Blocking** No
**Claim pointer** Real saliva and public WGS results are presented as validation.
**Evidence pointer** Saliva concordance has three usable paired samples. Isolate self-calls are nonindependent, and absent-isolate tests show conspecific misassignment.
**Concern** Real-data evidence is directional and small.
**Why it matters** Concordance and self-recovery do not establish independent accuracy.
**Resolution test** Report sample-aware confidence intervals, add independent validation, and consistently use compatibility or hypothesis-generation language for self-calls.

#### R2-M7 [reproducibility]

**Severity** Major
**Blocking** No
**Claim pointer** Public and reference-panel analyses use BioProject-level accessions and reference panels.
**Evidence pointer** Several panel scripts query dynamic NCBI listings, body-site manifests use local volume paths, and complete run-level or database checksums are not provided for every panel.
**Concern** Reference-panel composition and selected runs are not fully versioned.
**Why it matters** Panel composition is a central determinant of recall and abundance.
**Resolution test** Provide pinned accession lists, SRR run manifests, database commands, database hashes, binary commits, and environment locks for every panel.

### Minor Comments

- R2-m1. Define leave-one-out and leave-one-timepoint-out saliva protocols separately.
- R2-m2. Distinguish distinct panel markers from raw marker observations in public-cohort summaries.
- R2-m3. Rephrase the dependency claim as no third-party Rust crates and list external system commands such as `gzip`.
- R2-m4. Add a limitation that the 15-species simulation does not measure sequencing-error robustness.
- R2-m5. State whether species-level saliva p values are nominal or multiplicity adjusted.
- R2-m6. Remove the author-facing bibliography note and complete metadata.
- R2-m7. Simplify dense Figure 9 and Figure 12 layouts and label measured versus projected evidence.
- R2-m8. Reorder supplementary figure descriptions or explain the S3, S1, S2 ordering.

### Assessment against Nature-style criteria

- Originality. The native 2bRAD-M strain-profiling application is potentially novel and useful.
- Scientific importance. Potentially high if benchmark and reproducibility issues are fixed.
- Interdisciplinary readership. Broad in principle, although current figures and methods are dense.
- Technical soundness. Weakest area because of inconsistent metrics, versioning, and benchmark provenance.
- Readability for nonspecialists. Motivation is clear, but the primary manuscript is long and repetitive.

### Recommendation posture

Major revision is required. The study should report one prespecified configuration, pinned software, deterministic benchmarks, complete manifests, and measured runtime comparisons clearly separated from projections.

---

## Reviewer 3

### Overall assessment

The manuscript presents a potentially useful and well-motivated bridge from native BcgI 2bRAD-M libraries to scalable shotgun strain profiling. The biological motivation is strong and the manuscript is often explicit that resolution is cluster-level and panel-dependent.

However, the case is not established as written. Several headline mock statements do not match the displayed figures or metric table. Cohort evidence is small, selected, and partly circular. The saliva cross-modal comparison is directional rather than a validated reciprocal truth comparison. The manuscript is too long and dense for a broad readership, and several figure legends do not identify the operating points shown.

### Who would be interested in the results, and why

Microbiome researchers interested in low-biomass or host-rich specimens, oral-microbiome and clinical-metagenomics groups, bioinformatics method developers, and infant or body-site microbiome researchers would be interested if the claims are narrowed to the demonstrated settings.

### Major strengths

1. The paper addresses a real application gap.
2. The dual-input design is biologically coherent.
3. Reference completeness and panel design are treated as central biological issues.
4. The manuscript is transparent that short-read profiling usually resolves clusters, not isolate-level strains.
5. The evaluation spans simulations, mocks, saliva, public WGS subsets, cohort panels, and runtime benchmarks.

### Major Concerns

#### R3-M1 [claim-moderation]

**Severity** Major
**Blocking** Yes
**Claim pointer** The Results report favourable precision, recall, F1, and AUPR values for native low-biomass and clean shotgun mocks.
**Evidence pointer** `figures/numbered/Fig06_native_2brad_mocks.png`, `figures/numbered/Fig12_wms_toolcompare.png`, and `data/fig6_fig12_metrics.tsv` display lower values for the cited default rows.
**Concern** Prose reports more favourable operating points than the cited figures and source table.
**Why it matters** These are central accuracy claims.
**Resolution test** Reconcile every numeric claim with one selected source row and metric. State whether values are raw, thresholded at 1e-4, AUPR, or filtered.

#### R3-M2 [figures-and-tables]

**Severity** Major
**Blocking** Yes
**Claim pointer** Figures 6 and 12 are described as three-tool or default-method row-per-sample comparisons.
**Evidence pointer** The rendered figures display `AUPR`, `P@1e-4`, `R@1e-4`, `F1@1e-4`, and additional rows labelled `WMS (S2bS-port layers)` without adequate definition.
**Concern** The displayed variants and thresholds are not defined in the main narrative or legends.
**Why it matters** Different thresholds and variants materially change precision and F1.
**Resolution test** Label tool, database size, variant, and threshold for every row. Define all abbreviated metrics or remove secondary variants from the main figure.

#### R3-M3 [claim-moderation]

**Severity** Major
**Blocking** No
**Claim pointer** The manuscript invokes saliva, tumour, FFPE tissue, and skin as clinical niches.
**Evidence pointer** Native 2bRAD biological evidence is limited to ATCC mocks, eight saliva subjects, and four exploratory oral samples. No skin, tumour, or FFPE data are analysed.
**Concern** The demonstrated biological scope is broader than the tested materials.
**Why it matters** Host fraction, biomass, diversity, marker density, and degradation differ across these settings.
**Resolution test** Restrict demonstrated scope to mocks, saliva, and exploratory oral samples. Present tumour, FFPE, and skin only as motivation or future work.

#### R3-M4 [clinical-validity]

**Severity** Major
**Blocking** Yes for cohort-level biological interpretation, no for runtime feasibility
**Claim pointer** The Abstract describes public WGS cohorts, self-calls, and persistent strain tracking.
**Evidence pointer** The analyses include one pregnancy subject, one infant over three weeks, six isolate self-read sets, and no independent ground truth for the persistent strain.
**Concern** The evidence does not support cohort-level biological generalisation.
**Why it matters** The title and abstract imply more breadth than the selected subsets support.
**Resolution test** Label all real WGS analyses as exploratory demonstrations. Add SRR accessions, depth, body site, selection rationale, and isolate provenance. Present persistence as within-panel reidentification unless independently validated.

#### R3-M5 [claim-moderation]

**Severity** Major
**Blocking** No
**Claim pointer** The Abstract states that native BcgI 2bRAD-M confirmed 65/65 shotgun calls and recovered 128 to 163 additional strains per sample.
**Evidence pointer** The concordance set has three usable paired samples, shotgun R1 was a prefix subsample, and 2bRAD-only calls lack independent truth labels.
**Concern** Concordance is directional and does not validate additional native-only calls.
**Why it matters** Readers may interpret reference-panel calls as validated discoveries.
**Resolution test** State sample count and prefix-subsampling in the Abstract and Results. Replace `recovered` with `candidate strain-cluster calls` unless independent validation is provided.

#### R3-M6 [writing-clarity]

**Severity** Major
**Blocking** No
**Claim pointer** The manuscript is intended for broad biological readers while also serving as a technical benchmark paper.
**Evidence pointer** The front matter is incomplete, public-cohort Methods are duplicated, future-work text repeats, and Figure 9 has many panels not described by the caption.
**Concern** The manuscript is unfinished, repetitive, and too dense for nonspecialists.
**Why it matters** Important results are difficult to identify.
**Resolution test** Consolidate duplicates, complete front matter, shorten the main text, and make each figure panel self-identifying with concise takeaway captions.

### Minor Comments

- R3-m1. Distinguish 97 percent strain-level PCoA 1-NN accuracy from 100 percent leave-one-timepoint-out host ID and correct the cited panels.
- R3-m2. Reconcile 130 to 146 with 121 to 146 for community throughput.
- R3-m3. Add visible panel letters and captions for all plots in Figure 9.
- R3-m4. Include the actual Table S3 with sample metadata and availability status.
- R3-m5. Add SRR accessions, depth, body site, and selection criteria to Table 4.
- R3-m6. Complete authors, affiliations, funding, data availability, and reference verification.

### Assessment against Nature-style criteria

- Originality. Potentially strong for native 2bRAD-M strain profiling.
- Scientific importance. Potentially high after benchmark inconsistencies are corrected.
- Interdisciplinary readership. Broad problem, but current presentation is too specialist.
- Technical soundness. Mixed because numeric discrepancies prevent acceptance of the benchmark narrative.
- Readability for nonspecialists. Weak because the manuscript is overlong, repetitive, and figure-dense.

### Recommendation posture

Major revision. The central concept is promising, but the current benchmark, cohort, and readability issues prevent support in the present form.

---

## Cross-review synthesis

This synthesis was performed only after the three reviewer reports were frozen.

### Consensus strengths

1. The biological motivation is strong. Low-biomass, high-host specimens and community-scale shotgun cohorts are important and distinct use cases.
2. The dual-input design is coherent. Native BcgI 2bRAD-M and in-silico-digested shotgun reads enter one sparse marker framework.
3. The evaluation scope is broad, spanning simulations, reference degradation, mocks, saliva, public WGS subsets, isolate panels, and runtime scaling.
4. The manuscript is candid that strain output is cluster-level and that reference-panel design constrains interpretation.
5. The figure manifest and many result tables provide a useful starting point for reproducibility.

### Consensus blocking concerns

The following issues are shared by at least two reviewers and block the current central claims.

1. **Metric and variant inconsistency.** R1, R2, and R3 found that headline mock results do not consistently match the default rows shown in Figures 6 and 12 and the rows in `data/fig6_fig12_metrics.tsv`. The manuscript mixes raw metrics, 1e-4-threshold metrics, AUPR, trace-gap variants, and port-layer variants without a prespecified primary configuration.
2. **Reproducibility of the central benchmark.** R1 and R2 found missing or incomplete raw per-sample StrainScan outputs, absent `scratchpad/eval` drivers, machine-specific paths, unstable seeds, and incomplete figure regeneration.
3. **Software provenance.** R2 identified the lack of a pinned Strain2bScan release, source commit, compiler version, binary checksum, and per-benchmark binary identity. R1 raised the same broader reproducibility concern.
4. **Runtime claim provenance.** R1 and R2 found that community-throughput comparisons include projected StrainScan values and that the text conflicts with the tracked parallel-scaling table and figure.
5. **Statistical treatment of real-data claims.** R1 and R2 found small dependent samples, missing uncertainty, conflation of saliva validation protocols, and a p value that does not model within-sample dependence.
6. **Circular real-data validation language.** R1, R2, and R3 found that 65/65 saliva concordance and 6/6 isolate self-calls should not be described as independent validation.

### Other shared major concerns

- The primary assembled manuscript is stale relative to the modular methods and contains duplicated Methods text.
- Reference panels and public run selections require complete accession, run-level, database-hash, and software-version manifests.
- The biological scope exceeds the tested materials, especially for tumour, FFPE, and skin applications.
- Figures 6, 9, and 12 need clearer panel labels, operating thresholds, variant definitions, and caption-to-panel mapping.

### Where emphasis differs

Reviewer 1 focused most strongly on simulator mismatch, incomplete comparator artifacts, and statistical inference. Reviewer 2 focused most strongly on software release provenance, stale assembly, deterministic benchmarking, and separation of measured from projected runtime. Reviewer 3 focused most strongly on biological scope, cohort interpretation, figure readability, and claim moderation. These emphases are complementary rather than contradictory.

### Priority revision checklist

1. **Freeze one primary configuration.** State the exact Strain2bScan source commit, binary hash, compiler, command line, thresholds, reference-panel hash, and database hash for every benchmark.
2. **Predefine metrics.** Choose one primary detection threshold and one primary abundance-swept metric. Keep raw and thresholded metrics separate in every abstract, table, figure, and legend.
3. **Regenerate Figures 6 and 12.** Remove secondary port, trace-gap, or layer rows from the main figures or move them to a clearly labelled sensitivity figure. Explain every remaining row, variant, threshold, and metric.
4. **Reconcile every number.** Write a script that extracts headline numbers from the final result tables and fails if Abstract, Results, legends, or figures disagree.
5. **Restore benchmark artifacts.** Commit or archive complete per-sample StrainScan and Strain2bScan outputs, truth mappings, command lines, versions, seeds, environment locks, and the exact scoring scripts.
6. **Fix the simulator description.** State the actual read length, error model, insert size, seeds, and software versions. Audit or regenerate reads to match one specification.
7. **Correct runtime claims.** Label projected and measured values separately. Reconcile 121 to 146 versus 130 to 146 and 4.6/5.8 versus 8.5/6.5.
8. **Improve statistics.** Add paired intervals or tests for benchmark differences, classification intervals for saliva, sample-aware inference for call-level abundance comparisons, and multiplicity handling for species-level tests.
9. **Moderate real-data claims.** Use compatibility or concordance language for 65/65 and 6/6. State sample counts, prefix subsampling, and accession-level metadata.
10. **Reassemble the manuscript.** Remove duplicate Methods sections, reconcile k-mer mode, eliminate repeated future-work text, complete front matter, and regenerate from modular files with a freshness check.
11. **Rewrite figures for broad readers.** Add visible panel letters, define all abbreviations, state sample sizes, mark measured versus projected results, and move technical secondary panels to supplements.
12. **Complete submission materials.** Add authors, affiliations, funding, data-use statements, full Supplementary Tables S3 and S4, verified references, and code archival identifiers.

### Minor revision checklist

- Distinguish generic leave-one-out from leave-one-timepoint-out saliva results.
- Use `strain-resolved cluster` or `panel unit` rather than `strain` when output is merged.
- Reconcile distinct marker counts versus marker observations in cohort tables.
- Rephrase the dependency claim as no third-party Rust crates and list external commands.
- Reorder Supplementary Figures S1 through S4 or explain the current order.
- Remove the working-bibliography author note and verify all reference metadata.

### Broad-interest and significance readout

The problem is important and the dual-mode concept is potentially valuable to microbiome and clinical-metagenomics readers. The strongest possible contribution is not simply a faster StrainScan port. It is the demonstration that native 2bRAD-M data can support strain-resolved inference in high-host, low-biomass specimens while also scaling to community-level shotgun cohorts. That contribution is plausible, but the current evidence chain does not yet support it at a Nature-style standard.

### Most important issues before a strong case is established

1. One immutable software configuration and one primary benchmark configuration.
2. Complete, independently runnable benchmark artifacts for the StrainScan head-to-head.
3. Consistency between every headline number, figure row, source table, and legend.
4. Correct treatment of measured, projected, and emulated runtime claims.
5. Dependence-aware statistical reporting for saliva and cohort results.
6. Moderated language for concordance, self-recovery, and cohort-level interpretation.

---

## Risk and unsupported claims

### Unsupported novelty claims

The manuscript does not make an unqualified first-method claim, but broad readers may infer clinical or anatomical generality from tumour, FFPE, and skin examples. These should remain motivation or future work.

### Significance claims requiring moderation

- `confirmed` for 65/65 saliva calls should be `concordant` unless independent truth exists.
- `recovered` for 128 to 163 additional calls should be `candidate strain-cluster calls`.
- `recovered all six expected isolate self-calls` should always state that reads used for assembly were also used for recovery.
- `persistent strain tracking` should state that persistence is observed against a cohort-derived panel without independent strain truth.
- `cohort` should be qualified as exploratory application subsets unless a larger, pre-defined cohort analysis is added.

### Missing controls or validations

- Independent isolate read sets or spike-ins not used to build the panel.
- Native same-environment StrainScan runs for the 55-species throughput comparison.
- Complete comparator run artifacts for Tables 1 through 3.
- Clean-sample false-positive reporting for every tool and mock.
- Uncertainty and dependence-aware inference for saliva and cohort calls.

### Readability claims not assessable here

The final readability of a typeset, edited version cannot be assessed from Markdown alone. The current source is clearly too long and dense for a broad Nature-style readership.

### Reliance on partial material

The review did not perform external literature searches or rerun the full benchmarks. It assessed repository consistency, tracked evidence, code-to-text correspondence, and internal numerical provenance at commit `b3161a1`.
