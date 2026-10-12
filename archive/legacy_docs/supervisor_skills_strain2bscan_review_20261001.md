# Supervisor-Skills review: Strain2bScan

Date: 2026-10-01
Input: `manuscript/full_manuscript.md`, modular manuscript files, and result tables under `results/`.
Skills applied: `idea-evaluator`, `benchmark-paper-template`, `pre-submission-reviewer`, and
`paper-polish` principles.

## Executive verdict

**Accept with revisions.** The paper has a clear and valuable methodological niche: it moves strain-level
profiling into native low-biomass/high-host 2bRAD-M data while retaining a scalable shotgun mode. The
scientific significance is high, but not because the work is a paradigm shift. Its value is practical: it
turns an underused reduced-representation sequencing regime into a usable strain-resolution tool, and it
removes the per-species database-scanning bottleneck that limits community-scale strain profiling.

Before the current revision, the largest manuscript risk was quoting stale saliva metrics from the old
`support` abundance field. That has been corrected in the main manuscript and figure text. The next major
constraint is that the newest public-cohort analyses are application-oriented and lack exhaustive strain
truth; they should remain framed as validation/application subsets rather than a population-level
benchmark resource.

## Scientific significance and value

| Dimension | Assessment | Evidence |
|---|---|---|
| Higher | Enables strain-level profiling in native low-biomass/high-host 2bRAD-M specimens, a setting where shotgun strain tools lose signal. | ATCC 99% host results; real saliva concordance and low-abundance recovery. |
| Faster | Digest-once profiling makes multi-species cost approximately independent of species count. | 4–105× faster profiling than StrainScan; 130–146× faster on 55-species communities. |
| Stronger | Matches StrainScan precision, improves median recall on the matched benchmark, and preserves both detection and abundance at high host contamination. | Matched 15-species benchmark and ATCC host-contamination ladder. |
| Cheaper | Sparse markers reduce database construction and profiling cost. | 249–614× faster database construction; 43–138× lighter build; 343 MiB peak RSS across 22 real WGS libraries. |
| Broader | Supports native 2bRAD-M and in-silico-digested shotgun, spanning clinical low-biomass and cohort-scale settings. | Mocks, real saliva, public pregnancy/vaginal/meconium/preterm-infant WGS subsets. |

**Value statement.** The work is most valuable where existing strain profilers are weakest: high-host,
low-biomass clinical specimens and many-sample cohorts where full k-mer indexing or repeated per-species
querying becomes impractical. It also provides a useful evaluation axis—input modality, marker density,
reference completeness, host fraction, community size and cohort scale—that future strain profilers can
reuse.

**Benchmark-pillar assessment.**

1. **Research gap:** clear. The paper identifies both the scalability ceiling of full-k-mer strain tools
   and the inability of host-dominated shotgun to recover low-abundance strains.
2. **Construction pipeline:** strong for controlled panels, simulations and mocks. The new isolate-derived
   cohort panel adds a useful real-data path, but panel assembly and filtering criteria should eventually
   be standardized.
3. **Evaluation framework:** strong. It spans accuracy, depth sensitivity, reference completeness, enzyme
   density, host contamination, abundance accuracy, runtime and memory.
4. **Empirical findings:** strong. Results separate detection, abundance, scalability and reference-quality
   effects rather than reporting one aggregate score.
5. **Companion method:** Strain2bScan itself is the companion implementation, which is appropriate.

## Quality issues found and fixed

| Severity | Issue | Action |
|---|---|---|
| CRITICAL | The saliva section still quoted stale metrics from the old `support` field. | Replaced with corrected `sample_fraction` results: PERMANOVA R² 0.757 strain versus 0.755 species; leave-one-timepoint-out host ID 100% versus 78.1%. |
| MAJOR | Abstract did not state the new public-cohort results. | Added the 22-library runtime, 343 MiB memory, 6/6 isolate self-calls and persistent *B. bifidum* LHCA82 tracking. |
| MAJOR | Data availability omitted the four new BioProjects. | Added PRJNA288562, PRJNA1517970, PRJNA1191223 and PRJNA1191225. |
| MAJOR | The benchmark contribution was implicit. | Added an explicit benchmark-scope paragraph separating ground-truth, application and scalability claims. |
| MAJOR | Figure 7 text used stale saliva metrics. | Updated results, figures and detailed legend to corrected values. |
| MINOR | Isolate self-validation could be read as independent validation. | Clarified that it is a panel/read compatibility check because reads generated the assemblies. |
| MINOR | Abstract used broad “recovered all isolates” language. | Changed to “expected self-call for all six study isolates”. |

## Remaining work

1. Generate a body-site-specific panel for PRJNA1517970. The current no-call result is panel-limited and
   should not be interpreted as absence of microbial strains.
2. Perform leave-one-isolate-out validation so isolate classification is non-circular.
3. Deposit the six isolate assemblies and add accession numbers.
4. Reduce reliance on em dashes throughout the manuscript. The full draft still contains many instances;
   replace sentence-connecting em dashes with commas, colons, parentheses or separate sentences.
5. Run a final consistency pass over `thesis_chapter.md`, which still contains the stale saliva metrics.
6. Add external validation on a cohort with isolate or long-read ground truth.

## Recommendation

**Needs minor-to-moderate revision before submission.** The method and evidence base are strong, but the
public-cohort analyses should remain clearly framed as application subsets, and the remaining figure and
manuscript text should be made fully consistent with corrected abundance fields.

## Revision implemented on 2026-10-02

1. Regenerated Figure 7 and companion metrics from corrected `sample_fraction` tables.
2. Built and tested a 13-species, 146-genome body-site panel for PRJNA1517970. All seven selected
   libraries remained negative; a one-marker diagnostic test in the largest meconium library found only
   three *Cutibacterium acnes* markers.
3. Added leave-one-isolate-out *E. coli* tests. When the true isolate was absent, reads were assigned to
   a merged cluster containing the two closest relatives rather than inventing a false third cluster.
4. Prepared a submission-ready package for six isolate assemblies with manifests and SHA-256 checksums.
5. Removed all sentence-connecting em dashes from the main manuscript and cleaned resulting punctuation.
6. Updated Results, Methods, Discussion, Tables 4 and 7, Figure 7 text and the thesis chapter to the
   corrected saliva metrics and new real-cohort analyses.

Remaining non-manuscript action: obtain NCBI accession numbers after submitting the prepared isolate
assembly package.
