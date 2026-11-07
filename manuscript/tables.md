# Strain2bScan; Manuscript tables (simulated head-to-head)

**Table 1. Strain2bScan versus StrainScan on the 15-species simulated benchmark.** Accuracy is the median over completed run-level paired samples (different-cluster mixtures for 14 species; same-cluster mixtures for near-clonal *M. tuberculosis*; k = 2/3/5; depths 0.5–10×); each tool was scored in its own 0.95-similarity cluster space. Cluster-count columns are each tool's own 0.95 partitions of the common species pool. StrainGST F1 uses all 15 species-specific tests in its own 0.90-reference space and is not pooled with the 0.95-space columns. Strain2bScan values are from source commit `6ae760c2f3bd2f39efc0e6aa06c53a84366cb16e`; source table: `results/strainscan_rerun/current_single_persample.tsv`. Build cost is the cross-platform timing benchmark (Strain2bScan native arm64; StrainScan `linux/amd64` under emulation). *S. enterica* produced no completed StrainScan calls; *K. pneumoniae* was successfully built and profiled in the final cross-platform analysis.

| Species | Paired n | Pool genomes | S2B clusters | SS clusters | S2B P | S2B R | S2B F1 | SS P | SS R | SS F1 | Separate StrainGST F1 (0.90-space)† | S2B build | SS build |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| *A. muciniphila* | 14 | 50 | 24 | 24 | 0.900 | 0.900 | 0.900 | 1.000 | 1.000 | 1.000 | 1.000 | 2.8 s / 0.16 GB | 17 min / 16 GB |
| *C. difficile* | 15 | 47 | 27 | 24 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.000 | 3.2 s / 0.28 GB | 30 min / 14 GB |
| *C. acnes* | 15 | 43 | 8 | 8 | 1.000 | 0.800 | 0.800 | 1.000 | 0.800 | 0.889 | 1.000 | 1.5 s / 0.17 GB | 11 min / 8 GB |
| *E. coli* | 15 | 50 | 30 | 27 | 1.000 | 1.000 | 1.000 | 1.000 | 0.800 | 0.889 | 1.000 | 4.3 s / 0.34 GB | 43 min / 28 GB |
| *F. nucleatum* | 15 | 25 | 19 | 19 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.000 | 0.7 s / 0.10 GB | 5 min / 8 GB |
| *K. pneumoniae* | 15 | 47 | 28 | 27 | 0.333 | 0.333 | 0.286 | 1.000 | 0.800 | 0.889 | 1.000 | 5.1 s / 0.40 GB | DNF (timing run) |
| *L. plantarum* | 15 | 50 | 39 | 39 | 1.000 | 1.000 | 1.000 | 1.000 | 0.800 | 0.889 | 1.000 | 5.0 s / 0.24 GB | 24 min / 18 GB |
| *M. tuberculosis* | 15 | 29 | 1 | 1 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.6 s / 0.25 GB | 16 min / 15 GB |
| *P. dorei* | 10 | 15 | 10 | 10 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.000 | 1.0 s / 0.25 GB | 5 min / 13 GB |
| *P. gingivalis* | 15 | 43 | 27 | 31 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.000 | 2.3 s / 0.15 GB | 19 min / 11 GB |
| *P. copri* | 15 | 19 | 19 | 19 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.000 | 1.8 s / 0.19 GB | 8 min / 26 GB |
| *S. enterica* | 0 | 45 | 19 | 12 | -- | -- | -- | -- | -- | -- | 1.000 | 4.2 s / 0.37 GB | 34 min / 16 GB |
| *S. aureus* | 15 | 48 | 13 | 12 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.000 | 1.7 s / 0.17 GB | 17 min / 12 GB |
| *S. epidermidis* | 15 | 50 | 9 | 9 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 0.800 | 1.6 s / 0.14 GB | 15 min / 8 GB |
| *S. pneumoniae* | 15 | 50 | 38 | 38 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.000 | 1.5 s / 0.16 GB | 11 min / 11 GB |

**Median across 14 species with completed paired runs:** Strain2bScan P 1.000 / R 1.000 / F1 1.000; StrainScan P 1.000 / R 0.667 / F1 0.800. Species medians conceal a precision/recall trade-off for *K. pneumoniae*. Mean paired differences and 95% intervals are in `results/strainscan_rerun/current_uncertainty_summary.tsv`.

**Table 2. Accuracy and same-container per-sample cost by sequencing depth (single-species medians).** Accuracy uses final Strain2bScan and StrainScan on the 204 completed paired runs from `results/strainscan_rerun/current_single_by_depth.tsv`; StrainGST accuracy uses its separate 0.90-reference space from `results/strainge_rerun/simulation_single_by_depth.tsv` and is not pooled. Cost is the median of the same-container subset in `results/sim_headtohead_profile_sameenv.tsv`; both tools ran as linux/amd64 containers. Cost n is reported separately from accuracy n.

| Depth (×) | Accuracy n | S2B P/R/F1 | StrainScan P/R/F1 | StrainGST P/R/F1† | Cost n | S2B time / RSS | StrainScan time / RSS |
|--:|--:|:--:|:--:|:--:|--:|:--:|:--:|
| 0.5 | 40 | 1.000/0.600/0.708 | 1.000/0.667/0.800 | 1.000/1.000/1.000 | 37 | 0.160 s / 21.1 MB | 5.300 s / 830.9 MB |
| 1 | 41 | 1.000/1.000/1.000 | 1.000/1.000/1.000 | 1.000/1.000/1.000 | 38 | 0.240 s / 29.9 MB | 3.475 s / 830.8 MB |
| 3 | 41 | 1.000/1.000/1.000 | 1.000/0.600/0.750 | 1.000/1.000/1.000 | 38 | 0.560 s / 60.0 MB | 4.165 s / 830.9 MB |
| 5 | 41 | 1.000/1.000/1.000 | 1.000/0.667/0.800 | 1.000/1.000/1.000 | 38 | 0.890 s / 90.8 MB | 4.790 s / 830.9 MB |
| 10 | 41 | 1.000/1.000/1.000 | 1.000/1.000/1.000 | 1.000/1.000/1.000 | 38 | 1.655 s / 160.7 MB | 6.915 s / 831.8 MB |

Strain2bScan reached median recall 1.0 from 1× onward; StrainScan was non-monotonic and did not exceed 0.800 median F1 at 0.5× or intermediate depths. The depth table aggregates heterogeneous strain mixtures and is not a causal depth-onset estimate.

**Table 3. Multi-species community accuracy and profiling cost.** Accuracy is the mean of four samples per depth. Strain2bScan and StrainScan use 0.95-cluster spaces from `results/strainscan_rerun/current_multi_by_depth.tsv` and `results/strainscan_rerun/multi_by_depth.tsv`; StrainGST uses its separate 0.90-reference space from `results/strainge_rerun/simulation_multi_by_depth.tsv`. Costs are medians from the timing tables (`results/sim_headtohead_multi.tsv` and `results/strainge_rerun/simulation_run_efficiency.tsv`) and are not pooled with accuracy. Strain2bScan profiles each community in one digest-once pass; StrainScan and StrainGST costs are medians of per-sample totals across species searches.

| Community depth | S2B P/R/F1 | StrainScan P/R/F1 | StrainGST P/R/F1† | S2B time / RSS | StrainScan time / RSS | StrainGST time / RSS |
|---|:--:|:--:|:--:|:--:|:--:|:--:|
| low | 0.748/0.667/0.699 | 0.938/0.702/0.798 | 0.918/1.000/0.944 | 0.96 s / 311 MB | 100.4 s / 1112 MB | 163.8 s / 2.923 GB |
| med | 0.838/0.838/0.837 | 0.954/0.903/0.927 | 0.901/1.000/0.934 | 4.33 s / 670 MB | 228.2 s / 1696 MB | 237.0 s / 7.896 GB |
| high | 0.806/0.806/0.806 | 0.972/0.958/0.965 | 0.908/1.000/0.951 | 8.70 s / 1119 MB | 397.9 s / 2028 MB | 350.0 s / 10.705 GB |

Of 180 species-by-community opportunities (15 databases × 12 samples), 158 StrainScan final reports were generated; absent reports were treated as no detection. The per-sample coverage audit is in `results/strainscan_rerun/multi_persample.tsv`.

**Table 4. Public real-metagenome application subsets profiled with Strain2bScan.** These were
informative application subsets rather than complete cohort analyses. The generic panel was the legacy
20-species MSA database; calls used all-enzyme mode with `--min-species-markers 20
--min-species-detect 2 --min-support 2 --min-coverage 0.01 --min-abundance 0`.

| Project | Design and selected libraries | Calls | Observation |
|---|---|--:|---|
| PRJNA288562 | Pregnancy subject T23; saliva, vaginal swab and distal gut at GD84 and GD273 (6 WGS libraries) | 34 calls in 5/6 libraries | Saliva was cluster-rich; four *Neisseria* and four *Schaalia* clusters were shared across timepoints. Gut *B. adolescentis* C4 persisted, whereas three other *Bifidobacterium* clusters and three *E. coli* clusters were GD84-only. |
| PRJNA1517970 | Preterm-birth vaginal/meconium subset plus extraction blank (7 WGS libraries) | 0 calls with the MSA panel; 0 calls with a 13-species body-site panel | The blank had 340 distinct markers, and selected libraries had 7,394–86,980. At a one-marker diagnostic gate, one meconium library showed only three *C. acnes* markers, supporting specificity while also indicating incomplete niche coverage. |
| PRJNA1191223 | Preterm infant P08 stool at W1, W2 and W3 (3 WGS libraries) | 9 generic-panel calls | Generic panel showed weekly turnover: one *S. aureus* cluster (W1), two *E. coli* clusters (W2), and six *E. faecalis* clusters (W3). |
| PRJNA1191225 | Six preterm-infant isolate WGS read sets used for cohort-panel validation | 6/6 expected self-calls | The cohort-specific panel recovered all expected *E. coli* clusters or Bifidobacterium species units. |

**Table 5. Isolate assemblies and cohort-specific panel validation.** Isolates were assembled with
SPAdes `--isolate`. *E. coli* assemblies were clustered at 0.95 similarity; each Bifidobacterium
species was built as a single-genome database. All isolate read sets recovered the expected panel unit.

| Isolate | Species | Contigs | Assembly (Mb) | N50 (kb) | Panel unit | Panel markers | Self-call |
|---|---|--:|--:|--:|---|--:|:--:|
| LHCA45 | *Escherichia coli* | 185 | 4.96 | 240.0 | C0 | 7,786 unique | ✓ |
| LHCA56 | *Escherichia coli* | 448 | 5.15 | 222.6 | C1 | 5,466 unique | ✓ |
| LHCA72 | *Escherichia coli* | 457 | 5.04 | 159.2 | C2 | 5,665 unique | ✓ |
| LHCA43 | *Bifidobacterium longum* | 164 | 2.36 | 55.6 | LHCA43 | 18,426 | ✓ |
| LHCA81 | *Bifidobacterium breve* | 90 | 2.39 | 228.3 | LHCA81 | 19,214 | ✓ |
| LHCA82 | *Bifidobacterium bifidum* | 180 | 2.30 | 75.6 | LHCA82 | 17,985 | ✓ |

**Table 6. P08 weekly stool profiling with the cohort-specific isolate panel.** Breadth is database
coverage; abundance is within-species abundance in Strain2bScan output.

| Week | Species | Panel unit | Breadth | Depth (×) | Abundance | Fraction of sample |
|--:|---|---|--:|--:|--:|--:|
| W1 | *Bifidobacterium bifidum* | LHCA82 | 0.764 | 236.7 | 1.000 | 0.6456 |
| W2 | *Bifidobacterium bifidum* | LHCA82 | 0.773 | 213.9 | 1.000 | 0.5834 |
| W2 | *Escherichia coli* | C0\|C1\|C2 | 0.026 | 0.029 | 1.000 | 0.000140 |
| W3 | *Bifidobacterium bifidum* | LHCA82 | 0.792 | 238.2 | 1.000 | 0.6099 |

The generic MSA panel detected the W2 *E. coli* signal as two low-abundance generic clusters
(fraction 6.1 × 10⁻⁵ and 5.7 × 10⁻⁵) but did not contain the persistent *B. bifidum* LHCA82 unit. This
contrast illustrates why body-site- or cohort-specific panels are needed for biological interpretation.


**Table 7. Leave-one-isolate-out behaviour for the three *Escherichia coli* panel isolates.** For each test, the
held-out isolate was absent from the database. Reads were profiled against a panel containing the other
two *E. coli* assemblies and the three Bifidobacterium species controls.

| Held-out isolate | *E. coli* genomes in panel | Observed *E. coli* unit | Breadth | Depth (×) | Fraction of sample | Interpretation |
|---|--:|---|--:|--:|--:|---|
| LHCA45 | 2 | C0\|C1 | 0.799 | 7.87 | 0.612 | Reads from the absent isolate were assigned to merged relatives rather than creating a false third cluster. |
| LHCA56 | 2 | C0\|C1 | 0.911 | 10.55 | 0.649 | Same conspecific-assignment behaviour. |
| LHCA72 | 2 | C0\|C1 | 0.891 | 7.25 | 0.729 | Same conspecific-assignment behaviour. |

This test shows that closed panels can misattribute a truly absent conspecific strain to near relatives.
It therefore supports the paper's conservative treatment of the low-coverage W2 *E. coli* signal as an
unresolved `C0|C1|C2` unit rather than assigning it to one isolate.

**Table 8. Comparator and assay scope.** This table distinguishes tools executed in the benchmark from
related tools discussed for context. It is a scope statement, not a performance comparison. Native 2bRAD
denotes BcgI-derived experimental libraries; shotgun denotes conventional WMS reads, which Strain2bScan
can also process after in-silico digestion.

| Tool | Executed here | Input scope | Typical resolution | Scope and role |
|---|:--:|---|---|---|
| Strain2bScan | Yes | Native BcgI 2bRAD reads and in-silico-digested shotgun reads | Strain clusters | Method under test; one marker framework supports both reduced-representation and shotgun input. |
| StrainScan v1.0.14 | Yes | Shotgun reads against per-species reference databases | Strain clusters | Direct reference-guided comparator on shotgun reads (Fig 7, Fig 9, Supplementary Fig S5; Tables 1–3). |
| inStrain v1.9.1 (mocks) and v1.11.0 (simulations) | Yes | Shotgun reads aligned to dereplicated mock references or full species simulation pools | Strain, SNV and genome coverage | Direct read-level shotgun comparator (Fig 7, Fig 9, Fig S4; simulation outputs in `results/instrain_sim225/`). |
| StrainGE / StrainGST / StrainGR | StrainGST only | Shotgun WMS reads against StrainGE 0.90-reference databases; native BcgI 2bRAD is not documented input | Strain level, reference guided | StrainGST completed 225 single-species simulations, 12 multi-species communities and four ATCC WMS mocks (Tables 1–3 and 9); StrainGR was not run. |
| 2bRAD-M / Fast2bRAD-M | No | Native BcgI 2bRAD reads | Mainly species | Same restriction-assay family; provides the species-level upstream context for Strain2bScan. |
| sylph | No | Shotgun sketches | Mainly species and ANI | Context for rapid species-level profiling; not a strain-resolved benchmark here. |
| StrainPhlAn | No | Shotgun marker genes | Strain types | Context for marker-gene strain profiling; not evaluated here. |

**Table 9. StrainGST (StrainGE 1.3.9) accuracy and efficiency.** Simulations are the matched
`_rep1_` subset and all 12 multi-species communities. Mocks are the primary WMS replicate used in
Figure 9; MSA-1002 is the 99%-host sample. Accuracy is in StrainGE's 0.90-reference space; mock calls
use abundance >= 1e-4. Wall time is one sample k-merization plus all per-species StrainGST searches; RSS
is the maximum across those stages. StrainGR was not run.

| Dataset | n | Precision | Recall | F1 | Median wall time/sample (s) | Median peak RSS/sample (GB) |
|---|--:|--:|--:|--:|--:|--:|
| Simulated single-species | 225 | 1.000 | 1.000 | 1.000 | 18.66 | 1.30 |
| Simulated multi-species | 12 | 0.901 | 1.000 | 0.941 | 236.98 | 7.90 |
| MSA-1002, 99% host | 1 | 0.700 | 0.700 | 0.700 | 1007.22 | 14.36 |
| MSA-1003 | 1 | 0.464 | 0.650 | 0.542 | 2158.75 | 21.52 |
| MSA-1005 | 1 | 0.333 | 0.833 | 0.476 | 1784.50 | 20.10 |
| MSA-1007 | 1 | 0.429 | 1.000 | 0.600 | 1555.76 | 17.69 |

**Table S1. Exploratory native BcgI 2bRAD-M oral-library profiling.** Four libraries from the
PRJNA1131785 exploratory oral cohort were profiled against the same 19-species panel used for the saliva
analysis. No case/control labels were available in public metadata. Values are from
`results/clinical_exploratory.tsv`.

| Sample | Distinct markers | Species resolved | Strain calls | Runtime (s) | Top species |
|---|--:|--:|--:|--:|---|
| S_0313206 | 271,322 | 15 | 115 | 2.3 | *Neisseria subflava*; *Haemophilus parainfluenzae*; *Actinomyces odontolyticus* |
| S_1413231 | 288,019 | 17 | 129 | 1.8 | *Neisseria subflava*; *Streptococcus mitis*; *Actinomyces odontolyticus* |
| S_8123213 | 370,379 | 16 | 149 | 2.0 | *Actinomyces odontolyticus*; *Neisseria subflava*; *Haemophilus parainfluenzae* |
| S_8221232 | 488,180 | 17 | 158 | 1.8 | *Actinomyces odontolyticus*; *Neisseria subflava*; *Rothia mucilaginosa* |

**Table S2. Failure-aware sensitivity for all 225 matched single-species entries.** Absent StrainScan final reports were retained and scored as zero detection. Failure-as-absence assigns zero recall to absent reports and uses a species-cluster bootstrap; precision and F1 are not assigned to no-report failures. Source: `results/strainscan_rerun/current_single_persample.tsv` and `results/strainscan_rerun/current_uncertainty_summary.tsv`.

| Analysis | n | No-calls | Mean recall difference (95% CI) | No-call rate |
|---|--:|--:|:--:|--:|
| Completed pairs | 204 | 0 | +0.063 (−0.037 to 0.133) | 0% |
| Failure as absence | 225 | 21 | +0.127 (0.007 to 0.248) | 9.33% |


**Table S3. Public-run manifest and exploratory profiling outcomes.** SRA spots, bases and instrument metadata were retrieved from NCBI on 2026-10-11. Marker and call counts are from the generic MSA screen; body-site-panel calls are shown for PRJNA1517970. P08 run-to-week mapping follows SRA submitter titles. Source: `results/realworld_cohort_benchmark/public_run_manifest.tsv`.


| Run | Project | Role/site | SRA spots | SRA bases | MSA markers | MSA strain calls | Body-site calls |
|---|---|---|--:|--:|--:|--:|--:|
| SRR32121885 | PRJNA1191223 | preterm infant P08 stool W1 | 4,862,900 | 1,215,725,000 | 748,209 | 1 | -- |
| SRR32121884 | PRJNA1191223 | preterm infant P08 stool W2 | 4,684,342 | 1,171,085,500 | 694,111 | 2 | -- |
| SRR32121883 | PRJNA1191223 | preterm infant P08 stool W3 | 4,403,456 | 1,100,864,000 | 830,986 | 6 | -- |
| SRR32076429 | PRJNA1191225 | Bifidobacterium bifidum LHCA82 isolate WGS | 265,039 | 78,810,882 | 129,681 | 0 | -- |
| SRR32076422 | PRJNA1191225 | Bifidobacterium breve LHCA81 isolate WGS | 292,272 | 87,006,759 | 142,707 | 0 | -- |
| SRR32076444 | PRJNA1191225 | Bifidobacterium longum LHCA43 isolate WGS | 117,024 | 34,917,748 | 65,183 | 0 | -- |
| SRR32076385 | PRJNA1191225 | Escherichia coli LHCA45 isolate WGS | 530,692 | 157,695,345 | 219,245 | 0 | -- |
| SRR32076393 | PRJNA1191225 | Escherichia coli LHCA56 isolate WGS | 694,301 | 206,230,492 | 242,121 | 0 | -- |
| SRR32076407 | PRJNA1191225 | Escherichia coli LHCA72 isolate WGS | 405,951 | 120,852,478 | 183,313 | 0 | -- |
| SRR40410776 | PRJNA1517970 | extraction blank | 969 | 181,397 | 464 | 0 | 0 |
| SRR40410528 | PRJNA1517970 | meconium | 130,092 | 35,779,278 | 116,410 | 0 | 0 |
| SRR40410657 | PRJNA1517970 | meconium | 63,704 | 12,967,250 | 31,171 | 0 | 0 |
| SRR40410851 | PRJNA1517970 | meconium | 37,274 | 8,265,014 | 23,718 | 0 | 0 |
| SRR40410567 | PRJNA1517970 | vagina | 29,852 | 5,373,843 | 9,693 | 0 | 0 |
| SRR40410734 | PRJNA1517970 | vagina | 123,942 | 24,263,782 | 24,067 | 0 | 0 |
| SRR40410753 | PRJNA1517970 | vagina | 71,724 | 13,985,625 | 14,802 | 0 | 0 |
| SRR6748033 | PRJNA288562 | subject T23 distal gut GD273 | 3,803,560 | 1,120,072,205 | 1,595,769 | 2 | -- |
| SRR6748028 | PRJNA288562 | subject T23 distal gut GD84 | 2,804,411 | 824,219,750 | 1,918,932 | 7 | -- |
| SRR6748088 | PRJNA288562 | subject T23 saliva GD273 | 771,556 | 226,981,334 | 876,200 | 11 | -- |
| SRR6748097 | PRJNA288562 | subject T23 saliva GD84 | 377,872 | 111,169,454 | 533,573 | 13 | -- |
| SRR6747958 | PRJNA288562 | subject T23 vagina GD273 | 636,817 | 186,911,987 | 226,983 | 1 | -- |
| SRR6748207 | PRJNA288562 | subject T23 vagina GD84 | 243,062 | 71,385,301 | 45,805 | 0 | -- |

**Table S4. Gate-matched full-depth saliva comparison across all 32 PRJNA1131785 pairs.** WMS and native BcgI 2bRAD libraries were profiled with the same 19-species oral panel, species gate (20 panel-specific markers and two detected markers), support threshold (two markers) and coverage threshold (0.01). Cluster calls were matched by species and cluster label. Differences in sequencing depth, read layout and representation prevent interpretation as accuracy validation. Source: `results/wms32_rebuild/per_sample_summary.tsv`.

| Sample | WMS calls | Native calls | Shared | WMS-only | Native-only | WMS species |
|---|--:|--:|--:|--:|--:|--:|
| S11-1 | 89 | 36 | 24 | 65 | 12 | 16 |
| S11-11 | 77 | 52 | 22 | 55 | 30 | 19 |
| S11-12 | 83 | 71 | 26 | 57 | 45 | 15 |
| S11-14 | 136 | 25 | 10 | 126 | 15 | 18 |
| S11-2 | 73 | 56 | 12 | 61 | 44 | 16 |
| S11-6 | 72 | 35 | 14 | 58 | 21 | 16 |
| S11-7 | 104 | 18 | 7 | 97 | 11 | 18 |
| S11-9 | 39 | 23 | 13 | 26 | 10 | 11 |
| S13-1 | 82 | 53 | 17 | 65 | 36 | 17 |
| S13-11 | 85 | 52 | 10 | 75 | 42 | 14 |
| S13-12 | 49 | 61 | 17 | 32 | 44 | 13 |
| S13-14 | 105 | 22 | 5 | 100 | 17 | 17 |
| S13-2 | 84 | 56 | 16 | 68 | 40 | 14 |
| S13-6 | 82 | 30 | 14 | 68 | 16 | 17 |
| S13-7 | 114 | 33 | 17 | 97 | 16 | 18 |
| S13-9 | 78 | 20 | 9 | 69 | 11 | 17 |
| S17-1 | 83 | 35 | 19 | 64 | 16 | 17 |
| S17-11 | 71 | 67 | 15 | 56 | 52 | 13 |
| S17-12 | 63 | 37 | 21 | 42 | 16 | 16 |
| S17-14 | 121 | 25 | 7 | 114 | 18 | 18 |
| S17-2 | 56 | 67 | 18 | 38 | 49 | 15 |
| S17-6 | 63 | 25 | 12 | 51 | 13 | 17 |
| S17-7 | 80 | 27 | 18 | 62 | 9 | 15 |
| S17-9 | 47 | 25 | 16 | 31 | 9 | 15 |
| S9-1 | 81 | 59 | 22 | 59 | 37 | 18 |
| S9-11 | 70 | 71 | 20 | 50 | 51 | 15 |
| S9-12 | 45 | 103 | 27 | 18 | 76 | 10 |
| S9-14 | 109 | 23 | 6 | 103 | 17 | 17 |
| S9-2 | 55 | 65 | 16 | 39 | 49 | 15 |
| S9-6 | 38 | 27 | 20 | 18 | 7 | 12 |
| S9-7 | 116 | 52 | 15 | 101 | 37 | 17 |
| S9-9 | 37 | 17 | 10 | 27 | 7 | 11 |

**Table S5. Rhesus macaque 148-sample cohort and monkey MAG panel.** The cohort contains 128 GI spatial
samples (66 luminal and 62 mucosal) plus 20 body controls. The monkey-specific panel used 315 selected
study MAGs across 18 species and 303 0.95 strain-cluster units. Counts use the same BcgI gates listed in
Methods. Source: `results/monkey_spatial_148/cohort_panel_summary.tsv`.

| Metric | Value |
|---|--:|
| Full cohort samples | 148 |
| GI spatial samples | 128 |
| GI luminal samples | 66 |
| GI mucosal samples | 62 |
| Body-control samples | 20 |
| Excluded vendor extra samples | 4 |
| Panel species | 18 |
| Selected study MAGs | 315 |
| Strain-cluster units | 303 |
| Panel markers before cross-species filtering | 30,008 |
| Usable markers after cross-species filtering | 29,162 |
| GI samples with ≥1 2bRAD strain call | 113/128 (88.3%) |
| Luminal samples with ≥1 call | 56/66 (84.8%) |
| Mucosal samples with ≥1 call | 57/62 (91.9%) |
| Body controls with ≥1 call | 2/20 (10.0%) |
| Total body-control calls | 3 |

**Table S6. Species prevalence in the rhesus macaque GI cohort.** Prevalence is the fraction of samples in
the indicated group with at least one strain call for the species. Source:
`results/monkey_spatial_148/species_GI_control_prevalence.tsv`.

| Species | GI | Luminal | Mucosal | Body control |
|---|--:|--:|--:|--:|
| *Helicobacter_B_sp003364135* | 0.766 | 0.727 | 0.806 | 0.100 |
| *UBA4372_sp900543815* | 0.664 | 0.727 | 0.597 | 0.000 |
| *Cryptobacteroides_sp000433355* | 0.656 | 0.773 | 0.532 | 0.050 |
| *Succinivibrio_sp000431835* | 0.578 | 0.652 | 0.500 | 0.000 |
| *WG_1_sp900539665* | 0.555 | 0.667 | 0.435 | 0.000 |
| *Oribacterium_sp900768305* | 0.500 | 0.530 | 0.468 | 0.000 |
| *CAG_873_sp001701165* | 0.492 | 0.591 | 0.387 | 0.000 |
| *Lachnospira_eligens_A* | 0.477 | 0.530 | 0.419 | 0.000 |
| *Treponema_D_succinifaciens* | 0.445 | 0.545 | 0.339 | 0.000 |
| *Treponema_D_berlinense* | 0.445 | 0.545 | 0.339 | 0.000 |
| *Ligilactobacillus_animalis* | 0.422 | 0.545 | 0.290 | 0.000 |
| *Cryptobacteroides_sp900546445* | 0.422 | 0.485 | 0.355 | 0.000 |
| *Ornithospirochaeta_sp934247875* | 0.414 | 0.500 | 0.323 | 0.000 |
| *Onthomorpha_sp900545215* | 0.414 | 0.470 | 0.355 | 0.000 |
| *UBA636_sp900546285* | 0.383 | 0.455 | 0.306 | 0.000 |
| *UBA1724_sp900548225* | 0.367 | 0.424 | 0.306 | 0.000 |
| *CAG_349_sp003539515* | 0.344 | 0.424 | 0.258 | 0.000 |
| *Ruminococcus_E_sp934476515* | 0.312 | 0.379 | 0.242 | 0.000 |

**Table S7. Representative LI-enriched strain clusters.** LI enrichment required peak LI prevalence ≥ 0.3
and region specificity ≥ 0.8. Ten high-confidence examples are shown. Source:
`results/monkey_spatial_148/LI_enriched_cluster_examples.tsv`.

| Species | Cluster | GI prevalence | Luminal prevalence | Mucosal prevalence | Peak LI prevalence | Region specificity |
|---|---|--:|--:|--:|--:|--:|
| *Treponema_D_succinifaciens* | C25 | 0.328 | 0.409 | 0.242 | 0.750 | 1.000 |
| *WG_1_sp900539665* | C22 | 0.328 | 0.439 | 0.210 | 0.714 | 0.936 |
| *Treponema_D_succinifaciens* | C2 | 0.289 | 0.424 | 0.145 | 0.643 | 0.963 |
| *Succinivibrio_sp000431835* | C6 | 0.266 | 0.364 | 0.161 | 0.589 | 0.960 |
| *Cryptobacteroides_sp900546445* | C1 | 0.258 | 0.364 | 0.145 | 0.571 | 0.959 |
| *Treponema_D_succinifaciens* | C15 | 0.242 | 0.394 | 0.081 | 0.554 | 1.000 |
| *Treponema_D_succinifaciens* | C5 | 0.242 | 0.333 | 0.145 | 0.554 | 1.000 |
| *UBA636_sp900546285* | C12 | 0.234 | 0.303 | 0.161 | 0.536 | 1.000 |
| *Oribacterium_sp900768305* | C4 | 0.250 | 0.379 | 0.113 | 0.536 | 0.917 |
| *Cryptobacteroides_sp000433355* | C6 | 0.258 | 0.439 | 0.065 | 0.536 | 0.880 |

**Table S8. Matched mucosal/luminal strain clusters after Bonferroni correction.** Samples were matched
within monkey, GI region and position; absent cluster calls were treated as zero within-species abundance.
Only clusters with Bonferroni q < 0.05 are shown. Source:
`results/monkey_spatial_148/mucosal_luminal_bonferroni.tsv`.

| Species | Cluster | Matched pairs | Mucosal positive | Luminal positive | p | Bonferroni q |
|---|---|--:|--:|--:|--:|--:|
| *Treponema_D_succinifaciens* | C25 | 50 | 15 | 27 | 6.0 × 10⁻⁶ | 9.6 × 10⁻⁴ |
| *UBA4372_sp900543815* | C13 | 50 | 4 | 24 | 1.8 × 10⁻⁵ | 2.8 × 10⁻³ |
| *Treponema_D_succinifaciens* | C6 | 50 | 7 | 23 | 2.7 × 10⁻⁵ | 4.1 × 10⁻³ |
| *Treponema_D_succinifaciens* | C22 | 50 | 8 | 22 | 3.1 × 10⁻⁵ | 4.7 × 10⁻³ |
| *UBA4372_sp900543815* | C5 | 50 | 7 | 32 | 4.2 × 10⁻⁵ | 6.5 × 10⁻³ |
| *Cryptobacteroides_sp900546445* | C13 | 50 | 4 | 24 | 4.4 × 10⁻⁵ | 6.7 × 10⁻³ |
| *Cryptobacteroides_sp900546445* | C3 | 50 | 4 | 27 | 8.7 × 10⁻⁵ | 1.3 × 10⁻² |
| *WG_1_sp900539665* | C20 | 50 | 16 | 30 | 1.2 × 10⁻⁴ | 1.8 × 10⁻² |
| *Cryptobacteroides_sp900546445* | C14 | 50 | 4 | 22 | 1.2 × 10⁻⁴ | 1.8 × 10⁻² |
| *Cryptobacteroides_sp000433355* | C3 | 50 | 5 | 28 | 2.1 × 10⁻⁴ | 3.1 × 10⁻² |
| *Treponema_D_succinifaciens* | C10 | 50 | 9 | 20 | 2.1 × 10⁻⁴ | 3.3 × 10⁻² |
| *Treponema_D_succinifaciens* | C12 | 50 | 4 | 21 | 2.1 × 10⁻⁴ | 3.3 × 10⁻² |
| *Treponema_D_berlinense* | C10 | 50 | 4 | 17 | 2.9 × 10⁻⁴ | 4.5 × 10⁻² |

**Table S9. Paired luminal 2bRAD/shotgun concordance on the monkey cohort.** Both modalities used the
same BcgI marker space and cluster gates. Concordance is restricted to 51 luminal samples with calls in
both modalities; differences in depth, representation and host content prevent accuracy interpretation.
Source: `results/monkey_spatial_148/paired_modality_concordance_summary.tsv`.

| Metric | Mean | Median |
|---|--:|--:|
| 2bRAD strain calls | 72.2 | 84 |
| Shotgun strain calls | 106.0 | 123 |
| Shared strain calls | 64.0 | 81 |
| Strain-call Jaccard | 0.420 | 0.605 |
| Species Jaccard | 0.735 | 0.889 |

Dominant-cluster agreement was 234/580 comparisons (40.3%).

**Table S10. Representative LI strain-cluster associations with bile acids and functions.** LI strain
composition was residualized against monkey identity and LI position, then associated with bile-acid class
composition or HUMAnN pathway composition using Spearman correlation and Benjamini–Hochberg FDR. These are
sample-level ecological associations, not gene-to-cluster assignments. Source:
`results/monkey_spatial_148/bileacid_coupling_key_results.tsv` and
`results/monkey_spatial_148/function_coupling_key_results.tsv`.

| Strain cluster | Feature | rho | p | BH q |
|---|---|--:|--:|--:|
| *Cryptobacteroides_sp000433355* C1 | Conjugated/high-mass bile acids | −0.943 | 1.8 × 10⁻¹¹ | 1.0 × 10⁻⁸ |
| *Cryptobacteroides_sp000433355* C5 | Cholate-like bile acids | −0.935 | 6.6 × 10⁻¹¹ | 1.4 × 10⁻⁸ |
| *Treponema_D_succinifaciens* C25 | Conjugated/high-mass bile acids | −0.923 | 3.6 × 10⁻¹⁰ | 5.2 × 10⁻⁸ |
| *Succinivibrio_sp000431835* C25 | Cholate-like bile acids | −0.909 | 1.9 × 10⁻⁹ | 1.4 × 10⁻⁷ |
| *UBA4372_sp900543815* C3 | Conjugated/high-mass bile acids | −0.895 | 8.1 × 10⁻⁹ | 2.7 × 10⁻⁷ |
| *UBA636_sp900546285* C11 | Dihydroxy bile acids | 0.648 | 8.2 × 10⁻⁴ | 2.5 × 10⁻³ |
| *UBA4372_sp900543815* C4 | Pyruvate fermentation to butanoate | 0.926 | 1.2 × 10⁻¹¹ | 4.0 × 10⁻¹⁰ |
| *UBA4372_sp900543815* C4 | Acetyl-CoA fermentation to butanoate II | 0.905 | 2.2 × 10⁻¹⁰ | 4.9 × 10⁻⁹ |
| *Oribacterium_sp900768305* C2 | Pyruvate fermentation to butanoate | 0.955 | 3.7 × 10⁻¹⁴ | 2.9 × 10⁻¹² |
| *Cryptobacteroides_sp000433355* C5 | Pyruvate fermentation to propanoate I | 0.981 | 1.5 × 10⁻¹⁸ | 5.6 × 10⁻¹⁶ |
| *UBA4372_sp900543815* C3 | Succinate fermentation to butanoate | −0.951 | 1.0 × 10⁻¹³ | 7.1 × 10⁻¹² |
| *Treponema_D_succinifaciens* C2 | N-acetylneuraminate degradation | 0.964 | 2.2 × 10⁻¹⁵ | 2.6 × 10⁻¹³ |
| *Treponema_D_succinifaciens* C2 | N-acetylglucosamine/mannosamine/neuraminate degradation | 0.965 | 1.8 × 10⁻¹⁵ | 2.1 × 10⁻¹³ |
| *Succinivibrio_sp000431835* C4 | N-acetylglucosamine/mannosamine/neuraminate degradation | 0.789 | 2.0 × 10⁻⁶ | 1.2 × 10⁻⁵ |

**Table S11. Strain-sharing groups in the 32-sample saliva time-course.** Rates are pairwise and were
calculated only over species with calls in both samples. Dominant-cluster sharing is the primary
StrainPhlAn-compatible metric; any-cluster sharing is the coexistence-aware secondary metric. Source:
`results/strain_sharing_20261105/saliva_strain_sharing_group_summary.tsv`.

| Comparison | Pairs | Mean any-cluster rate | Mean dominant-cluster rate | Pooled any-cluster rate |
|---|--:|--:|--:|--:|
| Same subject, between timepoints | 48 | 0.801 | 0.538 | 0.800 |
| Different subject, same timepoint | 112 | 0.574 | 0.216 | 0.578 |
| Different subject, different timepoints | 336 | 0.566 | 0.218 | 0.567 |

**Table S12. Strain-sharing groups in the 128-sample rhesus macaque GI spatial cohort.** Rates are pairwise
and restricted to GI samples. Pooled rates divide shared species-pairs by all analyzable species-pairs.
Source: `results/strain_sharing_20261105/monkey_strain_sharing_group_summary.tsv`.

| Comparison | Analyzable pairs | Mean any-cluster rate | Mean dominant-cluster rate | Pooled any-cluster rate |
|---|--:|--:|--:|--:|
| Within monkey, between positions, same region/compartment | 260 | 0.649 | 0.379 | 0.819 |
| Within monkey, mucosal/luminal, same region | 297 | 0.680 | 0.251 | 0.838 |
| Between monkeys, same age, same region/compartment | 289 | 0.535 | 0.159 | 0.728 |
| Between ages, corresponding region/compartment | 589 | 0.495 | 0.157 | 0.679 |
| Between monkeys, between regions, same compartment | 1,321 | 0.336 | 0.100 | 0.457 |
| Between compartments | 2,222 | 0.382 | 0.102 | 0.546 |

**Table S13. Old/young rhesus strain sharing at corresponding GI sites.** Corresponding-site pairs matched
GI region and compartment across O/Y monkey IDs. Region codes are retained from cohort metadata; rates use
only pairs with at least one species detected in both samples. Source:
`results/strain_sharing_20261105/monkey_between_age_corresponding_site_summary.tsv`.

| Region | Compartment | Pairs | Analyzable pairs | Mean any-cluster rate | Mean dominant-cluster rate |
|---|---|--:|--:|--:|--:|
| E | Luminal | 6 | 4 | 0.275 | 0.194 |
| K | Luminal | 2 | 2 | 0.500 | 0.125 |
| LI | Luminal | 196 | 196 | 0.870 | 0.222 |
| LI | Mucosal | 196 | 196 | 0.368 | 0.071 |
| SI | Luminal | 72 | 72 | 0.408 | 0.253 |
| SI | Mucosal | 143 | 119 | 0.149 | 0.131 |
| S | Luminal | 2 | 0 | -- | -- |
| S | Mucosal | 6 | 0 | -- | -- |

**Table S14. ATCC strain identifiers represented in the four DNA mock communities.** Each row maps the species-level colour key used in mock figures to the ATCC strain/reference genome and its 0.95-cluster representative. For MSA-1003, only staggered taxonomic abundance was provided; sequence abundance is therefore not shown. Source: `results/mock_benchmark_f26f234/mock_strain_identifiers.tsv`.

| Mock | Species | ATCC strain/reference | Truth abundance | 0.95 cluster |
|---|---|---|--:|---|
| MSA1002 | *Acinetobacter baumannii* | ATCC Acinetobacter baumannii ATCC 17978 | 0.06079 | C0 |
| MSA1002 | *Bacillus pacificus* | ATCC Bacillus pacificus ATCC 10987 | 0.08117 | C12 |
| MSA1002 | *Bifidobacterium adolescentis* | ATCC Bifidobacterium adolescentis ATCC 15703 | 0.03116 | C18 |
| MSA1002 | *Cereibacter sphaeroides* | ATCC Cereibacter sphaeroides ATCC 17029 | 0.06694 | C24 |
| MSA1002 | *Clostridium beijerinckii* | ATCC Clostridium beijerinckii ATCC 35702 | 0.0896 | C28 |
| MSA1002 | *Cutibacterium acnes* | ATCC Cutibacterium acnes ATCC 11828 | 0.03725 | C35 |
| MSA1002 | *Deinococcus radiodurans* | ATCC Deinococcus radiodurans ATCC BAA 816 | 0.04892 | C40 |
| MSA1002 | *Enterococcus faecalis* | ATCC Enterococcus faecalis ATCC 47077 | 0.04084 | C41 |
| MSA1002 | *Escherichia coli* | ATCC Escherichia coli ATCC 700926 | 0.06924 | C47 |
| MSA1002 | *Helicobacter pylori* | ATCC Helicobacter pylori ATCC 700392 | 0.02487 | C56 |
| MSA1002 | *Lactobacillus gasseri* | ATCC Lactobacillus gasseri ATCC 33323 | 0.02765 | C58 |
| MSA1002 | *Neisseria meningitidis* | ATCC Neisseria meningitidis ATCC BAA 335 | 0.03347 | C78 |
| MSA1002 | *Phocaeicola vulgatus* | ATCC Phocaeicola vulgatus ATCC 8482 | 0.077 | C84 |
| MSA1002 | *Porphyromonas gingivalis* | ATCC Porphyromonas gingivalis ATCC 33277 | 0.03579 | C90 |
| MSA1002 | *Pseudomonas paraeruginosa* | ATCC Pseudomonas paraeruginosa ATCC 9027 | 0.09507 | C100 |
| MSA1002 | *Schaalia odontolytica* | ATCC Schaalia odontolytica ATCC 17982 | 0.0357 | C105 |
| MSA1002 | *Staphylococcus aureus* | ATCC Staphylococcus aureus ATCC BAA 1556 | 0.0436 | C110 |
| MSA1002 | *Staphylococcus epidermidis* | ATCC Staphylococcus epidermidis ATCC 12228 | 0.03842 | C115 |
| MSA1002 | *Streptococcus agalactiae* | ATCC Streptococcus agalactiae ATCC BAA 611 | 0.03221 | C121 |
| MSA1002 | *Streptococcus mutans* | ATCC Streptococcus mutans ATCC 700610 | 0.0303 | C129 |
| MSA1003 | *Acinetobacter baumannii* | ATCC Acinetobacter baumannii ATCC 17978 | -- | C0 |
| MSA1003 | *Bacillus pacificus* | ATCC Bacillus pacificus ATCC 10987 | -- | C12 |
| MSA1003 | *Bifidobacterium adolescentis* | ATCC Bifidobacterium adolescentis ATCC 15703 | -- | C18 |
| MSA1003 | *Cereibacter sphaeroides* | ATCC Cereibacter sphaeroides ATCC 17029 | -- | C24 |
| MSA1003 | *Clostridium beijerinckii* | ATCC Clostridium beijerinckii ATCC 35702 | -- | C28 |
| MSA1003 | *Cutibacterium acnes* | ATCC Cutibacterium acnes ATCC 11828 | -- | C35 |
| MSA1003 | *Deinococcus radiodurans* | ATCC Deinococcus radiodurans ATCC BAA 816 | -- | C40 |
| MSA1003 | *Enterococcus faecalis* | ATCC Enterococcus faecalis ATCC 47077 | -- | C41 |
| MSA1003 | *Escherichia coli* | ATCC Escherichia coli ATCC 700926 | -- | C47 |
| MSA1003 | *Helicobacter pylori* | ATCC Helicobacter pylori ATCC 700392 | -- | C56 |
| MSA1003 | *Lactobacillus gasseri* | ATCC Lactobacillus gasseri ATCC 33323 | -- | C58 |
| MSA1003 | *Neisseria meningitidis* | ATCC Neisseria meningitidis ATCC BAA 335 | -- | C78 |
| MSA1003 | *Phocaeicola vulgatus* | ATCC Phocaeicola vulgatus ATCC 8482 | -- | C84 |
| MSA1003 | *Porphyromonas gingivalis* | ATCC Porphyromonas gingivalis ATCC 33277 | -- | C90 |
| MSA1003 | *Pseudomonas paraeruginosa* | ATCC Pseudomonas paraeruginosa ATCC 9027 | -- | C100 |
| MSA1003 | *Schaalia odontolytica* | ATCC Schaalia odontolytica ATCC 17982 | -- | C105 |
| MSA1003 | *Staphylococcus aureus* | ATCC Staphylococcus aureus ATCC BAA 1556 | -- | C110 |
| MSA1003 | *Staphylococcus epidermidis* | ATCC Staphylococcus epidermidis ATCC 12228 | -- | C115 |
| MSA1003 | *Streptococcus agalactiae* | ATCC Streptococcus agalactiae ATCC BAA 611 | -- | C121 |
| MSA1003 | *Streptococcus mutans* | ATCC Streptococcus mutans ATCC 700610 | -- | C129 |
| MSA1005 | *Acinetobacter johnsonii* | ATCC Acinetobacter johnsonii ATCC 17909 | 0.227 | C6 |
| MSA1005 | *Corynebacterium striatum* | ATCC Corynebacterium striatum ATCC 6940 | 0.1841 | C29 |
| MSA1005 | *Cutibacterium acnes* | ATCC Cutibacterium acnes ATCC 11828 | 0.1564 | C35 |
| MSA1005 | *Micrococcus luteus* | ATCC Micrococcus luteus ATCC 4698 | 0.1541 | C72 |
| MSA1005 | *Staphylococcus epidermidis* | ATCC Staphylococcus epidermidis ATCC 12228 | 0.1613 | C115 |
| MSA1005 | *Streptococcus mitis* | ATCC Streptococcus mitis ATCC 49456 | 0.117 | C127 |
| MSA1007 | *Gardnerella vaginalis* | ATCC Gardnerella vaginalis ATCC 14019 | 0.1576 | C50 |
| MSA1007 | *Lactobacillus gasseri* | ATCC Lactobacillus gasseri ATCC 33323 | 0.1753 | C58 |
| MSA1007 | *Lactobacillus jensenii* | ATCC Lactobacillus jensenii ATCC 25258 | 0.1577 | C64 |
| MSA1007 | *Metamycoplasma hominis* | ATCC Metamycoplasma hominis ATCC 23114 | 0.0633 | C68 |
| MSA1007 | *Prevotella bivia* | ATCC Prevotella bivia ATCC 29303 | 0.242 | C96 |
| MSA1007 | *Streptococcus agalactiae* | ATCC Streptococcus agalactiae ATCC BAA 611 | 0.2042 | C121 |

The source table also contains representative-genome accessions and taxonomic abundance.


**Table S15. Composition of the 2,025-sample single-species simulation pool.** Each species contained 2, 3
or 5 co-present truth strains from the same or different 0.95 clusters at five per-strain depths and five
replicates. Cells show sample counts. *M. tuberculosis* contributed only same-cluster 5-strain samples and
*P. copri* contributed only different-cluster 5-strain samples. Source:
`results/sim2025_current_summary/sample_composition.tsv`.

| Species | Design | Truth strains | 0.5× | 1× | 3× | 5× | 10× | Total |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| *Akkermansia muciniphila* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Clostridioides difficile* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Cutibacterium acnes* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Escherichia coli* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Fusobacterium nucleatum* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Klebsiella pneumoniae* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Lactiplantibacillus plantarum* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Mycobacterium tuberculosis* | same only | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Phocaeicola dorei* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Porphyromonas gingivalis* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Prevotella copri* | different only | 2/3 | 5 each | 5 each | 5 each | 5 each | 5 each | 100 |
| *Salmonella enterica* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Staphylococcus aureus* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Staphylococcus epidermidis* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| *Streptococcus pneumoniae* | same/different | 2/3/5 | 5 each | 5 each | 5 each | 5 each | 5 each | 150 |
| **Total** |  |  | **405** | **405** | **405** | **405** | **405** | **2,025** |
