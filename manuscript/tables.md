# Strain2bScan; Manuscript tables (simulated head-to-head)

**Table 1. Current Strain2bScan regression versus the reproducible StrainScan rerun on the 15-species simulated benchmark.** Accuracy is the median over completed run-level paired samples (different-cluster mixtures for 14 species; same-cluster mixtures for near-clonal *M. tuberculosis*; k = 2/3/5; depths 0.5–10×); each tool was scored in its own 0.95-similarity cluster space. Cluster-count columns are each tool's own 0.95 partitions of the common species pool. Current Strain2bScan values are from source commit `6ae760c2f3bd2f39efc0e6aa06c53a84366cb16e`; current source: `results/strainscan_rerun/current_paired_persample.tsv`. Build cost is the archived timing benchmark (Strain2bScan native arm64; StrainScan `linux/amd64` under emulation). *S. enterica* produced no completed StrainScan calls in the rerun; *K. pneumoniae*, which did not finish in the archived timing run, was successfully built and profiled in the reproducible rerun.

| Species | Paired n | Pool genomes | S2B clusters | SS clusters | S2B P | S2B R | S2B F1 | SS P | SS R | SS F1 | S2B build | SS build |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| *A. muciniphila* | 14 | 50 | 24 | 24 | 0.900 | 0.900 | 0.900 | 1.000 | 1.000 | 1.000 | 2.8 s / 0.16 GB | 17 min / 16 GB |
| *C. difficile* | 15 | 47 | 27 | 24 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 3.2 s / 0.28 GB | 30 min / 14 GB |
| *C. acnes* | 15 | 43 | 8 | 8 | 1.000 | 0.800 | 0.800 | 1.000 | 0.800 | 0.889 | 1.5 s / 0.17 GB | 11 min / 8 GB |
| *E. coli* | 15 | 50 | 30 | 27 | 1.000 | 1.000 | 1.000 | 1.000 | 0.800 | 0.889 | 4.3 s / 0.34 GB | 43 min / 28 GB |
| *F. nucleatum* | 15 | 25 | 19 | 19 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 0.7 s / 0.10 GB | 5 min / 8 GB |
| *K. pneumoniae* | 15 | 47 | 28 | 27 | 0.333 | 0.333 | 0.286 | 1.000 | 0.800 | 0.889 | 5.1 s / 0.40 GB | DNF (timing run) |
| *L. plantarum* | 15 | 50 | 39 | 39 | 1.000 | 1.000 | 1.000 | 1.000 | 0.800 | 0.889 | 5.0 s / 0.24 GB | 24 min / 18 GB |
| *M. tuberculosis* | 15 | 29 | 1 | 1 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.6 s / 0.25 GB | 16 min / 15 GB |
| *P. dorei* | 10 | 15 | 10 | 10 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.0 s / 0.25 GB | 5 min / 13 GB |
| *P. gingivalis* | 15 | 43 | 27 | 31 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 2.3 s / 0.15 GB | 19 min / 11 GB |
| *P. copri* | 15 | 19 | 19 | 19 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.8 s / 0.19 GB | 8 min / 26 GB |
| *S. enterica* | 0 | 45 | 19 | 12 | -- | -- | -- | -- | -- | -- | 4.2 s / 0.37 GB | 34 min / 16 GB |
| *S. aureus* | 15 | 48 | 13 | 12 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.7 s / 0.17 GB | 17 min / 12 GB |
| *S. epidermidis* | 15 | 50 | 9 | 9 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.6 s / 0.14 GB | 15 min / 8 GB |
| *S. pneumoniae* | 15 | 50 | 38 | 38 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.5 s / 0.16 GB | 11 min / 11 GB |

**Median across 14 species with completed paired runs:** Strain2bScan P 1.000 / R 1.000 / F1 1.000; StrainScan P 1.000 / R 0.667 / F1 0.800. Species medians conceal a precision/recall trade-off for *K. pneumoniae*. Mean paired differences and 95% intervals are in `results/strainscan_rerun/current_uncertainty_summary.tsv`.

**Table S4. Failure-aware sensitivity for all 225 matched single-species entries.** Absent StrainScan final reports were retained and scored as zero detection. Mean paired differences use the same species-cluster bootstrap as the completed-pair analysis. Source: `results/strainscan_rerun/current_single_persample.tsv` and `results/strainscan_rerun/current_uncertainty_summary.tsv`.

| Analysis | n | No-calls | Mean recall difference (95% CI) | Mean F1 difference (95% CI) |
|---|--:|--:|:--:|:--:|
| Completed pairs | 204 | 0 | +0.063 (−0.037 to 0.133) | +0.005 (−0.091 to 0.069) |
| Failure as absence | 225 | 21 | +0.127 (−0.020 to 0.176) | +0.076 (−0.075 to 0.120) |

**Table 2. Accuracy and archived per-sample cost vs sequencing depth (single-species rerun; medians).** Current Strain2bScan accuracy is from `results/strainscan_rerun/current_single_persample.tsv` on the 204 completed StrainScan pairs. Cost comes from the archived same-environment emulation subset (Strain2bScan `linux/amd64` vs StrainScan).

| Depth (×) | n | S2B P/R/F1 | SS P/R/F1 | S2B time/mem | SS time/mem |
|--:|--:|:--:|:--:|:--:|:--:|
| 0.5 | 40 | 1.000/0.600/0.708 | 1.000/0.667/0.800 | 0.16 s / 21 MB | 5.30 s / 831 MB |
| 1 | 41 | 1.000/1.000/1.000 | 1.000/1.000/1.000 | 0.24 s / 30 MB | 3.48 s / 831 MB |
| 3 | 41 | 1.000/1.000/1.000 | 1.000/0.600/0.750 | 0.56 s / 60 MB | 4.17 s / 831 MB |
| 5 | 41 | 1.000/1.000/1.000 | 1.000/0.667/0.800 | 0.89 s / 91 MB | 4.79 s / 831 MB |
| 10 | 41 | 1.000/1.000/1.000 | 1.000/1.000/1.000 | 1.66 s / 161 MB | 6.92 s / 832 MB |

Strain2bScan reached median recall 1.0 from 1× onward; StrainScan was non-monotonic and did not exceed 0.800 median F1 at 0.5× or intermediate depths. The depth table aggregates heterogeneous strain mixtures and is not a causal depth-onset estimate.

**Table 3. Multi-species community accuracy (reproducible rerun; medians) and archived profiling cost.** Current Strain2bScan values are from `results/strainscan_rerun/current_multi_persample.tsv`. Four samples per depth. Strain2bScan profiles each community in one digest-once pass; StrainScan has no multi-species mode, so archived cost is the sum over per-species runs.

| Community depth | S2B P/R/F1 | SS P/R/F1 | S2B time/mem | SS time/mem |
|---|:--:|:--:|:--:|:--:|
| low | 0.748/0.667/0.699 | 0.938/0.702/0.798 | 1.0 s / 311 MB | 100 s / 1112 MB |
| med | 0.838/0.838/0.837 | 0.954/0.903/0.927 | 4.3 s / 670 MB | 228 s / 1696 MB |
| high | 0.806/0.806/0.806 | 0.972/0.958/0.965 | 8.7 s / 1119 MB | 398 s / 2028 MB |

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
| StrainScan v1.0.14 | Yes | Shotgun reads against per-species reference databases | Strain clusters | Direct reference-guided comparator on shotgun reads (Fig 10–11, Supplementary Fig S5; Tables 1–3). |
| inStrain v1.9.1 | Yes | Shotgun reads aligned to a dereplicated reference | Strain and SNV | Direct read-level shotgun comparator (Fig 11 and Fig S4). |
| StrainGE / StrainGST / StrainGR | StrainGST only | Shotgun WMS reads against StrainGE 0.90-reference databases; native BcgI 2bRAD is not documented input | Strain level, reference guided | Related shotgun comparator rerun in Table 9; StrainGR was not run. |
| 2bRAD-M / Fast2bRAD-M | No | Native BcgI 2bRAD reads | Mainly species | Same restriction-assay family; provides the species-level upstream context for Strain2bScan. |
| sylph | No | Shotgun sketches | Mainly species and ANI | Context for rapid species-level profiling; not a strain-resolved benchmark here. |
| StrainPhlAn | No | Shotgun marker genes | Strain types | Context for marker-gene strain profiling; not rerun here. |

**Table 9. StrainGST (StrainGE 1.3.9) rerun accuracy and efficiency.** Simulations are the matched
`_rep1_` subset and all 12 multi-species communities. Mocks are the primary WMS replicate used in
Figure 11; MSA-1002 is the 99%-host sample. Accuracy is in StrainGE's 0.90-reference space; mock calls
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

**Table S3. Exploratory native BcgI 2bRAD-M oral-library profiling.** Four libraries from the
PRJNA1131785 exploratory oral cohort were profiled against the same 19-species panel used for the saliva
analysis. No case/control labels were available in public metadata. Values are from
`results/clinical_exploratory.tsv`.

| Sample | Distinct markers | Species resolved | Strain calls | Runtime (s) | Top species |
|---|--:|--:|--:|--:|---|
| S_0313206 | 271,322 | 15 | 115 | 2.3 | *Neisseria subflava*; *Haemophilus parainfluenzae*; *Actinomyces odontolyticus* |
| S_1413231 | 288,019 | 17 | 129 | 1.8 | *Neisseria subflava*; *Streptococcus mitis*; *Actinomyces odontolyticus* |
| S_8123213 | 370,379 | 16 | 149 | 2.0 | *Actinomyces odontolyticus*; *Neisseria subflava*; *Haemophilus parainfluenzae* |
| S_8221232 | 488,180 | 17 | 158 | 1.8 | *Actinomyces odontolyticus*; *Neisseria subflava*; *Rothia mucilaginosa* |

**Table S5. Public-run manifest and exploratory profiling outcomes.** SRA spots, bases and instrument metadata were retrieved from NCBI on 2026-10-11. Marker and call counts are from the generic MSA screen; body-site-panel calls are shown for PRJNA1517970. P08 run-to-week mapping follows SRA submitter titles. Source: `results/realworld_cohort_benchmark/public_run_manifest.tsv`.

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
