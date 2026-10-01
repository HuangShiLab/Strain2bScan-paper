# Strain2bScan; Manuscript tables (simulated head-to-head)

**Table 1. Strain2bScan vs StrainScan on the 15-species simulated benchmark, per species.** Single-species accuracy is the median over depth-matched paired samples (2/3/5-strain mixtures across the 0.5–10× ladder), each tool scored in its own cluster space. Database build cost is per species (Strain2bScan native, arm64; StrainScan `linux/amd64` under emulation). n = genomes in the pool.

| Species | n | S2B P | S2B R | S2B F1 | SS P | SS R | SS F1 | S2B build | SS build |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| *A. muciniphila* | 50 | 1.000 | 0.800 | 0.889 | 0.833 | 0.667 | 0.667 | 2.8 s / 0.16 GB | 17 min / 15 GB |
| *C. difficile* | 47 | 1.000 | 0.800 | 0.889 | 0.800 | 0.667 | 0.667 | 3.2 s / 0.28 GB | 30 min / 14 GB |
| *C. acnes* | 43 | 1.000 | 0.667 | 0.800 | 1.000 | 0.667 | 0.800 | 1.5 s / 0.16 GB | 11 min / 8 GB |
| *E. coli* | 50 | 1.000 | 1.000 | 1.000 | 1.000 | 0.600 | 0.667 | 4.3 s / 0.33 GB | 43 min / 28 GB |
| *F. nucleatum* | 25 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 0.7 s / 0.10 GB | 5 min / 8 GB |
| *L. plantarum* | 50 | 1.000 | 1.000 | 1.000 | 1.000 | 0.600 | 0.750 | 5.0 s / 0.23 GB | 24 min / 17 GB |
| *M. tuberculosis* | 29 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.6 s / 0.25 GB | 16 min / 15 GB |
| *P. dorei* | 15 | 1.000 | 0.667 | 0.800 | 0.833 | 0.583 | 0.667 | 1.0 s / 0.25 GB | 5 min / 13 GB |
| *P. gingivalis* | 43 | 1.000 | 0.667 | 0.800 | 1.000 | 0.667 | 0.800 | 2.3 s / 0.15 GB | 19 min / 10 GB |
| *P. copri* | 19 | 1.000 | 1.000 | 1.000 | 1.000 | 0.667 | 0.800 | 1.8 s / 0.18 GB | 8 min / 25 GB |
| *S. enterica* | 45 | 1.000 | 1.000 | 1.000 | 0.833 | 0.800 | 0.800 | 4.2 s / 0.36 GB | 34 min / 16 GB |
| *S. aureus* | 48 | 1.000 | 0.667 | 0.800 | 1.000 | 0.500 | 0.667 | 1.7 s / 0.16 GB | 17 min / 12 GB |
| *S. epidermidis* | 50 | 1.000 | 0.667 | 0.800 | 1.000 | 0.600 | 0.750 | 1.6 s / 0.13 GB | 15 min / 8 GB |
| *S. pneumoniae* | 50 | 1.000 | 0.667 | 0.800 | 1.000 | 0.600 | 0.750 | 1.5 s / 0.15 GB | 11 min / 11 GB |

**Median (14 resolvable species):** Strain2bScan P 1.00 / R 0.80 / F1 0.89; StrainScan P 1.00 / R 0.67 / F1 0.75. Build speed-up 249–614×; build memory 43–138× lighter.

**Table 2. Accuracy and per-sample cost vs sequencing depth (single-species, 14 species, 204 paired samples; medians).** Profile time/memory in the same emulated container (Strain2bScan `linux/amd64` vs StrainScan).

| Depth (×) | n | S2B P/R/F1 | SS P/R/F1 | S2B time/mem | SS time/mem |
|--:|--:|:--:|:--:|:--:|:--:|
| 0.5 | 40 | 1.000/0.500/0.667 | 1.000/0.667/0.800 | 0.16 s / 21 MB | 5.3 s / 831 MB |
| 1 | 41 | 1.000/0.667/0.800 | 1.000/0.333/0.500 | 0.24 s / 30 MB | 3.5 s / 831 MB |
| 3 | 41 | 1.000/1.000/1.000 | 1.000/0.600/0.750 | 0.56 s / 60 MB | 4.2 s / 831 MB |
| 5 | 41 | 1.000/1.000/1.000 | 1.000/0.667/0.800 | 0.89 s / 91 MB | 4.8 s / 831 MB |
| 10 | 41 | 1.000/1.000/1.000 | 1.000/1.000/0.889 | 1.66 s / 161 MB | 6.9 s / 832 MB |

**Table 3. Multi-species community profiling (4 samples/depth, matched to the 14 species with StrainScan databases; medians).** Strain2bScan profiles each community in one digest-once pass; StrainScan (no multi-species mode) profiles once per species, so its cost is the sum over 14 databases.

| Community depth | S2B P/R/F1 | SS P/R/F1 | S2B time/mem | SS time/mem |
|---|:--:|:--:|:--:|:--:|
| low | 0.926/0.678/0.782 | 0.898/0.767/0.827 | 1.0 s / 311 MB | 100 s / 1112 MB |
| med | 0.863/0.853/0.869 | 0.911/0.856/0.895 | 4.3 s / 670 MB | 228 s / 1696 MB |
| high | 0.773/0.872/0.819 | 0.814/0.972/0.895 | 8.7 s / 1119 MB | 398 s / 2028 MB |

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
