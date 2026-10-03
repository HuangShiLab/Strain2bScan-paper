# Benchmark datasets, results, and conclusions (summary)

All simulated benchmarks use the **corrected enzyme table** (Fast2bRAD-M tag lengths, single-strand
scan + canonical hash) — tags are interoperable with Fast2bRAD-M / `2bRADExtraction.pl`. Reference
genomes are real (NCBI/ENA, pinned in `data/accessions/`); simulated reads are 150 bp, error-free,
50% reverse-complemented, log-normal per-strain depth ≥1× (except the depth sweep). Default enzyme
set = 14 (BcgI + 13), unless noted. Metrics: precision/recall at 0.01 presence, Bray–Curtis abundance
error, wall-clock time, peak RSS. 16-core arm64 macOS.

## Master table

| # | Experiment | Reference panel (real genomes) | Samples | What it tests | Headline result | Conclusion |
|---|---|---|---|---|---|---|
| 1 | **Performance / head-to-head** | *C. acnes* 64 | 5 strain mixtures | per-sample speed & memory vs StrainScan | **0.86 s / 78 MB** vs StrainScan 7.06 s / 828 MB | **~8× faster, ~11× lighter** |
| 2 | **Cross-species** | *C. acnes* 64, *S. aureus* 60, *S. epidermidis* 60 | 5 each | does accuracy generalize across species | **precision 1.0** all three; recall 0.75 / 0.79 / 1.0; Bray–Curtis 0.24 / 0.33 / 0.02 | profiles cleanly across species; no over-detection |
| 3 | **Species expansion** (vs *real* StrainScan on its own DBs) | *A. muciniphila* 40, *P. copri* 40, *M. tuberculosis* 40 | 5 each | head-to-head accuracy on 3 more species | Strain2bScan **P=1.0** all; StrainScan **DNF on *M. tuberculosis*** (>3.3 h, >25 GB) | matches StrainScan precision, ~17–23× faster; StrainScan has a scalability ceiling on near-clonal species |
| 4 | **Reference panel size** | *P. copri* nested 40 ⊂ 80 ⊂ 112 | same 5 across sizes | robustness to DB size / clustering correctness | **P=1.0** at every size; clusters 23 / 43 / **51** | robust to panel size; clustering (112→51) exactly matches StrainScan's own *P. copri* DB |
| 5 | **Multi-species community** | **55 species × ~4 strains = 218 genomes** | 30 (12 species each) | scaling to a complex community | species gradient **flat** (~3.3 s, 10→55 sp); samples **linear** (~2.6 s); gate precision **0.96→1.0** | digest-once architecture scales; **~132× at 100 samples, ~146× at 200/500** vs per-species StrainScan |
| 6 | **Depth sensitivity** | *C. acnes* 64 (1 target strain) | 0.1× / 0.5× / 1× / 5× | low-depth detection vs StrainScan | detects at **0.5×**, misses 0.1× — **matches StrainScan** | no low-depth penalty |
| 7 | **Enzyme count (the 2bRAD knob)** | 14 resolvable species of the 15-species pool | 5, ladder 1/2/4/8/14 enzymes | more enzymes → better? optimum across species | **P=1.0 at every count incl. BcgI-alone**; median recall 0.50→0.72→1.00 | BcgI alone works; **~4-enzyme sweet spot**; more enzymes add recall/marker density, not clusters |
| 8 | **Reference-genome quality** | 15-species pool | fixed reads, degraded references | effect of incomplete/contaminated references | Jaccard precision 1.0→0.84→0.71, recall 0.96→0.80→0.74 at 100→90→70%; `--containment` restores precision/recall to ~0.92 at 90% | robust to moderate degradation; containment fixes subset fragmentation |
| 9 | **2bRAD-M (BcgI) data path** | *P. copri* 40 (BcgI-only DB) | reads = the BcgI tags (50% RC) | does the native 2bRAD-M workflow work | truth clusters detected exactly, **0 false positives** | native BcgI 2bRAD-M libraries profile correctly (unique capability) |

## Real-data benchmarks (not simulated)

| # | Experiment | Data | What it tests | Headline result | Conclusion |
|---|---|---|---|---|---|---|
| 10 | **ATCC DNA mocks (native 2bRAD)** | MSA-1002/1003/1005/1007 whole-cell mocks, native BcgI 2bRAD-M; host-DNA ladder and DNA-input titration | real DNA, known strain truth, host/low-input robustness | precision 1.0, recall 19/20 at 0.01 ng, full 20/20 recall at ≥0.1 ng; precision/recall ≈0.95–1.0 at 99 % host | native 2bRAD-M strain profiling works on real DNA under heavy host contamination |
| 11 | **Real saliva** | 8 subjects × 4 timepoints, native BcgI 2bRAD-M + paired WMS (PRJNA1131785) | individual discrimination, temporal stability, 2bRAD↔shotgun concordance | strain host-ID **100 %** vs species **78.1 %**; PERMANOVA R² 0.757 vs 0.755; native 2bRAD recovers 128–163 additional strains/sample | strain signatures are individual-specific and stable; native 2bRAD-M recovers low-abundance strains host-limited WMS misses |
| 12 | **Public WGS cohorts** | PRJNA288562, PRJNA1517970, PRJNA1191223, PRJNA1191225 subsets | real-world application: longitudinal, multi-site, isolate-derived panels | 22 libraries profiled in 100.6 s / 343 MiB; cohort-specific isolate panel 6/6 self-calls; persistent *B. bifidum* tracked across 3 weeks | Strain2bScan runs on real WGS cohorts; cohort-specific panels recover expected isolates |

## Simulated-benchmark caveat

The controlled accuracy, scaling, depth-sensitivity, enzyme-sweep and reference-quality experiments
above use **simulated reads from the same genomes that build the DB** (closed-world, error-free).
They rigorously test the algorithm but do not test: real sequencing error, unknown/absent strains
(open-world false-positive control), real community complexity/abundance skew, or biological utility
on real samples. The real-data benchmarks (ATCC mocks, saliva, public WGS cohorts) close this gap
for the specific claims they address.

## Cross-cutting conclusions

1. **Accurate:** precision 1.0 across every species and condition tested; abundance error low
   (Bray–Curtis 0.02–0.33). Recall is the species-dependent axis and tracks *genuine*
   intra-species diversity (near-clonal *M. tuberculosis* → low; diverse species → 0.9–1.0).
2. **Fast & light:** ~8× faster / ~11× lighter per sample than StrainScan, and — because a sample
   is digested once and matched against every species DB — **~130–150× faster at community scale**
   (55 species × hundreds of samples), where StrainScan must re-run per species.
3. **Interoperable & 2bRAD-native:** tags match Fast2bRAD-M exactly (species layer → strain layer
   compose), BcgI alone resolves strains, and native BcgI 2bRAD-M data profiles correctly.
4. **Honest limits:** near-clonal species are intrinsically coarse (recall, not precision); very
   degraded references (<~70% completeness) break detection; real metagenome strain calling can be
   limited by BcgI marker density (see `results/real_data_strain_benchmark.md`).
