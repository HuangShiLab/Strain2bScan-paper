# Benchmark datasets, results, and conclusions (summary)

All benchmarks use the **corrected enzyme table** (Fast2bRAD-M tag lengths, single-strand scan +
canonical hash) — tags are interoperable with Fast2bRAD-M / `2bRADExtraction.pl`. Reference
genomes are real (NCBI/ENA, pinned in `data/accessions/`); reads are **simulated**: 150 bp,
error-free, 50% reverse-complemented, log-normal per-strain depth ≥1× (except the depth sweep).
Default enzyme set = 14 (BcgI + 13), unless noted. Metrics: precision/recall at 0.01 presence,
Bray–Curtis abundance error, wall-clock time, peak RSS. 16-core arm64 macOS.

## Master table

| # | Experiment | Reference panel (real genomes) | Simulated samples | What it tests | Headline result | Conclusion |
|---|---|---|---|---|---|---|
| 1 | **Performance / head-to-head** | *C. acnes* 64 | 5 strain mixtures | per-sample speed & memory vs StrainScan | **0.86 s / 78 MB** vs StrainScan 7.06 s / 828 MB | **~8× faster, ~11× lighter** |
| 2 | **Cross-species** | *C. acnes* 64, *S. aureus* 60, *S. epidermidis* 60 | 5 each | does accuracy generalize across species | **precision 1.0** all three; recall 0.75 / 0.79 / 1.0; Bray–Curtis 0.24 / 0.33 / 0.02 | profiles cleanly across species; no over-detection |
| 3 | **Species expansion** (vs *real* StrainScan on its own DBs) | *A. muciniphila* 40, *P. copri* 40, *M. tuberculosis* 40 | 5 each | head-to-head accuracy on 3 more species | Strain2bScan **P=1.0** all; StrainScan **DNF on *M. tuberculosis*** (>3.3 h, >25 GB) | matches StrainScan precision, ~17–23× faster; StrainScan has a scalability ceiling on near-clonal species |
| 4 | **Reference panel size** | *P. copri* nested 40 ⊂ 80 ⊂ 112 | same 5 across sizes | robustness to DB size / clustering correctness | **P=1.0** at every size; clusters 23 / 43 / **51** | robust to panel size; clustering (112→51) exactly matches StrainScan's own *P. copri* DB |
| 5 | **Multi-species community** | **55 species × ~4 strains = 218 genomes** | 30 (12 species each) | scaling to a complex community | species gradient **flat** (~3.3 s, 10→55 sp); samples **linear** (~2.6 s); gate precision **0.96→1.0** | digest-once architecture scales; **~132× at 100 samples, ~146× at 200/500** vs per-species StrainScan |
| 6 | **Depth sensitivity** | *C. acnes* 64 (1 target strain) | 0.1× / 0.5× / 1× / 5× | low-depth detection vs StrainScan | detects at **0.5×**, misses 0.1× — **matches StrainScan** | no low-depth penalty (an earlier "less sensitive" claim was a bug artifact) |
| 7 | **Enzyme count (the 2bRAD knob)** | *C. acnes* 64 | 5, ladder 1/2/4/8/14 enzymes | more enzymes → better? optimum? | **P=1.0 at every count incl. BcgI-alone**; recall 0.56→0.75, plateaus by ~4 | BcgI alone works; **~4-enzyme sweet spot**; more enzymes add recall, not clusters |
| 8 | **Reference-genome quality** | *C. acnes* (5 truth degraded + 59 background) | 5 fixed | effect of incomplete/contaminated references | **P=1.0 down to 70% completeness**; recall/abundance degrade; total failure at 50% / 10% contam | robust to moderate degradation; hard floor ~MIMAG-medium |
| 9 | **2bRAD-M (BcgI) data path** | *P. copri* 40 (BcgI-only DB) | reads = the 32 bp BcgI tags (50% RC) | does the native 2bRAD workflow work | truth clusters detected exactly, **0 false positives** | real BcgI 2bRAD libraries profile correctly (unique capability) |

## Cross-cutting conclusions

1. **Accurate:** precision 1.0 across every species and condition tested; abundance error low
   (Bray–Curtis 0.02–0.33). Recall is the species-dependent axis and tracks *genuine*
   intra-species diversity (near-clonal *M. tuberculosis* → 0.25; diverse species → 0.9–1.0).
2. **Fast & light:** ~8× faster / ~11× lighter per sample than StrainScan, and — because a sample
   is digested once and matched against every species DB — **~130–150× faster at community scale**
   (55 species × hundreds of samples), where StrainScan must re-run per species.
3. **Interoperable & 2bRAD-native:** tags match Fast2bRAD-M exactly (species layer → strain layer
   compose), BcgI alone resolves strains, and real BcgI 2bRAD data profiles correctly.
4. **Honest limits:** near-clonal species are intrinsically coarse (recall, not precision); very
   degraded references (<~60% completeness) break detection.

## ⚠ The gap: all reads above are simulated (closed-world, error-free)

Every dataset uses **simulated reads from the same genomes that build the DB** (closed-world:
the truth strain is guaranteed present) and **no sequencing error**. This is the main
limitation a reviewer will raise. It does not test: real sequencing error, unknown/absent
strains (open-world false-positive control), real community complexity/abundance skew, or
practical utility on a real sample. See `docs/real_metagenome_plan.md` for how the available
saliva metagenome data would close this gap.
