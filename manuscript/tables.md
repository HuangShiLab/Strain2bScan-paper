# Strain2bScan — Manuscript tables (simulated head-to-head)

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
