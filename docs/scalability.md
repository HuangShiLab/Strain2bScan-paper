# Scalability & parallelism — Strain2bScan

Evidence for the scale/parallelism claim. Engineering added: zero-dependency multi-threading
(`std::thread::scope`, env knob `STRAIN2BSCAN_THREADS`), MinHash-sketch clustering
(O(n²·m) → O(n²·k)), and a `multi-profile` mode (digest the sample once, match all species
DBs in parallel).

## 1. Parallel speedup (real C. acnes, 64 genomes)
| stage | 1 thread | 16 threads | speedup |
|---|---|---|---|
| DB build (cluster) | 24.8 s | 2.5 s | **9.7×** |
| profile (200k reads) | 3.67 s | 0.50 s | **7.3×** |

→ Strain2bScan profile **0.50 s** vs StrainScan **7.2 s** = **~14× faster**, ~8× less memory.

## 2. MinHash clustering: scalable *and* faithful
Replaces exact all-pairs Jaccard (O(n²·m)) with bottom-k MinHash sketches (O(n²·k), k=2000).
On the 64 real C. acnes genomes, MinHash yields a partition **identical to exact** (60
clusters; same non-singleton groups). So the scalable path does not change results.

## 3. DB build scaling (synthetic genomes, 16 threads)
| genomes | 50 | 100 | 200 | 400 |
|---|---|---|---|---|
| build time | 0.26 s | 0.53 s | 1.56 s | 5.07 s |
| peak RSS | 24 MB | 41 MB | 64 MB | 111 MB |

400 genomes also builds in 15.3 s on **1 thread** → 3.0× from threading here (clustering is
the serial-ish O(n²) tail; digestion parallelizes fully). Memory grows slowly.

## 4. Multi-species throughput — the architectural win
`multi-profile`: digest the sample **once**, then match the shared tag counts against every
per-species cluster DB **in parallel**. Synthetic 80-species panel (5 strains each), one
sample containing strains from 40 of the 80 species at ~2×:

| | result |
|---|---|
| time (16 threads) | **0.11 s** |
| time (1 thread) | 0.74 s |
| peak RSS | 33 MB |
| species detected | **40 / 40 present, 0 / 40 absent (no false positives)** |

The contrast with a full-k-mer profiler is structural: StrainScan re-counts k-mers (jellyfish,
~6–7 s) **per run**, so S species ≈ S × 6.5 s. Strain2bScan pays the read-digestion cost
**once** (~0.1–0.5 s) and each additional species is a cheap hash-set match — so per-species
marginal cost drops from seconds to ~milliseconds. This is where the hundreds-of-species
advantage comes from.

## Honest scope
- Synthetic genomes are small (~230 kb) — absolute times aren't comparable to the real
  2.5 Mb C. acnes; the **scaling behavior and architecture** are what's demonstrated.
- For *thousands* of strains in one species, the O(n²) pairwise step (even with MinHash)
  would need LSH bucketing; current target (hundreds of species × tens–hundreds of strains
  each) is well within the O(n²·k) regime.
- Still to do for the full claim: a hundreds-of-species × hundreds-of-samples benchmark on
  realistic genomes, and a like-for-like multi-species comparison against StrainScan run
  per-species (needs Linux for StrainScan build).
