# strain2bRAD — manuscript framing & benchmark plan

Assessment of the three intended claims and the experiments needed to support each, with
evidence gathered so far. Bottom line: **the story is sound and publishable as a methods/
application paper** ("a fast, lightweight Rust strain profiler using 2bRAD-reduced markers,
honest about where strain resolution is feasible"), provided the claims are framed as a
speed/memory win with a controlled sensitivity trade-off — not as an accuracy win.

---

## Claim 1 — Performance (runtime + memory) ≫ Python StrainScan at scale

**Does it make sense?** Yes, with one caveat to control. The speedup has **two confounded
sources**: (a) Rust vs Python, and (b) 2bRAD markers (~1–2% of k-mers) vs the full k-mer set.
A reviewer will ask which is responsible. Design the benchmark to separate them.

**Experiments**
- Primary (what users care about): strain2bRAD (Rust+2b) vs StrainScan (Python+full k-mer),
  identical genomes & samples. Report **wall-clock + peak RSS** for (i) DB build and (ii)
  per-sample profile, across **panel sizes** (10/50/100/500/all genomes) and **read depths**.
- Attribution: report the **marker-count reduction** (tags vs k-mers per genome) so the
  data-size win is explicit; ideally add a dense-k-mer Rust variant to isolate language vs
  marker effects. At minimum, state that both contribute.
- Threads: report single- and multi-threaded; StrainScan parallelizes DB build, so strain2bRAD
  must too for a fair large-scale claim.

**Evidence so far** (this machine, 16-core arm64; 64-genome C. acnes panel, 14 enzymes):
- DB build: **24.3 s, 159 MB peak RSS**, 60 clusters, 34.5 MB DB.
- Per-sample profile (200k reads): **3.6 s, 106 MB peak RSS**.
- StrainScan numbers: pending (see head-to-head section).

**Risk / gap**: "large scale" requires scaling work first — the build is **single-threaded**
and clustering is **O(n²·m)** (fine at 64 genomes, not at 500+). Add `rayon` + inverted-index
single-linkage before the headline large-panel benchmark, or the claim is vulnerable.

---

## Claim 2 — 2b-based k-mers give fast strain ID vs StrainScan's full k-mer set

**Does it make sense?** Yes as **"comparable accuracy at sufficient depth, much faster,"**
NOT as "more accurate." 2b markers are ~50–100× sparser, which costs sensitivity for
low-abundance / low-depth strains and for many co-present near-relatives.

**Experiments**
- **Accuracy vs depth sweep** (1×/5×/10×/30×/100×): find the depth where 2b matches full
  k-mer — this empirically defines "sufficient depth," the crux of the paper.
- **Accuracy vs mixture complexity** (1…N co-present strains; even vs skewed abundance).
- **Head-to-head with StrainScan** on identical samples (accuracy + resources together).
- **Multiple species** (C. acnes, S. aureus, S. epidermidis — all in your mock repo).
- Metrics: precision/recall/F1 (presence) + abundance error (L1, Bray–Curtis, JSD) at **both
  strain and cluster resolution** (report resolution = #clusters/#genomes alongside accuracy;
  otherwise coarse clustering inflates scores — see the BcgI artifact below).

**Evidence so far** (C. acnes, 64-genome panel, ~12× total, 14 enzymes, after the #1 fixes):
- Precision **0.90**, recall **0.56** over 16 truth strains in 5 samples.
- Abundance accurate for dominant strains (Bray–Curtis 0.03) and good for the 5-strain
  mixture after the unique-marker-depth estimator (Bray–Curtis 0.44 → **0.16**).
- Recall is capped by low-abundance minors (2–8% → ~1–4× per-strain depth) below the sparse
  marker sensitivity — exactly the trade-off to characterize with the depth sweep.

---

## Claim 3 — BcgI 2bRAD experimental data for strain diversity (not all species)

**Does it make sense?** Yes, as an **exploratory capability**, and the "not necessarily for
all species" framing is scientifically correct and a strength (honest reporting).

**Experiments**
- **Resolvability survey**: across many species, count cluster/strain-specific BcgI tags per
  cluster; report the fraction of species/clusters that clear a resolvability threshold.
- Demonstrate success on a species with sufficient BcgI markers, and the honest "not doable"
  output on one without (e.g., C. acnes).
- If available, run on **real BcgI 2bRAD reads** end-to-end.

**Evidence so far**: BcgI on C. acnes yields only ~16 coarse clusters from 64 genomes
(strains lump together; sample4's 5 strains all collapse into one cluster), so it is **not
strain-resolvable** for C. acnes — the tool reports this. (Caution: coarse clustering makes
cluster-level scores look perfect; always report resolution.)

---

## Positioning & gaps to publication
- **Niche**: fastest/lightest profiler via 2bRAD-reduced markers + a unique capability to
  ingest BcgI 2bRAD experimental data. Compare against StrainScan (primary), and cite
  Strainify / StrainGE / sylph for context.
- **Must-do before submission**: (1) parallelize + scale (rayon, inverted-index clustering)
  for the large-panel performance claim; (2) depth & complexity sweeps; (3) head-to-head vs
  StrainScan across ≥3 species; (4) proper Layer-1 CST and depth-aware detection to lift
  recall on minors. Items (1) and (4) are the highest-leverage.
- **Target venues**: Bioinformatics (Application Note if scope is tight), GigaScience,
  Microbiome, or mSystems.
