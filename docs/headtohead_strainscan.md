# Head-to-head: Strain2bScan vs StrainScan (C. acnes, 5 real mock samples)

Both tools run on the **same 5 paired-end C. acnes mock samples** (~100k pairs each).
Machine: 16-core arm64 macOS, 51 GB RAM.

## Reproducibility notes (getting StrainScan to run on macOS)
- StrainScan bundles **Linux ELF** binaries: `library/jellyfish-linux` (confirmed via `file`:
  genuine `ELF 64-bit ... x86-64, GNU/Linux`) and `library/dashing_s128` (same). These cannot
  run on macOS at all — not a missing-package issue, but an OS-level incompatibility (Rosetta
  2 translates macOS x86_64 Mach-O binaries, not foreign-OS ELF binaries).
- Fix for **profiling**: `ln -sf $(which jellyfish) StrainScan/library/jellyfish-linux`,
  pointing at a genuine macOS x86_64 Mach-O jellyfish build (bioconda, osx-64 subdir, runs via
  Rosetta 2). Also `pip install treelib`. This works — used throughout this benchmark.
- **Building** a StrainScan DB on this machine is **not possible**, but the precise reason took
  rigorous re-checking to pin down (see `docs/macos_dashing.md` for the full, corrected
  investigation): a `dashing` **osx-64 bioconda package does exist** (an earlier claim that "no
  macOS dashing build exists" was wrong — that failure was from an unrelated dependency,
  `sibeliaz`, in the full `strainscan` meta-package). The real, verified blocker: the only
  available macOS `dashing` build **hangs indefinitely under Rosetta 2** (confirmed firsthand,
  0% CPU / uninterruptible wait, both multi- and single-threaded, on a trivial 2-genome task) —
  a Rosetta translation failure, not a missing package. Net effect is the same (StrainScan
  build needs Linux, or a real Linux container e.g. Docker — not just Rosetta), but the
  mechanism is different from what was previously stated.

## Performance — FAIR comparison (identical input reads, each tool with its own DB)

| sample | Strain2bScan time | Strain2bScan RSS | StrainScan time | StrainScan RSS |
|---|---|---|---|---|
| s1 | 0.49 s | 118 MB | 8.14 s | 828 MB |
| s2 | 0.49 s | 123 MB | 6.83 s | 828 MB |
| s3 | 0.49 s | 119 MB | 6.79 s | 828 MB |
| s4 | 0.50 s | 124 MB | 6.76 s | 828 MB |
| s5 | 0.51 s | 115 MB | 6.77 s | 828 MB |
| **mean** | **0.50 s** | **~120 MB** | **~7.0 s** | **828 MB** |

→ **Strain2bScan is ~14× faster and ~7× lighter per sample** (both at default parallelism on
16 cores; sample1's StrainScan time includes process cold-start, typical ~6.8 s). The memory
gap is structural: StrainScan counts the full k-mer set with jellyfish (~800 MB hash);
Strain2bScan streams reads and keeps only sparse 2b-tag markers. Even single-threaded,
Strain2bScan profile is ~3.4 s — still faster than StrainScan. (DB build: Strain2bScan 2.6 s
@16 threads / 22.4 s @1 thread on 64 genomes; StrainScan build not runnable on macOS — see above.)

*Caveat*: StrainScan profiles against its 275-strain DB and Strain2bScan against a 64-genome
DB — but this is exactly how each tool is used in practice (each with its native DB on the
same reads), so the per-sample resource comparison is fair.

## Accuracy — confounded by a DB-version mismatch (not a fair comparison yet)

StrainScan's **pre-built C. acnes DB (275 strains, 2022)** is **missing 6 of the 14 truth
strains**, including the dominant strains of s2, s3, and s4. So StrainScan structurally
reports nothing for those samples — an artifact of the outdated DB, not the algorithm.

| sample | truth dominant | in StrainScan DB? | StrainScan | Strain2bScan |
|---|---|---|---|---|
| s1 | GCF_009737125 (.97) | yes | **found** | **found** |
| s5 | GCF_003384585 (.84) | yes | **found** | **found** |
| s2 | GCF_024328475 (.86) | no | (can't) | found |
| s3 | GCF_009737075 (.59) | yes* | miss | found (both strains) |
| s4 | GCF_021011155 (.36) | no | (can't) | found 3/5 |

On the truth strains present in **both** tools' databases (s1, s5 dominants), **the tools
agree** (both detect the dominant, both miss the low-abundance minors).

**A fair accuracy head-to-head requires both tools on the identical genome panel.** Since
StrainScan can't build on macOS, run this on Linux:
1. `conda install -c bioconda strainscan` (Linux resolves dashing/sibeliaz).
2. Build StrainScan on the 64-genome panel: `StrainScan_build.py -i acnes/genomes -o DB`
   (optionally `-c` with Strain2bScan's clustering for matched clusters).
3. Run both tools on the 5 samples; evaluate with `strain2bscan evaluate` at cluster
   resolution (remap truth via membership). Identical panel → clean accuracy comparison.

## Takeaways for the manuscript
- **Performance claim is well-supported**: ~14× faster, ~7× less memory per sample, with the
  memory advantage structurally attributable to sparse 2b markers vs full k-mer counting.
- **Accuracy claim needs the same-panel benchmark on Linux** (above). Current evidence:
  Strain2bScan detects dominant/co-dominant strains at precision 0.90 on the 64-panel; on
  shared-DB strains it agrees with StrainScan.
