# Can StrainScan build a database on macOS? A corrected, rigorously re-checked answer

**Short answer: no, in practice — but for a more precise reason than originally stated.**
StrainScan's build step (`StrainScan_build.py`) unconditionally calls `dashing` (even when a
custom clustering file `-c` is supplied — the code runs `Cluster.construct_matrix()`, which
invokes dashing, *before* checking for `-c`), and the only available macOS build of `dashing`
does not functionally work. Profiling (`StrainScan.py`) does **not** need dashing and works
fine on macOS once `jellyfish` is fixed (see `docs/headtohead_strainscan.md`).

## What was checked, first-hand, on this machine (16-core Apple Silicon, macOS 14.7)

1. **`uname -m` → `arm64`.**
2. **The bundled binaries are genuine Linux ELF, not just "wrong CPU arch."** `file` on
   `StrainScan/library/dashing_s128` reports `ELF 64-bit LSB executable, x86-64, ... GNU/Linux`.
   This is an OS-boundary mismatch: Rosetta 2 translates macOS x86_64 Mach-O binaries, not
   binaries built for a different kernel/ABI. No amount of Rosetta configuration fixes this.
3. **Correction to an earlier claim.** We had stated "dashing has no macOS build." This was
   **wrong** — it generalized from a failed install of the *full* `strainscan` bioconda
   meta-package (blocked by an unrelated dependency, `sibeliaz`, which has no osx-64 build) to
   `dashing` itself, without testing `dashing` in isolation.
4. **Retested `dashing` alone:** `CONDA_SUBDIR=osx-64 conda create -n dashing_test -c bioconda
   -c conda-forge dashing` **succeeds** — `dashing-0.4.0-hceddf3b_2` installs, and `file` on the
   resulting binary confirms a genuine `Mach-O 64-bit executable x86_64` (not mislabeled).
   A native arm64 build does not exist (default-subdir `conda search dashing` returns nothing
   matching, only an unrelated `python-dashing` TUI package).
5. **Actually running it is where it fails.** Invoked via StrainScan's exact call pattern
   (`dashing dist -p{N} -k31 -O distance_matrix.txt -o size_estimates.txt -Q path_file.txt -F
   path_file.txt`) against two small real genomes (~2.5 Mb each — should complete in well under
   a second): the process **hangs indefinitely** — confirmed via `ps` showing `0.0% CPU`,
   uninterruptible-wait state, for minutes, with **no output, no crash, no error**. Reproduced
   with both `-p4` (multi-threaded/OpenMP) and `-p1` (single-threaded) — ruling out an
   OpenMP-under-Rosetta-specific deadlock. `lsof` on the stuck process shows Rosetta's
   ahead-of-time translation cache (`/private/var/db/oah/.../*.aot`) actively populated for
   `libomp`, `libc++`, `libz`, and `dyld` — so the binary *did* load and *is* being translated;
   it simply never produces output or exits.

## Bottom line

| claim | status |
|---|---|
| "No macOS `dashing` package exists" | **Wrong** (corrected here) — an osx-64 package exists and installs |
| "StrainScan build works on macOS if you just install dashing" | **Also wrong** — it installs but hangs when run, under Rosetta 2, on this hardware |
| "StrainScan build needs Linux (or a real Linux environment) on this machine" | **Confirmed**, by direct execution, not by absence of a package |

## A viable next step, not yet pursued

Docker (29.2.1) is installed on this machine (daemon not running by default). A Linux
**container** — unlike Rosetta 2, which translates x86_64 macOS binaries — runs a real Linux
kernel (natively for arm64 Linux images via Apple's virtualization framework, or emulated via
QEMU for x86_64 Linux images). Either should let the bundled `dashing_s128` **Linux ELF**
binary execute correctly, sidestepping the Rosetta hang entirely, and would let the same-panel
StrainScan build-and-profile comparison (`scripts/run_headtohead_strainscan.py`) run locally
instead of on a separate Linux/HPC machine. Not attempted here — starting Docker Desktop and
pulling images is a larger, separate undertaking; flagged for a future session if a local
same-panel build comparison becomes a priority.
