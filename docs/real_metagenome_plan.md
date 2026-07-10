# Real metagenome validation (saliva) — do we need it, and how

## Do we need it? **Yes — it's the single most valuable addition.**

Every benchmark so far is **simulated, closed-world, error-free**: reads are generated from the
same genomes that build the database, so the truth strain is guaranteed present and there is no
sequencing error. That is exactly the setup a reviewer discounts. A real saliva metagenome adds
the three things simulation cannot:

1. **Real sequencing error** (the sparse 2bRAD marker set's error tolerance is untested);
2. **Open-world** — strains/species *not* in the DB (real false-positive control, which
   closed-world can't measure);
3. **Real community complexity and abundance skew**, plus a concrete demonstration of practical
   utility and of the speed/scale advantage on real data.

The oral/saliva microbiome is a good target: it is well characterized, moderately diverse, and
we already have DBs for many oral species in the 55-species panel (*Streptococcus mutans/
pneumoniae/pyogenes/agalactiae*, *Fusobacterium nucleatum*, *Porphyromonas gingivalis*,
*Actinomyces odontolyticus*, *Treponema denticola*, *Prevotella intermedia*, *Capnocytophaga
sputigena*, *Veillonella parvula*, *Haemophilus influenzae*, *Neisseria meningitidis*, …).

## The core challenge: no ground-truth strain labels

Real samples have no known strain composition, so accuracy can't be scored directly. Standard
ways to validate strain calls on real data (use one or more):

- **A. Concordance with StrainScan** (primary). Run both tools on the *same* saliva reads for the
  species where StrainScan has a database; show they call the same strains/clusters. Agreement
  with the reference method validates without ground truth. *Caveat:* StrainScan DB build needs
  Linux/`dashing`, and its 8 pre-built DBs are only partly oral (*P. copri* overlaps; most don't),
  so direct concordance is limited unless we build oral StrainScan DBs on Linux.
- **B. Species-level sanity** vs an established species profiler (Fast2bRAD-M species layer,
  MetaPhlAn, or Kraken2): the species Strain2bScan strain-profiles should match the species
  actually present.
- **C. Reproducibility** (if ≥2 samples per subject or technical replicates): strain calls should
  be consistent within subject/replicate and differ across subjects — a real, quantifiable signal
  with no ground truth needed.
- **D. Known markers / isolates** (if available): any cultured oral isolate or prior strain typing
  from the same subjects gives partial ground truth.

## The actual dataset (confirmed) — a strong, self-validating design

**Paired shotgun + BcgI-2bRAD saliva, 8 subjects × 4 timepoints/day (9AM/11AM/1PM/5PM) = 32
samples per data type (64 total).** Biological question: does **strain-level** diversity
distinguish individuals better than species-level? Prior lab result with **strain2bfunc**
(github.com/HuangShiLab/strain2bfunc): yes, for 7 species (higher R² / PERMANOVA at strain level).
strain2bfunc needs a 9.3 GB DB, ~30 GB RAM, ~40 min/sample — Strain2bScan is ~78 MB / ~1 s, so we
can **reproduce the finding at a fraction of the compute.**

Having **paired shotgun + 2bRAD of the same samples** is what makes this self-validating (no
external ground truth needed).

## Experiment plan

1. **Species selection.** Species-level profiling (2bRAD-M / Fast2bRAD-M) → the dominant oral
   species. Build Strain2bScan DBs for them (many already in the 55-panel; extend to match the
   strain2bfunc species set for a like-for-like comparison; pin accessions).
2. **Strain profiling.** `multi-profile` all 64 samples (species gate on) → strain (cluster)
   abundance matrix, separately for shotgun and 2bRAD; record runtime/memory.
3. **Individual discrimination (headline).** Bray–Curtis on strain-level abundances → PCoA +
   **PERMANOVA with subject as factor (R²)**; compare to the same at **species level**. Expect
   strain-level R² > species-level (individuals more separable at strain level) — replicating
   strain2bfunc, now with a fast/light tool.
4. **Paired shotgun-vs-2bRAD concordance (validation #1).** Per sample, compare strain calls from
   shotgun (in-silico digest) vs real 2bRAD (BcgI). High concordance validates **both** the tool
   on real error-containing reads and the in-silico-digestion↔real-2bRAD equivalence — without any
   external ground truth.
5. **Temporal stability (validation #2).** Within-subject across the 4 timepoints: strain
   composition should be stable within-day and distinct between subjects (within/between variance,
   ICC).
6. **Efficiency on real data (headline #2).** 64 real samples' runtime/memory vs strain2bfunc
   (~40 min / 30 GB) and StrainScan where a DB exists — the speed/scale advantage, now on real data.

This one chapter closes the closed-world/error-free gap (real reads), gives open-world FP control
(real community), demonstrates the native-2bRAD capability (real 2bRAD), reproduces a known
biological result cheaply, and is internally validated by the paired design.

## What's needed to start

- **Data access:** FASTQ paths for the 64 samples. 64 saliva metagenomes may be large (tens–hundreds
  of GB) — can a subset (or all) be shared to this machine, or must this run on an HPC?
- **The species set:** the 7 species from the strain2bfunc result (for a like-for-like comparison),
  or let us select the dominant oral species from the data.
- **Sample sheet:** subject ID + timepoint + data-type per FASTQ (for PERMANOVA grouping).
