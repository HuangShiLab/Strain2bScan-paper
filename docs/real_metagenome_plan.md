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

## Two scenarios — which data is it?

**If the saliva data is shotgun metagenome** (most likely):
- Digitally digest the reads (all-16 or the ~8-enzyme set) and strain-profile the oral species.
- Demonstrates: works on real error-containing reads + real complexity; measure the speed/memory
  advantage on real data; validate by **concordance (A) + species sanity (B) + reproducibility (C)**.

**If the saliva data is BcgI 2bRAD (a wet-lab 2bRAD library)** — the *stronger* story:
- This directly demonstrates the **unique capability** no other strain tool has: strain profiling
  from real 2bRAD experimental data (paper claim 3). Reads are the BcgI tags; profile with the
  BcgI-only DB (verified path, `docs/benchmark_datasets.md` #9).
- Demonstrates the *low-input* economics of real 2bRAD (deep per-tag coverage at low sequencing
  cost) that in-silico shotgun digestion cannot show.

## Concrete plan (once the data type is known)

1. Assemble an **oral-species DB set** (extend the 55-panel with the dominant saliva taxa; pin
   accessions).
2. Run Strain2bScan `multi-profile` on each saliva sample (species gate on) → per-sample strain
   calls + abundances + wall-clock/memory.
3. **Validate:** concordance with StrainScan where possible (A), species-level agreement (B), and
   within-subject/replicate reproducibility (C).
4. Report: strains detected per sample, inter-sample/inter-subject variation (biologically
   interpretable), real-data runtime/memory, and — if 2bRAD data — the native-2bRAD demonstration.

## What we need from you to start

- **Data type:** shotgun metagenome, or BcgI 2bRAD library?
- **Scope:** how many samples; any subjects / timepoints / technical replicates?
- **Any orthogonal truth:** prior species/strain profiling, cultured isolates, or metadata?
- **Location/format:** FASTQ paths (and whether it can be shared to this machine, or must run on an HPC).
