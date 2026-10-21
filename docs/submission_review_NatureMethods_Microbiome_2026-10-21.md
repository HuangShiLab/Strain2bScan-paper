# Pre-submission review and resolution record

**Targets considered:** Nature Methods; Microbiome.
**Paper type:** method/tool paper with benchmarking and microbiome applications.
**Manuscript reviewed:** `manuscript/full_manuscript.md` after Results simplification and new benchmark reruns.

## Executive decision

The manuscript is close to submission-ready for **Microbiome** and is also a defensible methods/benchmark
submission for **Nature Methods**, provided the remaining administrative and running-analysis items below are
completed. The core claim is now scoped correctly: Strain2bScan extends reduced-representation sequencing to
strain-cluster resolution, with the strongest advantages in native low-biomass/high-host profiling, database
cost, multi-species shotgun runtime and 99%-host retention. It does not claim uniform superiority across all
clean mocks or all metrics.

## Issues resolved in this revision

### CRITICAL 1 — inStrain simulation comparison was incomplete
**Problem:** Earlier text included inStrain in Figure 8 but the repo only contained exploratory simulation
results.
**Resolution:** Completed all 225 matched single-species simulations with inStrain 1.11.0. Outputs:
`results/instrain_sim225/persample.tsv`, `by_depth.tsv`, `by_species.tsv`, `provenance.json`.
**Status:** Fixed.

### CRITICAL 2 — Figure 3 lacked per-sample provenance
**Problem:** The enzyme-sweep table contained only 70 aggregate rows and did not state the samples behind
each aggregate.
**Resolution:** Reran the complete enzyme ladder on all 225 matched simulations. The new table contains
1,125 rows (225 samples × five enzyme levels), plus database marker/cluster provenance.
**Status:** Fixed.

### CRITICAL 3 — paired saliva WMS data were incomplete
**Problem:** The published comparison used four prefix libraries although the BioProject contains 32 paired
WMS samples.
**Resolution:** Rebuilt all 32 samples into canonical gzipped mate files on HKU HPC. The 64 files passed
`gzip -t` (64 OK, 0 FAIL). Manifests and conversion evidence are in `results/wms32_rebuild/`.
**Status:** Fixed.

### MAJOR 1 — n=2025 result was historical
**Problem:** Figure 7 previously used archived results that could not be tied to the current software
version.
**Resolution:** Reran all 2,025 simulations with current-source Strain2bScan. Figure 7 and
`results/sim2025_current*` now use the current-source run.
**Status:** Fixed.

### MAJOR 2 — comparator reference spaces were easy to misread
**Problem:** Strain2bScan/StrainScan, StrainGST and inStrain do not emit calls in the same cluster space.
**Resolution:** Results, Figure 8 legend, Table 8 and Methods now state each tool's reference/cluster space
and explicitly say accuracies are not pooled.
**Status:** Fixed.

### MAJOR 3 — missing submission declarations
**Problem:** Competing interests, ethics/data-use statement, acknowledgements and benchmark software
versions were incomplete.
**Resolution:** Added Acknowledgements, Competing interests, Ethics, and explicit source-commit/tag
provenance for simulation and ATCC-mock benchmarks.
**Status:** Fixed.

### MINOR 1 — Results readability
**Resolution:** Results reduced from approximately 3,000 to approximately 1,250 words. Each benchmark now
states dataset → primary observation → bounded conclusion.
**Status:** Fixed.

### MINOR 2 — stale figure references
**Resolution:** Removed obsolete Figure 11 references and synchronized current Figure 8/10 numbering.
**Status:** Fixed.

## Remaining items before final submission

1. **Canonical 32-sample WMS profiling job:** Strain2bScan is currently running on HKU HPC over all 32
   complete paired WMS samples. This is not required for the existing prefix-subsample claim, but completing
   it will preempt a reviewer question about the full dataset.
2. **Archival code/data DOIs:** GitHub tags exist, but Zenodo or another DOI archive should be created for
   the exact manuscript, source-code and derived-output snapshots.
3. **Author-side metadata:** ORCID iDs and any journal-specific funding/declarations forms still need to be
   entered by the authors.
4. **Journal choice:** Microbiome is the strongest broad microbiome-method fit. Nature Methods is plausible
   if the submission emphasizes the benchmarking and algorithmic scaling advance; the abstract is already
   short enough to adapt easily.
