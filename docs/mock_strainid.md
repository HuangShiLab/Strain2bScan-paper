# Fig 12–13 — Strain-level ID + abundance on the ATCC MSA-1002/1003 DNA mock

**Question.** On a real 20-strain DNA mock, can Strain2bScan identify *which strain* of each species is
present (against a database seeded with same-species decoys) and quantify each species' abundance — from
both native BcgI 2bRAD (Part I, Fig 12) and in-silico-digested shotgun (Part II, Fig 13)?

## Database (combined 120-genome tree)
- 20 mock species; per species = the true ATCC genome + **5 conspecific decoy strains** (high-quality:
  CheckM completeness median 98.6 %, contamination ≤ 3.9 %; **within-species ANI > 95 %** verified by skani).
- One combined `strain2bscan cluster` over all 120 genomes → BcgI tree (`--containment`, for native 2bRAD)
  and all-enzyme tree (for shotgun). Truth = each species' ATCC cluster; a correct call requires picking
  the ATCC genome over 5 same-species decoys, and any non-ATCC cluster is a false positive.
- Profile with `--min-abundance 0` (keep every strain; the exposed floor set to 0 for strain-ID) and
  `--min-coverage 0.2` (precision via marker breadth).

## Results (`results/mock_strainid.tsv`)
- **Fig 12 native 2bRAD:** 19/20 strains at 90–99 % host (100 ng), 11/20 at 99.9 %; titration 0/16/19/19
  at 0.001/0.01/0.1/1 ng; ~0 false positives. The one miss, *L. gasseri*, has < 10 unique BcgI tags in the
  combined tree (single-enzyme limit, Fig 5) — still species-detected. `--containment` (Fig 4) is required:
  the > 95 %-ANI decoys are near-clonal and fragment the sparse BcgI tree, stripping unique markers.
- **Fig 13 shotgun:** 20/20 strains, 0 false positives, on both even and staggered mixes; the richer
  all-enzyme marker set avoids sparse-tree fragmentation.
- **Abundance** (both): per-species marker `depth` (∝ genome copy number; new output column) normalised
  across species recovers the mock design — uniform for even MSA-1002, ~3-orders-of-magnitude staggered for
  MSA-1003 (`results/mock_species_abundance.tsv`).

## Reproduce
Raw data, DBs and drivers under `../Strain2bScan-raw-data/MSA_demo/` (see its `db_strategy.md`):
genome re-selection, `strain2bscan cluster` combined trees, profiling with `--min-abundance 0
--min-coverage 0.2`. Figure: `scripts/plot_mock_strainid.py`.
