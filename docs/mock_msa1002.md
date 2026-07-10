# Real-data validation: ATCC MSA-1002 DNA mock community (native BcgI 2bRAD)

## Data
Figshare article **12272360** v8 (2B-RAD-M unique 2B tags database). `MOCK_MSA1002_2BRAD-M.rar`
(227 MB) = native **BcgI 2bRAD reads (32 bp tags)** of the ATCC MSA-1002 20-strain **even** mock,
at a DNA-input titration **0.001 / 0.01 / 0.1 / 1 / 100 ng**. Downloaded via
`https://ndownloader.figshare.com/files/25933027`; extracted with `bsdtar`. (Shotgun WMS, 16S, and
FFPE files are in the same article.)

## Reference panel
BcgI Strain2bScan DBs for the **20 MSA-1002 species** (13 from the existing 55-species panel; 7
fetched from ENA — *Bifidobacterium adolescentis, Deinococcus radiodurans, Phocaeicola vulgatus,
Lactobacillus gasseri, Bacillus cereus, Clostridium beijerinckii, Cereibacter sphaeroides*) **plus
42 decoy species** (62 total), so both recall and false positives are testable.
`scripts/dl_mock_extra.py` fetches the ENA genomes; DBs are `cluster --enzyme BcgI --similarity 0.95`.

## Result — `results/mock_msa1002_titration.tsv`, `figures/mock_msa1002_titration.*`
Strain2bScan `multi-profile --enzyme BcgI` (gate 50/5) on the real reads:

- **Precision 1.0 at every DNA input** — zero false species strain-resolved.
- **Recall 20/20 at 0.1 and 1 ng; 19/20 at 0.01 ng (56k reads); 6/20 at 0.001 ng (38k reads).**
- Out-of-panel relatives are correctly held in the **detected-not-resolvable** tier (2–15 per run),
  not counted as detections — the three-tier gate working on real open-world data.

This closes the "all-simulated" gap: Strain2bScan reproduces the Strain2bfunc MSA-1002 mock
validation (PPTX slide 26) on **real, error-containing 2bRAD data**, and demonstrates the
**low-biomass niche** (full recall at 0.1 ng, precision 1.0 down to 0.001 ng).

## Refinements still open
- **Cross-species abundance even-ness** (each strain ~5 %): `multi-profile` reports within-species
  cluster abundance; a cross-species abundance estimate + Bray–Curtis vs the even truth is a
  refinement for the final figure.
- **WMS / inStrain comparison** (the PPTX contrast): needs the shotgun files (2.6 GB) + inStrain
  (not installed). PPTX showed inStrain-WMS precision 0.4–0.5 vs 2bRAD 1.0.
- **Host-contamination series** (90/95/99 % human DNA): not in this archive; the PPTX had it.
- **0.001 ng recall** could improve with a lower gate; not pursued (precision was the priority).
