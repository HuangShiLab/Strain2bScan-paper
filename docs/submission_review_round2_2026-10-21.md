# Round-2 isolated review and revision record (2026-10-21)

## Review packet

The immutable packet was `manuscript/full_manuscript.md`, the numbered figures, table draft and
checksummed result tables. Three mutually blind reviewer contexts were run. Emphasis was, respectively:
significance and journal fit; benchmark validity and provenance; submission structure and readability.

## Blocking concerns and resolutions

| ID | Concern | Resolution |
|---|---|---|
| R1-M1 | High-host robustness was stated too generally. | All native host-ladder claims now specify MSA-1002 and distinguish the clean other mocks. |
| R1-M2 | The selected three-prefix saliva result conflicted with the full-depth rebuild. | The prefix result is labelled exploratory and non-depth-matched. A new gate-matched 32-pair native rerun and summary are now primary for the full rebuild. |
| R1-M3 | “Only tested tool” excluded StrainGST without an abundance-error comparison. | Removed. The claim now describes the primary shotgun comparators and explicitly separates tool-specific call units. |
| R2-M1 | Comparator reference spaces were conflated. | Methods, Results, Discussion and captions now define 0.95 clusters, StrainGST 0.90 references and inStrain genome references. Accuracy is descriptive and not pooled. |
| R2-M2 | The saliva comparison mixed R1 prefixes, read caps and thresholds. | The prefix comparison is bounded as exploratory. All 32 native and WMS libraries were reprofiled with the same species/support/coverage gates and are reported in Table S4. |
| R2-M3 | F1/provenance and Table 2 inconsistencies. | Corrected F1 values; froze the threshold-aware mock metric table; repaired all 1,125 enzyme-sweep prediction paths and added SHA-256 values plus a prediction archive. |
| R3-M1 | Tables 2–3 were malformed. | Rebuilt with explicit accuracy/cost columns, source-specific n, units, environments and separate StrainGST columns. |
| R3-M9 | Supplementary numbering was incomplete. | Supplementary tables renumbered S1–S4 and cited in Results. Figure S6 internal labels now match manuscript numbering. |

## Other revisions

- Added the available subject-aware saliva confidence intervals and the zero-inclusive strain-minus-species interval.
- Added species-level failures to the n=2025 summary, especially *Klebsiella pneumoniae* and *Salmonella enterica*.
- Qualified database-build cost as an archived cross-platform comparison and separated same-container profiling cost.
- Added the multi-species accuracy trade-off to the abstract.
- Corrected custom PERMANOVA provenance and restricted scaling language to the measured range.
- Added panel-level captions, rank-rank scatter wording, StrainGST nomenclature, DOIs, funding-role text and Table S4.
- Renamed final Figure 5 and 6 assets consistently with manuscript numbering.

## Current release gates

```bash
python3 scripts/assemble_manuscript.py
python3 scripts/check_manuscript_claims.py
python3 scripts/make_figure_manifest.py
git diff --check
python3 -m compileall -q -f scripts
```

All commands pass as of this record. The current DOCX is
`manuscript/Strain2bScan_manuscript_2026-10-21-submission.docx` and renders to 49 pages.

## Remaining author-side items

1. Supply the grant funder/programme name if required by the chosen journal.
2. Add ORCID IDs and journal administrative forms.
3. Mint a Zenodo/archival DOI for the final source and derived-output snapshot.
