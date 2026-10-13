# Submission checklist

## Completed

- [x] Authors, affiliations, corresponding author and email.
- [x] Funding award number.
- [x] Data and code availability with versioned repository/tag.
- [x] Figures 1–11 and Supplementary Figures S1–S5.
- [x] Tables 1–9 and Supplementary Tables S3–S5.
- [x] Methods for current Strain2bScan, StrainScan, inStrain and StrainGST comparators.
- [x] Failure-aware and database-resolution benchmark sensitivity.
- [x] Frozen artifact checksums and manuscript claim checks.
- [x] Rust formatting, clippy and test gates for the current benchmark snapshot.
- [x] Readable DOCX with embedded figures and tables.

## Author action required before submission

- [ ] Replace `[xxxx]` in the Author contributions section with the approved CRediT/text statement.
- [ ] Add ORCID identifiers if required by the target journal.
- [ ] Add competing-interest, ethics/data-use and acknowledgements statements required by the journal.
- [ ] Confirm the target journal's figure/table formatting and reference-style requirements.

## Recommended but not blocking

- [ ] Create a Zenodo DOI for the `paper-2026-10-12` tag.
- [ ] Run a native Linux/x86_64 StrainScan timing benchmark to remove the emulation caveat.
- [ ] Restore raw inStrain genome tables and run call-filter sensitivity.
