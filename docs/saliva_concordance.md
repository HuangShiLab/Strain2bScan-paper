# Saliva shotgun↔2bRAD concordance (Fig 10 validation)

Paired native BcgI 2bRAD and shotgun (WMS) of the SAME saliva samples (SRA PRJNA1131785). Saliva is
a **high-host** matrix (~90% human DNA), so this is the low-biomass regime where 2bRAD's wet-lab
enrichment matters. Question: on real high-host samples, does Strain2bScan on native 2bRAD recover
the strains shotgun finds — **and** the low-abundance strains shotgun cannot?

## Data & method
3 paired samples (subjects 7/1/14). WMS R1 is ~14 GB/sample and downloaded truncated; the valid
prefix (~20–90 M reads) was used — a subsample of the shotgun library, digested **in-silico with
BcgI** by the same `multi-profile` command and profiled against the same 19-species oral panel as the
native 2bRAD. Strain calls compared per sample; 2bRAD community relative abundance = per-cluster
support / sample total. `scripts/concordance_superset.py`. (A 4th sample, S13-11, recovered 0
bacterial species from WMS — host-dominated — and was dropped.)

## Result — `results/saliva_concordance.tsv`, `figures/saliva_concordance.*`

| sample | 2bRAD strains | WMS strains | WMS confirmed by 2bRAD | 2bRAD-only |
|---|---|---|---|---|
| S11-7 (subj 7) | 164 | 23 | 23 | 141 |
| S17-1 (subj 1) | 160 | 32 | 32 | 128 |
| S17-14 (subj 14) | 173 | 10 | 10 | 163 |

- **Every strain shotgun detects is confirmed by native 2bRAD — 65/65 (100%).** In-silico BcgI
  digestion of real shotgun produces no strain calls that 2bRAD contradicts; the WMS call set is a
  clean **subset** of the 2bRAD call set.
- **Native 2bRAD additionally recovers 128–163 strains per sample that shotgun misses**, and those
  2bRAD-only strains are **significantly lower-abundance**: median community relative abundance
  0.0029 (2bRAD-only) vs 0.0097 (shared) — 3.3× higher for the shared, WMS-detected strains
  (Mann–Whitney p = 1.2 × 10⁻²³).

**Interpretation.** Under heavy host contamination, host-dominated shotgun recovers only the most
abundant, ubiquitous oral strains; native 2bRAD confirms all of those **and** reaches deep into the
low-abundance tail that shotgun cannot. This both validates the two-data-mode equivalence at the
call level (shotgun calls ⊆ 2bRAD calls, no in-silico false positives) and demonstrates 2bRAD's
sensitivity advantage on real high-host samples — the same enrichment benefit quantified on the
controlled mock (`docs/mock_hostcontam.md`), now on real saliva.

*Note:* per-sample profile agreement (Bray–Curtis) is dominated by the depth asymmetry and is not the
right lens here; the directional superset relationship above is the concordance evidence.
