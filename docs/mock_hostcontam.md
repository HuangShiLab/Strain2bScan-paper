# DNA mock, host-contamination series — native 2bRAD vs shotgun (Fig 9 extension)

Extends the ATCC MSA-1002 mock validation (`docs/mock_msa1002.md`, Figshare titration) with the
**host-contamination axis** and the **shotgun (WMS) contrast**, both requested to close the "both
shotgun and 2bRAD, high-host" gap. Data: SRA **PRJNA1131785** (npj Biofilms Microbiomes 2025), the
`synthetic metagenome` samples — a designed spike of the 20-species even mock into human DNA.

## Data
Seven MSA-1002 DNA-mock runs from PRJNA1131785 (`source_material_id` decodes the design):

| data type | host DNA | runs |
|---|---|---|
| native BcgI 2bRAD | 90% | SRR29727660 (rep1), SRR29727659 (rep2) |
| native BcgI 2bRAD | 99% | SRR29727658 (rep1), SRR29727657 (rep2) |
| shotgun WMS (DNBSEQ) | 0% (control) | SRR29727654 |
| shotgun WMS (DNBSEQ) | 90% | SRR29727656 |
| shotgun WMS (DNBSEQ) | 99% | SRR29727655 |

Reference: the reconstructed **62-species BcgI DB** (55-panel + 7 extra = 20 MSA-1002 species + 42
decoys; `scripts/dl_genomes.py`, then `cluster --enzyme BcgI --similarity 0.95`). Profiled with
`multi-profile --enzyme BcgI --min-species-markers 50 --min-species-detect 5`; shotgun WMS is
digested **in-silico with BcgI** by the same command — a like-for-like comparison against native 2bRAD
on the identical DB. Species scored strain-resolved vs the 20-species truth.

## Result — `results/mock_hostcontam.tsv`, `figures/mock_hostcontam.*`

| data type | host % | sample markers | precision | recall | detected-not-resolved |
|---|---|---|---|---|---|
| 2bRAD | 90 | ~365k | **1.00** | **20/20** | 2 |
| 2bRAD | 99 | ~337k | **1.00** | **20/20** | 2 |
| WMS | 0 | 96k | 1.00 | 20/20 | 2 |
| WMS | 90 | 53k | 1.00 | 20/20 | 2 |
| WMS | 99 | 33k | 1.00 | **12/20** | 9 |

**Two findings:**
1. **Precision = 1.00 at every contamination level, for both data types** — zero false-positive strain
   calls even when 99% of the DNA is human. Out-of-panel relatives correctly stay in the
   detected-not-resolvable tier (the three-tier gate working on real, host-dominated data).
2. **Native 2bRAD is markedly more robust to host contamination than shotgun.** At 99% human DNA,
   2bRAD retains **full 20/20 recall** while in-silico-digested shotgun drops to **12/20**. The
   mechanism is visible in the marker yield: 2bRAD delivers ~340–370k usable BcgI markers regardless
   of host fraction (the reduction happens at the wet-lab enzyme step, before host swamps the
   library), whereas digesting an already host-dominated shotgun library in silico yields only
   96k → 53k → 33k markers as host rises 0 → 90 → 99%. Fewer surviving markers → species fall below
   the resolve gate.

This reproduces the prior Strain2bfunc claim (2bRAD precision/recall ~1.0 at ≥90% human) with
Strain2bScan on real reads, quantifies it against shotgun on the same samples/DB, and gives the
low-biomass / high-host-contamination selling point a direct, mechanistic figure. (inStrain was not
run — the contrast here is native-2bRAD vs in-silico-digested-shotgun within one tool; an inStrain-WMS
comparison remains optional future work.)
