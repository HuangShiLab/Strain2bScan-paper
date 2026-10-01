# Strain2bScan real-data subset run

## Run scope

- 22 paired-end WGS samples from PRJNA288562, PRJNA1517970, PRJNA1191223, and PRJNA1191225.
- 20-species legacy MSA panel; enzyme mode `all`; gates: `--min-species-markers 20 --min-species-detect 2 --min-support 2 --min-coverage 0.01 --min-abundance 0`.

## Efficiency

- Wall time: **100.59 s** for 22 samples (**4.57 s/sample**).
- Peak RSS: **343 MiB** for the whole batch.
- Tiny vaginal sample single-profile test: **0.44 s**, **191 MiB**.

## Calls by project

| project      |   samples |   positive_samples |   strain_calls |   species_observed |
|:-------------|----------:|-------------------:|---------------:|-------------------:|
| PRJNA1191223 |         3 |                  3 |              9 |                  3 |
| PRJNA1191225 |         6 |                  0 |              0 |                nan |
| PRJNA1517970 |         7 |                  0 |              0 |                nan |
| PRJNA288562  |         6 |                  5 |             34 |                  6 |

## Per-sample calls

| project      | run         | condition                      |   distinct_markers |   strain_calls |   species_calls |   max_depth |   sum_sample_fraction |
|:-------------|:------------|:-------------------------------|-------------------:|---------------:|----------------:|------------:|----------------------:|
| PRJNA1191225 | SRR32076385 | Escherichia coli LHCA45        |             219245 |              0 |               0 |      0      |              0        |
| PRJNA1191225 | SRR32076393 | Escherichia coli LHCA56        |             242121 |              0 |               0 |      0      |              0        |
| PRJNA1191225 | SRR32076407 | Escherichia coli LHCA72        |             183313 |              0 |               0 |      0      |              0        |
| PRJNA1191225 | SRR32076422 | Bifidobacterium breve LHCA81   |             142707 |              0 |               0 |      0      |              0        |
| PRJNA1191225 | SRR32076429 | Bifidobacterium bifidum LHCA82 |             129681 |              0 |               0 |      0      |              0        |
| PRJNA1191225 | SRR32076444 | Bifidobacterium longum LHCA43  |              65183 |              0 |               0 |      0      |              0        |
| PRJNA1191223 | SRR32121883 | W1                             |             748209 |              1 |               1 |      0.0398 |              6.5e-05  |
| PRJNA1191223 | SRR32121884 | W2                             |             694111 |              2 |               1 |      0.0153 |              0.000118 |
| PRJNA1191223 | SRR32121885 | W3                             |             830986 |              6 |               1 |      0.0616 |              0.000413 |
| PRJNA1517970 | SRR40410528 |                                |             116410 |              0 |               0 |      0      |              0        |
| PRJNA1517970 | SRR40410567 |                                |               9693 |              0 |               0 |      0      |              0        |
| PRJNA1517970 | SRR40410657 |                                |              31171 |              0 |               0 |      0      |              0        |
| PRJNA1517970 | SRR40410734 |                                |              24067 |              0 |               0 |      0      |              0        |
| PRJNA1517970 | SRR40410753 |                                |              14802 |              0 |               0 |      0      |              0        |
| PRJNA1517970 | SRR40410776 |                                |                464 |              0 |               0 |      0      |              0        |
| PRJNA1517970 | SRR40410851 |                                |              23718 |              0 |               0 |      0      |              0        |
| PRJNA288562  | SRR6747958  | vagina_GD273                   |             226983 |              1 |               1 |      0.0101 |              0.000123 |
| PRJNA288562  | SRR6748028  | gut_GD84                       |            1918932 |              7 |               2 |      0.0499 |              0.000838 |
| PRJNA288562  | SRR6748033  | gut_GD273                      |            1595769 |              2 |               1 |      0.0518 |              0.00019  |
| PRJNA288562  | SRR6748088  | saliva_GD273                   |             876200 |             11 |               2 |      0.1056 |              0.008799 |
| PRJNA288562  | SRR6748097  | saliva_GD84                    |             533573 |             13 |               3 |      0.0928 |              0.010918 |
| PRJNA288562  | SRR6748207  | vagina_GD84                    |              45805 |              0 |               0 |      0      |              0        |

## Pregnancy subject T23: longitudinal/multi-site signal

- Gut: *Bifidobacterium adolescentis* C4 persisted from GD84 to GD273; three *E. coli* clusters and three other *Bifidobacterium* clusters were GD84-only.
- Saliva: four *Neisseria* and four *Schaalia* clusters were shared across GD84/GD273; oral profiles were body-site specific.
- Vagina: one low-depth *S. epidermidis* call at GD273; no call at GD84.

## Preterm infant P08: weekly turnover

- W1: one *S. aureus* cluster.
- W2: two low-abundance *E. coli* clusters (`sample_fraction` 5.7e-5 and 6.1e-5).
- W3: six *E. faecalis* clusters.

## PRJNA1517970 low-biomass specificity

- Six vaginal/meconium samples plus one extraction blank produced no strain calls under the 20-species panel.
- Distinct marker counts ranged from 464 (blank) to 116,410 (meconium), so this reflects panel/threshold limits, not runtime failure.
- The blank remained negative, supporting specificity under these gates.

## Isolate pilot / database mismatch evidence

- For *E. coli* isolate SRR32076407, the MSA panel detected **23,999/65,891 E. coli markers at 2.95x depth** with ultra-sensitive gates, but made **no strain call**.
- Runtime was **1.12 s** and peak RSS **205 MiB**.
- Interpretation: species/depth detection works, but study-specific strain resolution requires a cohort-specific isolate-derived panel.

## Caveats

- The MSA panel is not a full oral/vaginal/infant-gut panel; absent calls are not absence of microbes.
- Legacy DBs lack CST; results are exploratory and should not replace the full simulation/mock benchmark.
- Isolate validation requires assembling PRJNA1191225 isolates and rebuilding a cohort-specific panel.
