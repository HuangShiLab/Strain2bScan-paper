# Clinical oral 2bRAD cohort — exploratory profiling

The PRJNA1131785 clinical cohort (`S_########` aliases, 38 native-2bRAD saliva samples, incl. the
oral-cancer/ECC series) has **no case/control or tumour/normal labels in the public SRA/ENA
metadata** (samples are annotated only `isolation_source=saliva` + an opaque `source_material_id`).
A differential (tumour vs normal) analysis therefore needs the paper's supplementary sample table.
This is an **exploratory demonstration** that Strain2bScan resolves strain-level profiles on the
clinical cohort — not a differential test.

## Result — `results/clinical_exploratory.tsv`
4 clinical samples, `multi-profile` vs the 19-species oral panel (BcgI, gate 50/5):

| sample | markers | species resolved | strain calls | runtime | top species |
|---|---|---|---|---|---|
| S_0313206 | 271k | 15 | 115 | 2.3 s | Neisseria subflava, Haemophilus parainfluenzae, Actinomyces odontolyticus |
| S_1413231 | 288k | 17 | 129 | 1.8 s | Neisseria subflava, Streptococcus mitis, Actinomyces odontolyticus |
| S_8123213 | 370k | 16 | 149 | 2.0 s | Actinomyces odontolyticus, Neisseria subflava, Haemophilus parainfluenzae |
| S_8221232 | 488k | 17 | 158 | 1.8 s | Actinomyces odontolyticus, Neisseria subflava, Rothia mucilaginosa |

Strain2bScan resolves 15–17 oral species and 115–158 strain-level calls per clinical sample in ~2 s,
dominated by the expected oral commensals — confirming the tool is directly applicable to the
clinical cohort. **To pursue the oral-cancer direction (tumour vs normal strain differences), the
per-sample case/control labels from the paper's supplementary table are required.**
