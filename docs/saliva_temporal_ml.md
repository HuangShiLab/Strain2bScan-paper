# Saliva temporal stability + ML host-ID (Fig 10, panels D–E)

Complements the individual-discrimination result (`docs/saliva_individual_discrimination.md`) with the
two follow-ups the plan called for: is an individual's strain profile **stable across the day**, and
can we **classify the host** from strain features? Same 8 subjects × 4 timepoints (9AM/11AM/1PM/5PM)
native BcgI 2bRAD, oral panel, strain-level abundances. `scripts/temporal_ml_saliva.py`.

## Temporal stability — `results/saliva_temporal_ml.tsv`, `figures/saliva_temporal_ml.*`
Strain-level Bray–Curtis distance **within subject** (across the 4 times of day) vs **between
subjects**:

- within-subject mean **0.186** vs between-subject **0.444** (Mann–Whitney p = 5 × 10⁻²⁵).

An individual's salivary strain profile is ~2.4× more similar across the day than to anyone else's —
the profile is a **stable, person-specific signature**, not a transient snapshot. (Consistent with
timepoint being non-significant in the whole-community PERMANOVA, R² ≈ 0.05.)

## ML host identification — leave-one-timepoint-out
Train a subject classifier on 3 timepoints, predict the held-out 4th (1-NN on Bray–Curtis), cycling
all 4 folds; 8-way problem (chance 12.5%):

- **strain-level accuracy = 100%** (32/32), species-level 93.8%.

Predicting the held-out timepoint means the classifier never sees the test day for that subject, so
this is genuine generalization, not memorization. **Strain-level features identify the host with
perfect accuracy**, beating species-level — reproducing the Strain2bfunc result ("100% host-ID from
strain-level features") with the fast/light Strain2bScan.

Together with the PERMANOVA (strain R²=0.833>0.822) and the *Rothia mucilaginosa* R²=0.921, this
establishes salivary strain signatures as individual-specific and temporally stable, resolvable in
~1 s/sample from native 2bRAD.
