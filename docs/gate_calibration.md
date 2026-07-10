# Layer-1 species-gate calibration (breadth term)

**Question.** Which species should Strain2bScan strain-profile in each sample, and how to set the
gate — in particular the new breadth term `--min-species-marker-frac` (see
`src/main.rs::species_tier`)? The gate decides from **absolute species-specific marker evidence**
(never relative abundance):

```
resolve_gate = max(--min-species-markers, ceil(--min-species-marker-frac × total_specific))
detect_gate  = min(--min-species-detect, resolve_gate)
```

**Setup.** Real 55-species panel (218 genomes), 14-enzyme DBs. 20 simulated communities of 12
species each, at two depth regimes: **normal** (paper distribution, median ~2×, min 1×) and
**low** (median 0.62×, 70 % of species sub-1×). A species counts as a *positive call* only if it
is strain-**resolved** (tier `Resolved`) — exactly what the gate controls. Metrics aggregated over
all samples; `detect=10` throughout. Data: `results/gate_calibration_*.tsv`; figure
`figures/gate_calibration.{png,pdf}`; runners `scripts/calibrate_gate.py` (frac sweep, env:
`FLOOR`/`DETECT`/`FRACS`/`SAMPLES_DIR`/`TRUTH_DIR`) and `scripts/sim_lowdepth.py` (low-depth
samples); plot `scripts/plot_gate_calibration.py`.

## Results

| Pass | floor | depth | frac=0 | key finding |
|---|---|---|---|---|
| 1 | 200 | normal | **P=1.00, R=1.00** | P=1.0 at every frac; R=1.0 through frac≤0.06, then declines |
| 2 | 200 | low | **P=1.00, R=0.97** | P=1.0 at every frac; recall loss is all in sub-1× species (**≥1× recall = 0.98**) |
| 3 | 50 | normal | **P=0.97 (8 FP)**, R=1.00 | floor too low → leakage; **frac=0.02 restores P=1.00 at R=1.00** |
| 4 | 50 | low | **P=0.98 (4 FP)**, R=1.00 | frac=0.01 restores P=1.00 (R=0.95) |

## What the calibration shows

1. **Precision is set by the gate (floor/frac); recall is set by depth.** At the default
   **floor=200 the panel never leaks — species precision is 1.0 at both depths and every frac.**
   The 50–60 leakage species per run are correctly parked in *detected-not-resolvable*, not
   resolved.
2. **At floor=200 the breadth term only trades recall** (Fig, panel A): raising frac just pushes
   faint species from *resolved* into *detected-not-resolvable*. So at the default floor the
   optimal breadth is **frac = 0**.
3. **The breadth term's real job is scale insurance** (Fig, panel B). When the absolute floor is
   relaxed to 50, a fixed floor is outrun by cross-species leakage (precision 0.968, 8 false
   species at normal depth); a small **frac = 0.02 restores precision to 1.0 at no recall cost**,
   because the frac term raises the bar in proportion to each species' panel size — exactly where
   large-panel leakage lives. The same mechanism guards precision as a panel grows toward
   open-world size, where any fixed floor would eventually be outrun.
4. **A lower-floor + frac combo does *not* beat the default at equal precision.** (This corrected
   our initial guess.) At P=1.0, low-depth species recall is **0.967** for floor=200/frac=0 versus
   **0.954** for floor=50/frac=0.01 — the frac needed to rescue precision costs slightly more
   recall than the lower floor gains. floor=200/frac=0 is the better P=1.0 operating point here.

## Recommended defaults (validated)

- **`--min-species-markers 200`** and **`--min-species-marker-frac 0.0`** — the shipped defaults.
  Calibration validates them: species precision 1.0 at both depths, best recall among P=1.0
  configurations, and (from the depth benchmark) strain resolution for ≥0.5× species.
- **`--min-species-marker-frac ≈ 0.02–0.05`** is the recommended setting **only if the absolute
  floor is lowered or the panel is much larger than 55 species / open-world**, where a fixed floor
  can be outrun by leakage — it holds precision at 1.0 with minimal recall cost. Re-calibrate per
  enzyme set and depth with `calibrate_gate.py`.

## Practical guidance for low-abundance species

The honest three-tier behavior answers the low-abundance question directly: a faint species with
too few observed species-specific markers is reported **detected-not-resolvable** (species-level
present, no strain claim) rather than dropped or over-resolved. To resolve fainter species, lower
`--min-species-markers` and accept the precision/recall trade-off (optionally add a small
`--min-species-marker-frac` to hold precision); do **not** gate on relative abundance, which
conflates community composition with sequencing depth.
