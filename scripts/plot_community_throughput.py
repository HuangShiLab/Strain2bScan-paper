#!/usr/bin/env python3
"""Community-scale throughput: strain-profiling an S-species community across N samples.
Strain2bScan side is MEASURED (multi-species sample gradient: ~1.24 s/sample, flat in species,
+ ~20 s build for 40 species). StrainScan has no multi-species mode, so resolving S species
means running it once per species database per sample; we PROJECT that from its measured
per-run cost on C. acnes (~6.8 s). Clearly a projection for the StrainScan line."""
import os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

S = 40                    # species in the community
S2B_PER_SAMPLE = 1.24     # measured (Fig 2b, 40 species)
S2B_BUILD = 20            # measured, one-off
SS_PER_RUN = 6.8          # measured StrainScan per-run (C. acnes); one run per species per sample
N = [10, 50, 100, 200, 500]
s2b = [S2B_BUILD + n * S2B_PER_SAMPLE for n in N]
ss = [n * S * SS_PER_RUN for n in N]

fig, ax = plt.subplots(figsize=(6.8, 4.4))
ax.plot(N, ss, "o-", color="#8c8c8c", label=f"StrainScan (projected: {S}× runs/sample)")
ax.plot(N, s2b, "s-", color="#1f77b4", label="Strain2bScan (measured)")
ax.set_yscale("log"); ax.set_xlabel("# samples"); ax.set_ylabel("total wall-clock (s, log)")
ax.set_title(f"Community-scale strain profiling ({S} species)")
ax.grid(alpha=0.3, which="both")
ax.legend(loc="lower right")
for n, a, b in zip(N, s2b, ss):
    ax.annotate(f"{b/a:.0f}×", (n, (a*b)**0.5), fontsize=8, color="#d62728", ha="center")
# human-readable anchors
ax.annotate("~2.4 min", (100, s2b[2]), textcoords="offset points", xytext=(0, -14), fontsize=8, color="#1f77b4", ha="center")
ax.annotate("~7.6 h", (100, ss[2]), textcoords="offset points", xytext=(0, 8), fontsize=8, color="#8c8c8c", ha="center")
fig.tight_layout(); os.makedirs("figures", exist_ok=True)
fig.savefig("figures/community_throughput.png", dpi=150); fig.savefig("figures/community_throughput.pdf")
print("wrote figures/community_throughput.png + .pdf")
print(f"ratio ≈ S × (SS_per_run / S2B_per_sample) = {S} × {SS_PER_RUN/S2B_PER_SAMPLE:.1f} = {S*SS_PER_RUN/S2B_PER_SAMPLE:.0f}× (independent of N)")
