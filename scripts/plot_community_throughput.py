#!/usr/bin/env python3
"""Community-scale throughput: strain-profiling an S-species community across N samples.
Strain2bScan side is MEASURED end-to-end (results/multispecies_sample_gradient55.tsv: real
multi-profile wall-clock at N=10..500 samples against all 55 species databases; no formula,
no extrapolation). StrainScan has no multi-species mode, so resolving S species means running
it once per species database per sample; that side is PROJECTED from its measured per-run
cost on C. acnes (results/headtohead_performance.tsv, 6.8 s/run) -- the StrainScan line is
explicitly a projection, not a measurement (same-panel head-to-head: see
scripts/run_headtohead_strainscan.py, Linux-only)."""
import os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

def rows(path):
    return [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]

S = 55  # species in the panel (data/accessions/multispecies_55x4.tsv)
SS_PER_RUN = 6.8  # measured StrainScan per-run cost, C. acnes (results/headtohead_performance.tsv)

sg = [r for r in rows("results/multispecies_sample_gradient55.tsv") if r[0] != "n_samples"]
N = [int(r[0]) for r in sg]
s2b = [float(r[1]) for r in sg]           # measured total wall-clock
ss = [n * S * SS_PER_RUN for n in N]      # projected StrainScan (S runs/sample, no multi-species mode)

fig, ax = plt.subplots(figsize=(6.8, 4.4))
ax.plot(N, ss, "o-", color="#8c8c8c", label=f"StrainScan (projected: {S}× runs/sample)")
ax.plot(N, s2b, "s-", color="#1f77b4", label=f"Strain2bScan (measured, {S} real species)")
ax.set_yscale("log"); ax.set_xlabel("# samples"); ax.set_ylabel("total wall-clock (s, log)")
ax.set_title(f"Community-scale strain profiling ({S} real species)")
ax.grid(alpha=0.3, which="both")
ax.legend(loc="lower right")
for n, a, b in zip(N, s2b, ss):
    ax.annotate(f"{b/a:.0f}×", (n, (a*b)**0.5), fontsize=8, color="#d62728", ha="center")
ax.annotate(f"~{s2b[2]/60:.1f} min", (100, s2b[2]), textcoords="offset points", xytext=(0, -14), fontsize=8, color="#1f77b4", ha="center")
ax.annotate(f"~{ss[2]/3600:.1f} h", (100, ss[2]), textcoords="offset points", xytext=(0, 8), fontsize=8, color="#8c8c8c", ha="center")
fig.tight_layout(); os.makedirs("figures", exist_ok=True)
fig.savefig("figures/community_throughput.png", dpi=150); fig.savefig("figures/community_throughput.pdf")
print("wrote figures/community_throughput.png + .pdf")
for n, a, b in zip(N, s2b, ss):
    print(f"  N={n:>4}: Strain2bScan {a:.0f}s ({a/60:.1f} min) | StrainScan {b:.0f}s ({b/3600:.1f} h) | {b/a:.0f}x")
