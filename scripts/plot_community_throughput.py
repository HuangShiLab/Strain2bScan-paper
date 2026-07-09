#!/usr/bin/env python3
"""Community-scale throughput (55 real species, strand-fixed): strain-profiling across N samples.
Strain2bScan side is MEASURED end-to-end at N=10..200 (results/community_throughput.tsv, from
the real multi-profile sample-gradient); N=500 is a linear extrapolation of the steady-state
per-sample cost. StrainScan has no multi-species mode, so resolving S species means running it
once per species database per sample; that side is PROJECTED from its measured per-run cost
(6.8 s/run) -- explicitly a projection, not a measurement."""
import os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

def rows(path):
    return [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]

r = [x for x in rows("results/community_throughput.tsv") if x[0] != "n_samples"]
N = [int(x[0]) for x in r]
s2b = [float(x[1]) for x in r]
ss = [float(x[3]) for x in r]

fig, ax = plt.subplots(figsize=(6.8, 4.4))
ax.plot(N, ss, "o-", color="#8c8c8c", label="StrainScan (projected: 55× runs/sample)")
ax.plot(N, s2b, "s-", color="#1f77b4", label="Strain2bScan (measured, 55 real species)")
ax.set_yscale("log"); ax.set_xlabel("# samples"); ax.set_ylabel("total wall-clock (s, log)")
ax.set_title("Community-scale strain profiling (55 real species)")
ax.grid(alpha=0.3, which="both")
ax.legend(loc="lower right")
for n, a, b in zip(N, s2b, ss):
    ax.annotate(f"{b/a:.0f}×", (n, (a*b)**0.5), fontsize=8, color="#d62728", ha="center")
i100 = N.index(100)
ax.annotate(f"~{s2b[i100]/60:.1f} min", (100, s2b[i100]), textcoords="offset points", xytext=(0, -14), fontsize=8, color="#1f77b4", ha="center")
ax.annotate(f"~{ss[i100]/3600:.1f} h", (100, ss[i100]), textcoords="offset points", xytext=(0, 8), fontsize=8, color="#8c8c8c", ha="center")
fig.tight_layout(); os.makedirs("figures", exist_ok=True)
fig.savefig("figures/community_throughput.png", dpi=150); fig.savefig("figures/community_throughput.pdf")
print("wrote figures/community_throughput.png + .pdf")
for n, a, b in zip(N, s2b, ss):
    print(f"  N={n:>4}: Strain2bScan {a:.0f}s ({a/60:.1f} min) | StrainScan {b:.0f}s ({b/3600:.1f} h) | {b/a:.0f}x")
