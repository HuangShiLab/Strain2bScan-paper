#!/usr/bin/env python3
"""Fig 1 - performance: Strain2bScan vs StrainScan per-sample cost + parallel speedup."""
import os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, numpy as np
def rows(path): return [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]
hh = {}
for f in rows("results/headtohead_performance.tsv"):
    if f[0] == "mean": hh = dict(s2t=float(f[1]), s2m=float(f[2]), sst=float(f[3]), ssm=float(f[4]))
par = {}
for f in rows("results/parallel_and_build_scaling.tsv"):
    if f[0] in ("db_build", "profile_sample4"): par[f[0]] = (float(f[1]), float(f[2]), float(f[3]))

fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.2))
tools = ["Strain2bScan", "StrainScan"]
bars = a.bar([0, 1], [hh["s2t"], hh["sst"]], width=0.55, color=["#1f77b4", "#8c8c8c"])
a.set_xticks([0, 1]); a.set_xticklabels(tools); a.set_ylabel("profile time / sample (s)")
a.set_title(f"Per-sample profiling: ~{hh['sst']/hh['s2t']:.0f}× faster, ~{hh['ssm']/hh['s2m']:.0f}× lighter")
for x, t, m in [(0, hh["s2t"], hh["s2m"]), (1, hh["sst"], hh["ssm"])]:
    a.annotate(f"{t:.2f} s\n{m:.0f} MB", (x, t), ha="center", va="bottom", fontsize=9)
a.set_ylim(0, hh["sst"] * 1.25)

labels = ["DB build", "profile"]; t1 = [par["db_build"][0], par["profile_sample4"][0]]
t16 = [par["db_build"][1], par["profile_sample4"][1]]; sp = [par["db_build"][2], par["profile_sample4"][2]]
x = np.arange(2); w = 0.35
b.bar(x - w/2, t1, w, label="1 thread", color="#aec7e8")
b.bar(x + w/2, t16, w, label="16 threads", color="#1f77b4")
b.set_xticks(x); b.set_xticklabels(labels); b.set_ylabel("time (s)")
b.set_title("Parallel speedup (real C. acnes, 64 genomes)")
for i in x: b.annotate(f"{sp[i]:.1f}×", (i, t1[i]), ha="center", va="bottom", fontsize=10, color="#1f77b4")
b.legend()
fig.tight_layout(); os.makedirs("figures", exist_ok=True)
fig.savefig("figures/performance.png", dpi=150); fig.savefig("figures/performance.pdf")
print("wrote figures/performance.png + .pdf")
