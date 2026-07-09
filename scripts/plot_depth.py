#!/usr/bin/env python3
"""Fig 3 - depth sensitivity: detection vs per-strain sequencing depth, Strain2bScan vs
StrainScan (single C. acnes strain present in both DBs; digital shotgun, error-free)."""
import os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def rows(path): return [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]
d, s10, s3, ss = [], [], [], []
for f in rows("results/depth_sensitivity.tsv"):
    if f[0] == "depth_x": continue
    d.append(float(f[0])); s10.append(int(f[1])); s3.append(int(f[2])); ss.append(int(f[3]))

fig, ax = plt.subplots(figsize=(6.6, 4.2))
ax.step(d, ss, where="mid", marker="o", lw=2, color="#8c8c8c", label="StrainScan (full k-mer, -l low-depth mode)")
ax.step(d, s3, where="mid", marker="s", lw=2, color="#1f77b4", label="Strain2bScan (permissive floor)")
ax.step(d, s10, where="mid", marker="^", lw=2, ls="--", color="#2ca02c", label="Strain2bScan (default floor)")
ax.set_xscale("log"); ax.set_xticks(d); ax.set_xticklabels([f"{x}×" for x in d])
ax.set_yticks([0, 1]); ax.set_yticklabels(["missed", "detected"]); ax.set_ylim(-0.15, 1.15)
ax.set_xlabel("per-strain sequencing depth")
ax.set_title("Detection vs depth (C. acnes strain): both tools match at 0.5× (strand-fixed)")
ax.legend(loc="center left", fontsize=8); ax.grid(alpha=0.3, axis="x")
ax.axvspan(0.08, 0.4, color="orange", alpha=0.07)
ax.text(0.16, 0.55, "both miss\n< ~0.5×", fontsize=8, color="#cc7000", ha="center")
ax.text(2.0, 0.45, "all three curves coincide:\nStrain2bScan now = StrainScan", fontsize=7.5, color="#555", ha="center")
fig.tight_layout(); os.makedirs("figures", exist_ok=True)
fig.savefig("figures/depth_sensitivity.png", dpi=150); fig.savefig("figures/depth_sensitivity.pdf")
print("wrote figures/depth_sensitivity.png + .pdf")
