#!/usr/bin/env python3
"""Gate calibration figure: (A) at the default floor, the breadth term only trades recall while
precision stays 1.0; (B) relaxing the absolute floor is where the breadth term earns its keep."""
import os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt


def load(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]
    h = rows[0]
    return [dict(zip(h, r)) for r in rows[1:]]


def col(rows, k):
    return [float(r[k]) for r in rows]


norm = load("results/gate_calibration_normaldepth.tsv")
low = load("results/gate_calibration_lowdepth.tsv")
f50 = load("results/gate_calibration_floor50.tsv") if os.path.exists("results/gate_calibration_floor50.tsv") else None

fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.4))

# --- Panel A: floor=200, breadth term only trades recall (precision flat at 1.0) ---
fr = col(norm, "frac")
a.plot(fr, col(norm, "sp_recall"), "-o", color="#1f77b4", label="species recall (normal depth)")
a.plot(col(low, "frac"), col(low, "sp_recall"), "-o", color="#ff7f0e", label="species recall (low depth, med 0.6×)")
a.plot(col(low, "frac"), col(low, "strain_recall_ge1x"), "--^", color="#ff7f0e", ms=5,
       label="strain recall ≥1× (low depth)")
a.plot(fr, col(norm, "sp_precision"), ":s", color="#2ca02c", label="species precision (both) = 1.0")
a.axvspan(0, 0.06, color="#1f77b4", alpha=0.07)
a.text(0.03, 0.13, "free\nheadroom", ha="center", fontsize=8, color="#1f77b4")
a.set_xlabel("--min-species-marker-frac"); a.set_ylabel("score")
a.set_title("Floor=200: breadth term only trades recall\n(precision stays 1.0 at both depths)", fontsize=10.5)
a.set_ylim(0, 1.06); a.legend(fontsize=7.5, loc="lower left"); a.grid(alpha=0.25)

# --- Panel B: relax the floor (200 -> 50); breadth restores precision if leakage appears ---
if f50:
    fr5 = col(f50, "frac")
    b.plot(fr5, col(f50, "sp_precision"), "-o", color="#d62728", label="precision (floor 50)")
    b.plot(fr5, col(f50, "sp_recall"), "-o", color="#1f77b4", label="recall (floor 50)")
    b.plot(fr, col(norm, "sp_precision"), ":s", color="#2ca02c", label="precision (floor 200)")
    b.plot(fr, col(norm, "sp_recall"), ":^", color="#8c8c8c", label="recall (floor 200)")
    b.set_xlabel("--min-species-marker-frac"); b.set_ylabel("score")
    b.set_title("Relaxing the floor (200→50):\nwhere the breadth term earns its keep", fontsize=10.5)
    b.set_ylim(0, 1.06); b.legend(fontsize=7.5, loc="lower left"); b.grid(alpha=0.25)
else:
    b.text(0.5, 0.5, "floor=50 pass pending", ha="center", va="center")
    b.axis("off")

fig.tight_layout(); os.makedirs("figures", exist_ok=True)
fig.savefig("figures/gate_calibration.png", dpi=150); fig.savefig("figures/gate_calibration.pdf")
print("wrote figures/gate_calibration.png + .pdf")
