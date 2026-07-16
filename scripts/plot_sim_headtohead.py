#!/usr/bin/env python3
"""Fig 11 - systematic Strain2bScan vs StrainScan head-to-head on the 15-species simulated benchmark.
Six panels: (A-C) single-species precision/recall/F1 vs depth; (D) per-species DB build time;
(E) same-environment profile time vs depth; (F) multi-species profiling time per community sample.
Reads results/sim_headtohead_*.tsv; writes figures/sim_headtohead.{png,pdf}."""
import os, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np

S2B, SS = "#1f77b4", "#8c8c8c"        # house colours (blue = Strain2bScan, grey = StrainScan)
DEPTHS = ["0.5", "1", "3", "5", "10"]

def load(p): return list(csv.DictReader(open(p), delimiter="\t"))

bd   = {r["depth"]: r for r in load("results/sim_headtohead_single_by_depth.tsv")}
bc   = [r for r in load("results/sim_headtohead_build_cost.tsv") if r["ss_build_s"] != "NA"]
pt   = load("results/sim_headtohead_profile_sameenv.tsv")
mt   = load("results/sim_headtohead_multi.tsv")

dd = [float(d) for d in DEPTHS if d in bd]
def col(metric, tool): return [float(bd[d][f"{tool}_{metric}"]) for d in DEPTHS if d in bd]

fig = plt.figure(figsize=(13.5, 8))
gs = fig.add_gridspec(2, 3, hspace=0.40, wspace=0.30)

def panel_tag(ax, tag):
    ax.text(-0.14, 1.06, tag, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")

# ---- A-C: accuracy vs depth ----
for j, (metric, title) in enumerate(zip(["precision", "recall", "f1"], ["Precision", "Recall", "F1"])):
    ax = fig.add_subplot(gs[0, j])
    ax.plot(dd, col(metric, "s2b"), "o-", color=S2B, lw=2, ms=7, label="Strain2bScan")
    ax.plot(dd, col(metric, "ss"), "s--", color=SS, lw=2, ms=7, label="StrainScan")
    ax.set_xscale("log"); ax.set_xticks(dd); ax.set_xticklabels([str(x) for x in dd])
    ax.set_xlabel("per-strain sequencing depth (×)"); ax.set_ylabel(f"median {title.lower()}")
    ax.set_title(f"Single-species {title.lower()}"); ax.set_ylim(-0.03, 1.06); ax.grid(alpha=0.3)
    if j == 0: ax.legend(loc="lower right", frameon=True)
    panel_tag(ax, "ABC"[j])

# ---- D: build time per species (log, sorted) ----
axd = fig.add_subplot(gs[1, 0])
bc.sort(key=lambda r: float(r["ss_build_s"]))
y = np.arange(len(bc))
axd.barh(y - 0.2, [float(r["s2b_build_s"]) for r in bc], 0.4, color=S2B, label="Strain2bScan")
axd.barh(y + 0.2, [float(r["ss_build_s"]) for r in bc], 0.4, color=SS, label="StrainScan")
axd.set_yticks(y)
axd.set_yticklabels([r["species"].split("_")[0][0] + ". " + r["species"].split("_")[1] for r in bc], fontsize=7.5, style="italic")
axd.set_xscale("log"); axd.set_xlabel("DB build time (s, log scale)")
axd.set_title("Database build time / species"); axd.legend(fontsize=8, loc="lower right"); axd.grid(alpha=0.3, axis="x")
axd.margins(y=0.02)
panel_tag(axd, "D")

# ---- E: same-env profile time vs depth (log) ----
axe = fig.add_subplot(gs[1, 1])
pdd = [float(r["depth"]) for r in pt]
axe.plot(pdd, [float(r["s2b_emu_wall_s"]) for r in pt], "o-", color=S2B, lw=2, ms=7, label="Strain2bScan")
axe.plot(pdd, [float(r["ss_wall_s"]) for r in pt], "s--", color=SS, lw=2, ms=7, label="StrainScan")
axe.set_xscale("log"); axe.set_yscale("log"); axe.set_xticks(pdd); axe.set_xticklabels([str(x) for x in pdd])
axe.set_xlabel("per-strain sequencing depth (×)"); axe.set_ylabel("profile time / sample (s, log)")
axe.set_title("Profiling cost (same emulated env)"); axe.legend(fontsize=8, loc="upper left"); axe.grid(alpha=0.3)
panel_tag(axe, "E")

# ---- F: multi-species time / community sample (log bars) ----
axf = fig.add_subplot(gs[1, 2])
x = np.arange(len(mt)); w = 0.38
axf.bar(x - w/2, [float(r["s2b_wall_s"]) for r in mt], w, color=S2B, label="Strain2bScan (one pass)")
axf.bar(x + w/2, [float(r["ss_total_wall_s"]) for r in mt], w, color=SS, label="StrainScan ($\\Sigma$ per-species runs)")
axf.set_yscale("log"); axf.set_xticks(x); axf.set_xticklabels([r["depth"].replace("depth_", "") for r in mt])
axf.set_xlabel("community depth"); axf.set_ylabel("time / community sample (s, log)")
axf.set_title("Multi-species profiling time"); axf.legend(fontsize=7.5, loc="upper left"); axf.grid(alpha=0.3, axis="y")
for i, r in enumerate(mt):
    ratio = float(r["ss_total_wall_s"]) / float(r["s2b_wall_s"])
    axf.annotate(f"{ratio:.0f}×", (i, float(r["ss_total_wall_s"])), ha="center", va="bottom", fontsize=8, color="#333")
panel_tag(axf, "F")

os.makedirs("figures", exist_ok=True)
fig.savefig("figures/sim_headtohead.png", dpi=150, bbox_inches="tight")
fig.savefig("figures/sim_headtohead.pdf", bbox_inches="tight")
print("wrote figures/sim_headtohead.png + .pdf")
