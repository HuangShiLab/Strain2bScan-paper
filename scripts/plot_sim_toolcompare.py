#!/usr/bin/env python3
"""Fig12-style tool comparison for simulated communities: Strain2bScan vs StrainScan vs inStrain.

Inputs:
  results/sim_port_comparison.tsv
  figure_raw_data/sim_headtohead/strainscan_single_persample.tsv
  figure_raw_data/sim_headtohead/strainscan_multi_persample.tsv
  results/instrain_sim.tsv

Outputs:
  figures/sim_toolcompare.png
  figures/sim_toolcompare.pdf
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PAPER = Path(__file__).resolve().parent.parent
RES = PAPER / "results"
FIG = PAPER / "figures"
FIGRAW = PAPER / "figure_raw_data" / "sim_headtohead"

SINGLE_SAMPLES = [
    "Escherichia_coli__diff_k2_rep1_d5",
    "Escherichia_coli__same_k2_rep1_d5",
    "Cutibacterium_acnes__diff_k2_rep1_d5",
    "Akkermansia_muciniphila__same_k2_rep1_d5",
    "Staphylococcus_epidermidis__diff_k2_rep1_d5",
    "Prevotella_copri__diff_k2_rep1_d5",
]
MULTI_SAMPLES = ["depth_low_sample01", "depth_med_sample01", "depth_high_sample01"]

COLORS = {
    "Strain2bScan": "#1f77b4",
    "StrainScan": "#8c8c8c",
    "inStrain": "#d62728",
}


def load(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def get_s2b():
    rows = load(RES / "sim_port_comparison.tsv")
    out = {}
    for r in rows:
        if r["mode"] != "auto-depth":
            continue
        sample = r["sample"]
        if sample in SINGLE_SAMPLES or sample in MULTI_SAMPLES:
            out[sample] = {
                "precision": float(r["precision"]),
                "recall": float(r["recall"]),
                "F1": float(r["F1"]),
            }
    return out


def get_strainscan():
    single = load(FIGRAW / "strainscan_single_persample.tsv")
    multi = load(FIGRAW / "strainscan_multi_persample.tsv")
    out = {}
    for r in single:
        sample = r["sample"]
        if sample in SINGLE_SAMPLES:
            out[sample] = {
                "precision": float(r["precision"]),
                "recall": float(r["recall"]),
                "F1": float(r["f1"]),
            }
    for r in multi:
        sample = f"{r['depth']}_{r['sample']}"
        if sample in MULTI_SAMPLES:
            out[sample] = {
                "precision": float(r["precision"]),
                "recall": float(r["recall"]),
                "F1": float(r["f1"]),
            }
    return out


def get_instrain():
    if not (RES / "instrain_sim.tsv").exists():
        return {}
    rows = load(RES / "instrain_sim.tsv")
    out = {}
    for r in rows:
        sample = r["sample"]
        if sample in SINGLE_SAMPLES or sample in MULTI_SAMPLES:
            out[sample] = {
                "precision": float(r["precision"]),
                "recall": float(r["recall"]),
                "F1": float(r["F1"]),
            }
    return out


def plot_panel(ax, samples, s2b, ss, instrain, title):
    x = np.arange(len(samples))
    width = 0.25
    tools = ["Strain2bScan", "StrainScan", "inStrain"]
    metrics = ["precision", "recall", "F1"]
    metric_labels = ["Precision", "Recall", "F1"]
    offsets = np.linspace(-width, width, len(tools))

    for ti, tool in enumerate(tools):
        data = {"Strain2bScan": s2b, "StrainScan": ss, "inStrain": instrain}[tool]
        vals = {m: [data.get(s, {}).get(m, 0.0) for s in samples] for m in metrics}
        bottoms = np.zeros(len(samples))
        for mi, m in enumerate(metrics):
            ax.bar(x + offsets[ti], vals[m], width, bottom=bottoms,
                   color=COLORS[tool], alpha=0.6 + 0.2*mi,
                   edgecolor="white", linewidth=0.3)
            bottoms += np.array(vals[m])

    ax.set_xticks(x)
    short_labels = [s.replace("__", "\n").replace("_", " ") for s in samples]
    ax.set_xticklabels(short_labels, fontsize=6, rotation=30, ha="right")
    ax.set_ylim(0, 3.3)
    ax.set_ylabel("P + R + F1")
    ax.set_title(title)
    # custom legend
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=COLORS[t], edgecolor="white", label=t) for t in tools]
    ax.legend(handles=handles, loc="upper right", fontsize=7)


def main():
    s2b = get_s2b()
    ss = get_strainscan()
    instrain = get_instrain()

    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    plot_panel(axes[0], SINGLE_SAMPLES, s2b, ss, instrain,
               "Single-species simulated communities")
    plot_panel(axes[1], MULTI_SAMPLES, s2b, ss, instrain,
               "Multi-species simulated communities")

    fig.tight_layout()
    fig.savefig(FIG / "sim_toolcompare.png", dpi=200)
    fig.savefig(FIG / "sim_toolcompare.pdf")
    print(f"wrote {FIG / 'sim_toolcompare.png'} and .pdf")


if __name__ == "__main__":
    main()
