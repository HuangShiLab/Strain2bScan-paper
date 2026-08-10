#!/usr/bin/env python3
"""Fig12-style tool comparison for simulated communities.

Similar layout to Fig12_wms_toolcompare: one row per sample, grouped by
single-species vs multi-species, with stacked P/R/F1 bars per tool and a
right-hand metric table.

Inputs:
  results/sim_port_comparison.tsv
  figure_raw_data/sim_headtohead/strainscan_single_persample.tsv
  figure_raw_data/sim_headtohead/strainscan_multi_persample.tsv
  results/instrain_sim.tsv

Outputs:
  figures/sim_fig12_style.png
  figures/sim_fig12_style.pdf
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
    ("Escherichia_coli__diff_k2_rep1_d5", "E. coli diff-k2"),
    ("Cutibacterium_acnes__diff_k2_rep1_d5", "C. acnes diff-k2"),
    ("Staphylococcus_epidermidis__diff_k2_rep1_d5", "S. epidermidis diff-k2"),
    ("Prevotella_copri__diff_k2_rep1_d5", "P. copri diff-k2"),
]
MULTI_SAMPLES = [
    ("depth_low_sample01", "depth low"),
    ("depth_med_sample01", "depth med"),
    ("depth_high_sample01", "depth high"),
]

TOOLS = ["Strain2bScan", "StrainScan", "inStrain"]
COLORS = {
    "Strain2bScan": "#1f77b4",
    "StrainScan": "#8c8c8c",
    "inStrain": "#d62728",
}
METRIC_ALPHAS = {"precision": 0.55, "recall": 0.75, "F1": 0.95}
METRIC_NAMES = {"precision": "P", "recall": "R", "F1": "F1"}


def load(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def get_s2b():
    rows = load(RES / "sim_port_comparison.tsv")
    out = {}
    for r in rows:
        if r["mode"] != "auto-depth":
            continue
        sample = r["sample"]
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
        out[sample] = {
            "precision": float(r["precision"]),
            "recall": float(r["recall"]),
            "F1": float(r["f1"]),
        }
    for r in multi:
        sample = f"{r['depth']}_{r['sample']}"
        out[sample] = {
            "precision": float(r["precision"]),
            "recall": float(r["recall"]),
            "F1": float(r["f1"]),
        }
    return out


def get_instrain():
    rows = load(RES / "instrain_sim.tsv")
    out = {}
    for r in rows:
        sample = r["sample"]
        out[sample] = {
            "precision": float(r["precision"]),
            "recall": float(r["recall"]),
            "F1": float(r["F1"]),
        }
    return out


def format_val(v):
    if v is None:
        return "—"
    return f"{v:.2f}"


def main():
    s2b = get_s2b()
    ss = get_strainscan()
    ins = get_instrain()
    data_sources = {"Strain2bScan": s2b, "StrainScan": ss, "inStrain": ins}

    all_groups = [("Single-species", SINGLE_SAMPLES), ("Multi-species", MULTI_SAMPLES)]
    n_rows = sum(len(samples) for _, samples in all_groups)

    fig = plt.figure(figsize=(16, 0.85 * n_rows + 2.5))
    # main bar axes on the left, table axes on the right
    ax_bar = fig.add_axes([0.16, 0.08, 0.46, 0.84])
    ax_table = fig.add_axes([0.66, 0.08, 0.32, 0.84])
    ax_table.set_axis_off()

    y_positions = []
    y_labels = []
    group_lines = []
    current_y = 0

    # Table column positions: sample label + 3 tool blocks
    sample_col_w = 0.22
    tool_block_w = (1.0 - sample_col_w) / len(TOOLS)
    metric_w = tool_block_w / 3.0
    header_x = [sample_col_w + tool_block_w * i for i in range(len(TOOLS))]
    sub_x_offsets = [0, metric_w, 2 * metric_w]

    table_texts = []
    for group_name, samples in all_groups:
        # group label sits in the gap just above the first row of the group
        group_lines.append((current_y - 0.5, group_name))
        for sample_id, sample_label in samples:
            y_positions.append(current_y)
            y_labels.append(sample_label)
            row_text = [sample_label]
            for ti, tool in enumerate(TOOLS):
                vals = data_sources[tool].get(sample_id, {})
                row_text.append(format_val(vals.get("precision")))
                row_text.append(format_val(vals.get("recall")))
                row_text.append(format_val(vals.get("F1")))
            table_texts.append(row_text)
            current_y += 1
        current_y += 1.0  # gap between groups

    # --- bars ---
    bar_height = 0.22
    metric_order = ["precision", "recall", "F1"]
    for yi, sample_id in enumerate(y_positions):
        sample_key = y_labels[yi]  # not used; we need original id
    # rebuild mapping for bars
    sample_ids_flat = []
    for _, samples in all_groups:
        for sid, _ in samples:
            sample_ids_flat.append(sid)

    for yi, sample_id in enumerate(sample_ids_flat):
        y0 = y_positions[yi]
        for ti, tool in enumerate(TOOLS):
            vals = data_sources[tool].get(sample_id, {})
            x0 = 0
            for mi, metric in enumerate(metric_order):
                v = vals.get(metric, 0.0)
                ax_bar.barh(
                    y0 + (ti - 1) * bar_height,
                    v,
                    height=bar_height * 0.85,
                    left=x0,
                    color=COLORS[tool],
                    alpha=METRIC_ALPHAS[metric],
                    edgecolor="white",
                    linewidth=0.4,
                )
                x0 += v

    ax_bar.set_yticks(y_positions)
    ax_bar.set_yticklabels(y_labels, fontsize=8)
    ax_bar.set_xlim(0, 3.15)
    ax_bar.set_xlabel("P + R + F1", fontsize=9)
    ax_bar.set_title("Strain-level profiling on simulated communities", fontsize=12, pad=10)
    ax_bar.invert_yaxis()

    # group labels on the left of the y-axis (axes coords)
    y_min, y_max = ax_bar.get_ylim()
    for center_y, name in group_lines:
        # convert data y to axes y (inverted axis handled by ax.transData + axes y)
        ax_y = (center_y - y_min) / (y_max - y_min)
        ax_bar.text(
            -0.08, ax_y, name,
            ha="right", va="center", fontsize=9, fontweight="bold",
            transform=ax_bar.transAxes,
        )

    # legend
    from matplotlib.patches import Patch
    legend_handles = []
    for tool in TOOLS:
        for metric in metric_order:
            legend_handles.append(
                Patch(facecolor=COLORS[tool], alpha=METRIC_ALPHAS[metric],
                      edgecolor="white", label=f"{tool} {METRIC_NAMES[metric]}")
            )
    ax_bar.legend(handles=legend_handles, loc="lower right", fontsize=6,
                  ncol=3, frameon=True)

    # --- table ---
    # header
    ax_table.text(0.01, 1.02, "Sample", ha="left", va="bottom",
                  fontsize=8, fontweight="bold", transform=ax_table.transAxes)
    for ti, tool in enumerate(TOOLS):
        x = header_x[ti]
        ax_table.text(x + tool_block_w / 2, 1.02, tool, ha="center", va="bottom",
                      fontsize=8, fontweight="bold", transform=ax_table.transAxes)
        for mi, metric in enumerate(metric_order):
            ax_table.text(x + sub_x_offsets[mi] + metric_w / 2, 0.98, METRIC_NAMES[metric],
                          ha="center", va="top", fontsize=7,
                          transform=ax_table.transAxes)

    row_y_start = 0.94
    row_dy = 1.0 / (n_rows + 2)
    label_idx = 0
    y_idx = 0
    for group_name, samples in all_groups:
        for sid, _ in samples:
            y = row_y_start - y_idx * row_dy
            # alternate background within group
            if label_idx % 2 == 1:
                ax_table.axhspan(y - row_dy * 0.45, y + row_dy * 0.45,
                                 xmin=0, xmax=1, color="0.95", transform=ax_table.transAxes)
            # sample label
            ax_table.text(0.01, y, y_labels[label_idx], ha="left", va="center",
                          fontsize=7, transform=ax_table.transAxes)
            # values
            for ti, tool in enumerate(TOOLS):
                vals = data_sources[tool].get(sid, {})
                x_base = header_x[ti]
                for mi, metric in enumerate(metric_order):
                    v = vals.get(metric)
                    text = format_val(v)
                    ax_table.text(x_base + sub_x_offsets[mi] + metric_w / 2, y, text,
                                  ha="center", va="center", fontsize=7,
                                  transform=ax_table.transAxes)
            label_idx += 1
            y_idx += 1
        y_idx += 0.6  # visual gap (non-integer)

    fig.savefig(FIG / "sim_fig12_style.png", dpi=200)
    fig.savefig(FIG / "sim_fig12_style.pdf")
    print(f"wrote {FIG / 'sim_fig12_style.png'} and .pdf")


if __name__ == "__main__":
    main()
