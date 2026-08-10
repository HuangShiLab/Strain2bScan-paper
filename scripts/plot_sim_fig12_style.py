#!/usr/bin/env python3
"""Fig12-style tool comparison for simulated communities.

Layout mirrors Fig12_wms_toolcompare:
  - Left: stacked relative-abundance bars for Ground Truth and Strain2bScan
    prediction, one row per sample.
  - Right: horizontal stacked P/R/F1 bars for Strain2bScan, StrainScan and
    inStrain.

Note: Only Strain2bScan produces per-sample abundance profiles in this
benchmark. StrainScan and inStrain outputs are summarised by precision/recall/F1
on the right-hand side.

Inputs:
  figure_raw_data/sim_single_species/<species>/truth/<sample>.truth.tsv
  figure_raw_data/sim_multi_species/<depth>/truth/sample01.truth.tsv
  work/mock_retest/Strain2bScan-raw-data/sim_benchmark/results/{single,multi}/default/<sample>.pred
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
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PAPER = Path(__file__).resolve().parent.parent
RES = PAPER / "results"
FIG = PAPER / "figures"
FIGRAW = PAPER / "figure_raw_data" / "sim_headtohead"
TRUTH_SINGLE = PAPER / "figure_raw_data" / "sim_single_species"
TRUTH_MULTI = PAPER / "figure_raw_data" / "sim_multi_species"
PRED_ROOT = PAPER / "work" / "mock_retest" / "Strain2bScan-raw-data" / "sim_benchmark" / "results"

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
TOOL_COLORS = {
    "Strain2bScan": "#1f77b4",
    "StrainScan": "#8c8c8c",
    "inStrain": "#d62728",
}
METRIC_ORDER = ["precision", "recall", "F1"]
METRIC_ALPHAS = {"precision": 0.55, "recall": 0.75, "F1": 0.95}
METRIC_NAMES = {"precision": "P", "recall": "R", "F1": "F1"}

# tab20-based palette for truth/pred elements
PALETTE = plt.cm.tab20(np.linspace(0, 1, 20))[:, :3]
FP_COLOR = "#bdbdbd"  # grey for false positives / unmatched prediction entries


def load(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def load_truth_single(sample_id):
    """Return dict cluster -> relative_abundance for a single-species sample."""
    species = sample_id.split("__")[0]
    path = TRUTH_SINGLE / species / "truth" / f"{sample_id}.truth.tsv"
    out = {}
    with open(path) as fh:
        next(fh)  # header
        for ln in fh:
            p = ln.rstrip().split("\t")
            if len(p) >= 3:
                out[p[1]] = float(p[2])  # cluster -> relative_abundance
    return out


def load_truth_multi(sample_id):
    """Return dict species -> summed relative_abundance for a multi-species sample."""
    depth = sample_id.replace("_sample01", "")
    path = TRUTH_MULTI / depth / "truth" / "sample01.truth.tsv"
    out = defaultdict(float)
    with open(path) as fh:
        next(fh)
        for ln in fh:
            p = ln.rstrip().split("\t")
            if len(p) >= 4:
                species = p[0].replace("_", " ")
                out[species] += float(p[3])
    return dict(out)


def load_pred_single(sample_id):
    """Return dict cluster -> abundance for a Strain2bScan single-species pred."""
    path = PRED_ROOT / "single" / "default" / f"{sample_id}.pred"
    out = {}
    with open(path) as fh:
        next(fh)
        for ln in fh:
            p = ln.rstrip().split("\t")
            if len(p) >= 2:
                out[p[0]] = float(p[1])
    return out


def load_pred_multi(sample_id):
    """Return dict species -> global_abundance for a Strain2bScan multi-species pred."""
    path = PRED_ROOT / "multi" / "default" / f"{sample_id}.pred"
    out = {}
    with open(path) as fh:
        next(fh)
        for ln in fh:
            p = ln.rstrip().split("\t")
            if len(p) >= 7:
                species = p[0].replace("_", " ")
                out[species] = float(p[6])  # global_abundance
    return out


def build_color_map(all_entries):
    """Assign a stable color to every unique entry name."""
    out = {}
    for i, entry in enumerate(sorted(all_entries)):
        out[entry] = PALETTE[i % len(PALETTE)]
    return out


def get_s2b():
    rows = load(RES / "sim_port_comparison.tsv")
    out = {}
    for r in rows:
        if r["mode"] != "auto-depth":
            continue
        sample = r["sample"]
        out[sample] = {m: float(r[m]) for m in METRIC_ORDER}
    return out


def get_strainscan():
    single = load(FIGRAW / "strainscan_single_persample.tsv")
    multi = load(FIGRAW / "strainscan_multi_persample.tsv")
    out = {}
    for r in single:
        out[r["sample"]] = {m: float(r["precision" if m == "precision" else ("recall" if m == "recall" else "f1")]) for m in METRIC_ORDER}
    for r in multi:
        sample = f"{r['depth']}_{r['sample']}"
        out[sample] = {m: float(r["precision" if m == "precision" else ("recall" if m == "recall" else "f1")]) for m in METRIC_ORDER}
    return out


def get_instrain():
    rows = load(RES / "instrain_sim.tsv")
    out = {}
    for r in rows:
        out[r["sample"]] = {m: float(r[m]) for m in METRIC_ORDER}
    return out


def stacked_bar(ax, y, bar_height, entries, values, color_map, fp_color=FP_COLOR):
    """Draw a horizontal stacked bar. entries without a color use fp_color."""
    x0 = 0.0
    for entry, val in zip(entries, values):
        color = color_map.get(entry, fp_color)
        ax.barh(y, val, height=bar_height, left=x0, color=color, edgecolor="white", linewidth=0.3)
        x0 += val


def main():
    s2b_metrics = get_s2b()
    ss_metrics = get_strainscan()
    ins_metrics = get_instrain()
    all_metrics = {
        "Strain2bScan": s2b_metrics,
        "StrainScan": ss_metrics,
        "inStrain": ins_metrics,
    }

    all_groups = [("Single-species", SINGLE_SAMPLES), ("Multi-species", MULTI_SAMPLES)]

    # Load truth and pred profiles, and build a global color map
    truth_profiles = {}
    pred_profiles = {}
    all_entries = set()
    for group_name, samples in all_groups:
        for sample_id, _ in samples:
            if group_name == "Single-species":
                truth = load_truth_single(sample_id)
                pred = load_pred_single(sample_id)
            else:
                truth = load_truth_multi(sample_id)
                pred = load_pred_multi(sample_id)
            truth_profiles[sample_id] = truth
            pred_profiles[sample_id] = pred
            all_entries.update(truth.keys())
            all_entries.update(pred.keys())

    color_map = build_color_map(all_entries)

    n_rows = sum(len(samples) for _, samples in all_groups)
    fig = plt.figure(figsize=(16, 0.85 * n_rows + 2.5))

    # Left axes: truth + S2bS pred abundance bars
    ax_left = fig.add_axes([0.13, 0.08, 0.37, 0.84])
    # Right axes: P/R/F1 performance bars
    ax_right = fig.add_axes([0.54, 0.08, 0.20, 0.84])
    # Far-right axes: metric table
    ax_table = fig.add_axes([0.78, 0.08, 0.20, 0.84])
    ax_table.set_axis_off()

    y_positions = []
    y_labels = []
    group_lines = []
    current_y = 0

    for group_name, samples in all_groups:
        group_lines.append((current_y - 0.5, group_name))
        for sample_id, sample_label in samples:
            y_positions.append(current_y)
            y_labels.append(sample_label)
            current_y += 1
        current_y += 1.0

    sample_ids_flat = []
    for _, samples in all_groups:
        for sid, _ in samples:
            sample_ids_flat.append(sid)

    # --- left abundance bars ---
    bar_height = 0.28
    for yi, sample_id in enumerate(sample_ids_flat):
        y0 = y_positions[yi]
        truth = truth_profiles[sample_id]
        pred = pred_profiles[sample_id]

        # Truth bar (upper)
        t_entries = sorted(truth.keys(), key=lambda k: truth[k], reverse=True)
        t_vals = [truth[k] for k in t_entries]
        stacked_bar(ax_left, y0 - 0.15, bar_height, t_entries, t_vals, color_map)

        # S2bS pred bar (lower)
        p_entries = sorted(pred.keys(), key=lambda k: pred[k], reverse=True)
        p_vals = [pred[k] for k in p_entries]
        stacked_bar(ax_left, y0 + 0.15, bar_height, p_entries, p_vals, color_map)

    ax_left.set_yticks(y_positions)
    ax_left.set_yticklabels(y_labels, fontsize=8)
    ax_left.set_xlim(0, 1.05)
    ax_left.set_xlabel("Relative abundance", fontsize=9)
    ax_left.set_title("Ground truth vs Strain2bScan prediction", fontsize=11, pad=10)
    ax_left.invert_yaxis()

    # annotation for abundance bars
    ax_left.text(
        0.98, 0.98,
        "Upper bar = ground truth\nLower bar = Strain2bScan\nGrey = prediction-only (FP)",
        ha="right", va="top", fontsize=7,
        transform=ax_left.transAxes,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.8", alpha=0.9),
    )

    # group labels
    y_min, y_max = ax_left.get_ylim()
    for center_y, name in group_lines:
        ax_y = (center_y - y_min) / (y_max - y_min)
        ax_left.text(
            -0.10, ax_y, name,
            ha="right", va="center", fontsize=9, fontweight="bold",
            transform=ax_left.transAxes,
        )

    # --- right performance bars (P + R + F1) ---
    for yi, sample_id in enumerate(sample_ids_flat):
        y0 = y_positions[yi]
        for ti, tool in enumerate(TOOLS):
            vals = all_metrics[tool].get(sample_id, {})
            x0 = 0.0
            for mi, metric in enumerate(METRIC_ORDER):
                v = vals.get(metric, 0.0)
                ax_right.barh(
                    y0 + (ti - 1) * 0.22,
                    v,
                    height=0.18,
                    left=x0,
                    color=TOOL_COLORS[tool],
                    alpha=METRIC_ALPHAS[metric],
                    edgecolor="white",
                    linewidth=0.3,
                )
                x0 += v

    ax_right.set_yticks(y_positions)
    ax_right.set_yticklabels([])
    ax_right.set_xlim(0, 3.15)
    ax_right.set_xlabel("P + R + F1", fontsize=9)
    ax_right.set_title("Performance", fontsize=11, pad=10)
    ax_right.invert_yaxis()

    from matplotlib.patches import Patch
    # performance legend
    perf_handles = []
    for tool in TOOLS:
        for metric in METRIC_ORDER:
            perf_handles.append(
                Patch(facecolor=TOOL_COLORS[tool], alpha=METRIC_ALPHAS[metric],
                      edgecolor="white", label=f"{tool} {METRIC_NAMES[metric]}")
            )
    ax_right.legend(handles=perf_handles, loc="lower right", fontsize=6,
                    ncol=3, frameon=True)

    # --- far-right metric table ---
    sample_col_w = 0.34
    tool_block_w = (1.0 - sample_col_w) / len(TOOLS)
    metric_w = tool_block_w / 3.0
    header_x = [sample_col_w + tool_block_w * i for i in range(len(TOOLS))]
    sub_x_offsets = [0, metric_w, 2 * metric_w]

    ax_table.text(0.02, 1.02, "Sample", ha="left", va="bottom",
                  fontsize=7, fontweight="bold", transform=ax_table.transAxes)
    for ti, tool in enumerate(TOOLS):
        x = header_x[ti]
        ax_table.text(x + tool_block_w / 2, 1.02, tool, ha="center", va="bottom",
                      fontsize=7, fontweight="bold", transform=ax_table.transAxes)
        for mi, metric in enumerate(METRIC_ORDER):
            ax_table.text(x + sub_x_offsets[mi] + metric_w / 2, 0.98, METRIC_NAMES[metric],
                          ha="center", va="top", fontsize=6, transform=ax_table.transAxes)

    row_y_start = 0.94
    row_dy = 1.0 / (n_rows + 2)
    label_idx = 0
    y_idx = 0
    for group_name, samples in all_groups:
        for sid, _ in samples:
            y = row_y_start - y_idx * row_dy
            if label_idx % 2 == 1:
                ax_table.axhspan(y - row_dy * 0.45, y + row_dy * 0.45,
                                 xmin=0, xmax=1, color="0.95", transform=ax_table.transAxes)
            ax_table.text(0.02, y, y_labels[label_idx], ha="left", va="center",
                          fontsize=6, transform=ax_table.transAxes)
            for ti, tool in enumerate(TOOLS):
                vals = all_metrics[tool].get(sid, {})
                x_base = header_x[ti]
                for mi, metric in enumerate(METRIC_ORDER):
                    v = vals.get(metric)
                    text = f"{v:.2f}" if v is not None else "—"
                    ax_table.text(x_base + sub_x_offsets[mi] + metric_w / 2, y, text,
                                  ha="center", va="center", fontsize=6,
                                  transform=ax_table.transAxes)
            label_idx += 1
            y_idx += 1
        y_idx += 0.6

    fig.savefig(FIG / "sim_fig12_style.png", dpi=200)
    fig.savefig(FIG / "sim_fig12_style.pdf")
    print(f"wrote {FIG / 'sim_fig12_style.png'} and .pdf")


if __name__ == "__main__":
    main()
