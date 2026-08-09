#!/usr/bin/env python3
"""Fig6-style stacked-bar figure for simulated single/multi-species benchmark,
comparing ground truth vs default (auto-depth) vs --trace-gap predictions.

Inputs:
  work/mock_retest/Strain2bScan-port-results/sim/{single,multi}/auto-depth/*.pred
  work/mock_retest/Strain2bScan-raw-data/sim_benchmark/results/{single,multi}/tracegap/*.pred
  figure_raw_data/sim_single_species/<Species>/truth/<sample>.truth.tsv
  figure_raw_data/sim_multi_species/<depth>/truth/sample01.truth.tsv

Outputs:
  figures/sim_tracegap.png
  figures/sim_tracegap.pdf
"""
import csv
import re
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PAPER = Path(__file__).resolve().parent.parent
DEFAULT_ROOT = PAPER / "work" / "mock_retest" / "Strain2bScan-port-results" / "sim"
TRACE_ROOT = PAPER / "work" / "mock_retest" / "Strain2bScan-raw-data" / "sim_benchmark" / "results"
TRUTH_SINGLE = PAPER / "figure_raw_data" / "sim_single_species"
TRUTH_MULTI = PAPER / "figure_raw_data" / "sim_multi_species"
FIGDIR = PAPER / "figures"

COLORS = plt.cm.tab20(np.linspace(0, 1, 20))
FP_COLOR = "#1a1a1a"


def read_truth_single(path):
    rows = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 3:
                rows.append((p[1], float(p[2])))  # cluster, relative_abundance
    return dict(rows)


def read_truth_multi(path):
    rows = defaultdict(float)
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 4:
                rows[p[0]] += float(p[3])  # species, sum relative_abundance
    return dict(rows)


def read_pred_single(path):
    out = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2:
                try:
                    out[p[0]] = float(p[1])
                except ValueError:
                    continue
    return out


def read_pred_multi(path):
    out = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 9:
                try:
                    out[p[0]] = float(p[8])  # global_abundance
                except ValueError:
                    continue
    return out


def prf(pred, truth):
    pos = set(truth)
    det = set(pred)
    tp = len(pos & det)
    fp = len(det - pos)
    fn = len(pos - det)
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return prec, rec, f1


def collect_single():
    default_dir = DEFAULT_ROOT / "single" / "auto-depth"
    trace_dir = TRACE_ROOT / "single" / "tracegap"
    if not trace_dir.exists():
        raise FileNotFoundError(f"Trace-gap predictions not yet synced: {trace_dir}")

    samples = []
    for pred_path in sorted(default_dir.glob("*.pred")):
        sample = pred_path.stem
        species = sample.split("__")[0]
        truth_path = TRUTH_SINGLE / species / "truth" / f"{sample}.truth.tsv"
        trace_path = trace_dir / f"{sample}.pred"
        if not truth_path.exists() or not trace_path.exists():
            continue
        truth = read_truth_single(truth_path)
        default = read_pred_single(pred_path)
        trace = read_pred_single(trace_path)
        samples.append({
            "label": sample.replace("__", "\n").replace("_", " "),
            "short": sample.split("__")[0].replace("_", " ") + "\n" + sample.split("__")[1].replace("_", " ") if "__" in sample else sample,
            "truth": truth,
            "default": default,
            "trace": trace,
        })
    return samples


def collect_multi():
    default_dir = DEFAULT_ROOT / "multi" / "auto-depth"
    trace_dir = TRACE_ROOT / "multi" / "tracegap"
    if not trace_dir.exists():
        raise FileNotFoundError(f"Trace-gap predictions not yet synced: {trace_dir}")

    samples = []
    depth_order = ["depth_low", "depth_med", "depth_high"]
    for depth in depth_order:
        pred_path = default_dir / f"{depth}_sample01.pred"
        trace_path = trace_dir / f"{depth}_sample01.pred"
        truth_path = TRUTH_MULTI / depth / "truth" / "sample01.truth.tsv"
        if not all(p.exists() for p in (pred_path, trace_path, truth_path)):
            continue
        truth = read_truth_multi(truth_path)
        default = read_pred_multi(pred_path)
        trace = read_pred_multi(trace_path)
        samples.append({
            "label": depth.replace("depth_", "") + "\ndepth",
            "short": depth.replace("depth_", ""),
            "truth": truth,
            "default": default,
            "trace": trace,
        })
    return samples


def plot_group(ax, samples, title):
    # collect all true items across samples for consistent colors
    all_items = sorted({it for s in samples for it in s["truth"]})
    color_map = {it: COLORS[i % len(COLORS)] for i, it in enumerate(all_items)}

    n = len(samples)
    modes = ["truth", "default", "trace"]
    x = np.arange(n)
    width = 0.25
    offsets = np.linspace(-width, width, len(modes))

    for mi, mode in enumerate(modes):
        bottoms = np.zeros(n)
        # true items
        for item in all_items:
            vals = [s[mode].get(item, 0.0) for s in samples]
            if any(v > 0 for v in vals):
                ax.bar(x + offsets[mi], vals, width, bottom=bottoms,
                       color=color_map[item], edgecolor="white", linewidth=0.3)
                bottoms += np.array(vals)
        # false positives (anything not in truth)
        fp_vals = []
        for s in samples:
            pred = s[mode]
            true_set = set(s["truth"])
            fp = sum(v for k, v in pred.items() if k not in true_set)
            fp_vals.append(fp)
        if any(v > 0 for v in fp_vals):
            ax.bar(x + offsets[mi], fp_vals, width, bottom=bottoms,
                   color=FP_COLOR, edgecolor="white", linewidth=0.3)

    # precision/recall annotations
    for si, s in enumerate(samples):
        for mi, mode in enumerate(modes):
            p, r, f1 = prf(s[mode], s["truth"])
            ax.text(x[si] + offsets[mi], 1.02, f"P={p:.2f}\nR={r:.2f}",
                    ha="center", va="bottom", fontsize=5.5, rotation=0)

    ax.set_xticks(x)
    ax.set_xticklabels([s["short"] for s in samples], fontsize=7, rotation=30, ha="right")
    ax.set_ylim(0, 1.18)
    ax.set_ylabel("relative abundance")
    ax.set_title(title)
    ax.set_xlim(x[0] - 0.6, x[-1] + 0.6)


def main():
    single = collect_single()
    multi = collect_multi()

    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    plot_group(axes[0], single, "Single-species mock communities (truth / default / --trace-gap)")
    plot_group(axes[1], multi, "Multi-species mock communities (truth / default / --trace-gap)")

    # legend
    all_truth = sorted({it for s in single + multi for it in s["truth"]})
    handles = [plt.Rectangle((0, 0), 1, 1, color=plt.cm.tab20(i % 20)) for i, _ in enumerate(all_truth)]
    handles.append(plt.Rectangle((0, 0), 1, 1, color=FP_COLOR))
    labels = list(all_truth) + ["false positive"]
    fig.legend(handles, labels, loc="lower center", ncol=min(len(labels), 8),
               fontsize=7, title="strains / species", title_fontsize=8)

    fig.tight_layout(rect=[0, 0.08, 1, 1])
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / "sim_tracegap.png", dpi=200)
    fig.savefig(FIGDIR / "sim_tracegap.pdf")
    print(f"wrote {FIGDIR / 'sim_tracegap.png'} and .pdf")


if __name__ == "__main__":
    main()
