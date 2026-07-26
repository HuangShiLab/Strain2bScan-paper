#!/usr/bin/env python3
"""Fig 6 (2bRAD) and Fig 12 (WMS): one row per sample.
  left  = stacked per-genome abundance (truth + each tool), true genomes coloured, FP/decoy grey
  right = six metric bars: precision, recall, F1, AUPR, Bray-Curtis, L2  (tool vs ground truth)

Tool colours use the Okabe-Ito colourblind-safe set (<=4 series, the panel that must be CVD-safe).
Genome composition (20-28 categories, inherently beyond a CVD-safe palette) uses an ordered
qualitative ramp for the TRUE genomes and a single grey bucket for every false-positive/decoy
genome -- so precision loss reads as 'grey creeping into the bar'.
"""
import csv, json, os, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SP = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
RAW = "/Users/macstudio/Downloads/Strain2bScan-raw-data"
FIGDIR = f"{RAW}/MSA_demo"; os.makedirs(FIGDIR, exist_ok=True)

OKABE = {"Strain2bScan": "#0072B2", "Strain2bScan(120)": "#56B4E9",
         "StrainScan": "#E69F00", "inStrain": "#009E73", "truth": "#555555"}
INK = "#222222"; MUTED = "#777777"; GRID = "#DDDDDD"; FPGREY = "#BBBBBB"
METRICS = [("precision", "Precision"), ("recall", "Recall"), ("f1", "F1"),
           ("aupr", "AUPR"), ("bray_curtis", "Bray-Curtis"), ("l2", "L2")]

profiles = json.load(open(f"{SP}/all_profiles.json"))
rows = list(csv.DictReader(open(f"{SP}/all_metrics.tsv"), delimiter="\t"))
def M(sample, tool, variant):
    for r in rows:
        if r["sample"] == sample and r["tool"] == tool and r["variant"] == variant:
            return {k: float(r[k]) for k in ("precision", "recall", "f1", "aupr", "bray_curtis", "l2")}
    return None

def short(g):
    g = re.sub(r"__ATCC_.*", "", g); g = re.sub(r"__GCF_.*", "", g)
    return g.replace("_", " ")

def true_genomes_of(mock, variant_key):
    return [g for g, v in profiles.get(f"{mock_ref(mock)}|truth|{variant_key}", {}).items() if v > 0]

def mock_ref(mock):  # a sample whose truth we stored (any works; truth depends only on mock)
    for r in rows:
        if r["mock"] == mock:
            return r["sample"]
    return None

def genome_palette(true_gs):
    from matplotlib import colormaps
    cmap = colormaps.get_cmap("tab20")
    return {g: cmap(i % 20) for i, g in enumerate(sorted(true_gs, key=short))}

def draw_row(fig, gs_left, gs_right, sample, series, mock, gpal, show_metric_titles):
    """series: list of (label, tool, variant, profile_key). Truth first."""
    axL = fig.add_subplot(gs_left)
    truth_key = f"{sample}|truth|{'120' if any(s[2]=='120' for s in series) else '164'}"
    # build union of genomes to show: true genomes (coloured) + an FP bucket
    true_gs = set(g for g, v in profiles.get(truth_key, {}).items() if v > 0)
    y = np.arange(len(series))
    for i, (label, tool, variant, pkey) in enumerate(series):
        prof = profiles.get(pkey, {})
        tot = sum(prof.values()) or 1.0
        left = 0.0
        # true genomes first, in stable order, then FP bucket
        for g in sorted(true_gs, key=short):
            v = prof.get(g, 0.0) / tot
            if v > 0:
                axL.barh(i, v, left=left, color=gpal[g], edgecolor="white", linewidth=0.5, height=0.72)
                left += v
        fp = sum(v for g, v in prof.items() if g not in true_gs) / tot
        if fp > 0:
            axL.barh(i, fp, left=left, color=FPGREY, edgecolor="white", linewidth=0.5, height=0.72,
                     hatch="///")
    axL.set_yticks(y); axL.set_yticklabels([s[0] for s in series], fontsize=8, color=INK)
    axL.invert_yaxis(); axL.set_xlim(0, 1); axL.set_xticks([0, 0.5, 1.0])
    axL.tick_params(labelsize=7, color=MUTED); axL.set_axisbelow(True)
    for sp in ("top", "right", "left"): axL.spines[sp].set_visible(False)
    axL.spines["bottom"].set_color(MUTED)
    axL.set_ylabel(sample.replace("WMS_", "").replace("BcgI_", ""), fontsize=8,
                   rotation=0, ha="right", va="center", color=INK)

    # right: 6 metric mini-axes
    toolseries = [s for s in series if s[1] != "truth"]
    for j, (mk, mlabel) in enumerate(METRICS):
        ax = fig.add_subplot(gs_right[j])
        vals, cols = [], []
        for (label, tool, variant, pkey) in toolseries:
            m = M(sample, tool, variant)
            vals.append(m[mk] if m else 0.0)
            cols.append(OKABE.get(label, OKABE.get(tool, "#888")))
        xx = np.arange(len(vals))
        disp = [min(v, 1.0) for v in vals] if mk in ("bray_curtis", "l2") else vals
        ax.bar(xx, disp, color=cols, width=0.7, edgecolor="white", linewidth=0.5)
        for k, v in enumerate(vals):
            ax.text(k, min(v, 1.0) + 0.03, f"{v:.2f}", ha="center", va="bottom",
                    fontsize=5.5, color=INK)
        ax.set_ylim(0, 1.15); ax.set_xticks([])
        ax.set_yticks([0, 1] if j == 0 else [])
        ax.tick_params(labelsize=6, color=MUTED)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        ax.spines["left"].set_color(MUTED); ax.spines["bottom"].set_color(MUTED)
        if show_metric_titles:
            ax.set_title(mlabel, fontsize=7.5, color=INK, pad=3)

def build_figure(fig_samples, series_for, title, outbase):
    """fig_samples: list of (sample, mock). series_for(sample,mock)->series list."""
    n = len(fig_samples)
    fig = plt.figure(figsize=(13, 1.15 * n + 1.2))
    outer = fig.add_gridspec(n, 2, width_ratios=[3.1, 4.0], wspace=0.16,
                             hspace=0.55, left=0.13, right=0.985, top=0.93, bottom=0.06)
    # shared genome palette per mock
    palettes = {}
    all_true = set()
    for s, mock in fig_samples:
        vk = "120" if mock in ("MSA1002", "MSA1003") else "164"  # any; truth same genomes
        tg = set(g for g, v in profiles.get(f"{s}|truth|164", {}).items() if v > 0)
        palettes.setdefault(mock, set()).update(tg); all_true.update(tg)
    gpal = genome_palette(all_true)
    for r, (s, mock) in enumerate(fig_samples):
        gl = outer[r, 0]
        gr = outer[r, 1].subgridspec(1, 6, wspace=0.35)
        draw_row(fig, gl, gr, s, series_for(s, mock), mock, gpal, show_metric_titles=(r == 0))
    fig.suptitle(title, fontsize=13, fontweight="bold", color=INK, x=0.5, y=0.985)
    # tool legend
    handles = [Patch(facecolor=OKABE[t], label=t) for t in OKABE if t != "truth"]
    handles += [Patch(facecolor=FPGREY, hatch="///", label="false-positive / decoy genome")]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, 0.005))
    for ext in ("png", "pdf"):
        fig.savefig(f"{outbase}.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outbase}.png / .pdf  ({n} rows)")

# ---------------- Fig 12: WMS, 3 tools ----------------
def wms_series(sample, mock):
    ser = [("truth", "truth", "164", f"{sample}|truth|164"),
           ("Strain2bScan", "Strain2bScan", "164", f"{sample}|Strain2bScan|164")]
    # MSA1003 is the staggered mock where DB expansion costs precision -> show the 20-species tree too
    if mock == "MSA1003" and f"{sample}|Strain2bScan|120" in profiles:
        ser.append(("Strain2bScan(120)", "Strain2bScan", "120", f"{sample}|Strain2bScan|120"))
    ser.append(("StrainScan", "StrainScan", "-", f"{sample}|StrainScan|-"))
    if f"{sample}|inStrain|native" in profiles:
        ser.append(("inStrain", "inStrain", "native", f"{sample}|inStrain|native"))
    return [s for s in ser if s[3] in profiles]

# only WMS samples that have at least one shotgun competitor (StrainScan) -> drop bare rows
def has_competitor(s):
    return f"{s}|StrainScan|-" in profiles
WMS_FIG = [(r["sample"], r["mock"]) for r in rows
           if r["kind"] == "WMS" and r["tool"] == "Strain2bScan" and r["variant"] == "164"
           and has_competitor(r["sample"])]
build_figure(WMS_FIG, wms_series,
             "Figure 12 — Strain-level profiling on shotgun (WMS) mock communities",
             f"{FIGDIR}/fig12_wms_toolcompare")

# ---------------- Fig 6: 2bRAD, Strain2bScan (164 vs 120) ----------------
def brad_series(sample, mock):
    ser = [("truth", "truth", "164", f"{sample}|truth|164"),
           ("Strain2bScan", "Strain2bScan", "164", f"{sample}|Strain2bScan|164")]
    if f"{sample}|Strain2bScan|120" in profiles:
        ser.append(("Strain2bScan(120)", "Strain2bScan", "120", f"{sample}|Strain2bScan|120"))
    return [s for s in ser if s[3] in profiles]

BRAD_FIG = [(r["sample"], r["mock"]) for r in rows
            if r["kind"] == "2bRAD" and r["tool"] == "Strain2bScan" and r["variant"] == "164"]
build_figure(BRAD_FIG, brad_series,
             "Figure 6 — Strain-level profiling on native 2bRAD mock communities",
             f"{FIGDIR}/fig6_2brad")
