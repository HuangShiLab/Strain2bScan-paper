#!/usr/bin/env python3
"""Fig 6 (2bRAD) & Fig 12 (WMS) in the Strain2bFunc house layout.

One row per (tool, sample). Left: horizontal stacked per-genome abundance aligned to that row.
Right: metric columns as horizontal grey bars on the SAME row (value labelled), all pointing
'longer = better': Precision, Recall, F1, AUPR, BC_similarity (=1-BrayCurtis),
L2_similarity (=1-L2/sqrt2). Rows are grouped mock -> tool/data-type -> condition, with a
Ground-Truth reference row (profile only) at the top of each mock.
"""
import csv, json, os, re, math
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import colormaps

PAPER = Path(__file__).resolve().parent.parent
SP = PAPER / "data"
RAW = PAPER / "work" / "mock_retest" / "Strain2bScan-raw-data"
FIGDIR = PAPER / "figures"
os.makedirs(FIGDIR, exist_ok=True)

INK = "#222222"; MUTED = "#6b6b6b"; BARGREY = "#c9c9c9"; FPBLACK = "#1a1a1a"; GROUPLINE = "#888888"
METRICS = [("aupr", "AUPR"), ("precision_at_1e_4", "P@1e-4"),
           ("recall_at_1e_4", "R@1e-4"), ("f1_at_1e_4", "F1@1e-4"),
           ("bc_sim", "BC sim"), ("l2_sim", "L2 sim")]

profiles = json.load(open(SP / "fig6_fig12_profiles.json"))
rows = list(csv.DictReader(open(SP / "fig6_fig12_metrics.tsv"), delimiter="\t"))

def metric_row(sample, tool, variant):
    for r in rows:
        if r["sample"] == sample and r["tool"] == tool and r["variant"] == variant:
            l2, bc = float(r["l2"]), float(r["bray_curtis"])
            def _f(k):
                try:
                    return float(r[k])
                except (KeyError, ValueError):
                    return float("nan")
            return {"aupr": _f("aupr"),
                    "precision_at_1e_4": _f("precision_at_1e_4"),
                    "recall_at_1e_4": _f("recall_at_1e_4"),
                    "f1_at_1e_4": _f("f1_at_1e_4"),
                    "bc_sim": max(0.0, 1 - bc), "l2_sim": max(0.0, 1 - l2 / math.sqrt(2))}
    return None

def species_of(g):
    # 164-panel style: "Genus_species__ATCC_Genus_species_ATCC_12345"
    s = re.sub(r"__(ATCC_|GCF_).*", "", g)
    if s != g:
        return s
    # 120-panel style: "ATCC_Genus_species_ATCC_12345" (no species prefix)
    m = re.match(r"ATCC_(.+?)_ATCC_", g)
    if m:
        return m.group(1)
    return g

def atcc_token(g):
    """Return the terminal ATCC identifier, e.g. ATCC_17978 or ATCC_BAA_816."""
    if "_ATCC_" in g:
        return "ATCC_" + g.rsplit("_ATCC_", 1)[-1]
    m = re.search(r"ATCC_[A-Za-z0-9_]+$", g)
    return m.group(0) if m else None

def is_true_pred(g, true_set, true_tokens):
    """Match by exact name or by terminal ATCC token (handles 120 vs 164 naming)."""
    if g in true_set:
        return True
    tok = atcc_token(g)
    return tok is not None and tok in true_tokens

def sp_label(g):
    return species_of(g).replace("_", " ")

# fixed colour per ATCC target species across the whole figure
def build_palette(all_true):
    sps = sorted({species_of(g) for g in all_true})
    c20, c20b = colormaps.get_cmap("tab20"), colormaps.get_cmap("tab20b")
    pal = {}
    for i, s in enumerate(sps):
        pal[s] = c20(i) if i < 20 else c20b((i - 20) % 20)
    return pal

# ------- row model -------------------------------------------------------------
def truth_key(sample, variant="164"):  # per-variant truth; 164 is the richest/default
    return f"{sample}|truth|{variant}"

def rowspec(mock, group, cond, pkey, sample=None, tool=None, variant=None):
    m = metric_row(sample, tool, variant) if sample else None
    return dict(mock=mock, group=group, cond=cond, pkey=pkey, metrics=m)

def draw_figure(rowlist, title, outbase, legend_species):
    # y positions: top->down, extra gap between groups and bigger gap between mocks
    ys, y = [], 0.0
    prev_mock = prev_group = None
    for r in rowlist:
        if prev_mock is not None and r["mock"] != prev_mock:
            y += 1.4
        elif prev_group is not None and r["group"] != prev_group:
            y += 0.7
        ys.append(y); y += 1.0
        prev_mock, prev_group = r["mock"], r["group"]
    ys = [max(ys) - v for v in ys]          # flip so first row is on top
    H = max(ys) + 1
    all_true = set(g for r in rowlist for g, v in profiles.get(truth_key_of(r), {}).items() if v > 0)
    pal = build_palette(all_true)

    fig = plt.figure(figsize=(16.0, 0.34 * len(rowlist) + 2.2))
    gs = fig.add_gridspec(1, 1 + len(METRICS), width_ratios=[7.2] + [1.0] * len(METRICS),
                          wspace=0.28, left=0.24, right=0.86, top=0.94, bottom=0.06)
    axP = fig.add_subplot(gs[0, 0])
    axM = [fig.add_subplot(gs[0, 1 + j], sharey=axP) for j in range(len(METRICS))]

    for r, yy in zip(rowlist, ys):
        prof = profiles.get(r["pkey"], {})
        tot = sum(prof.values()) or 1.0
        tg = set(g for g, v in profiles.get(truth_key_of(r), {}).items() if v > 0)
        true_tokens = {atcc_token(g) for g in tg if atcc_token(g)}
        left = 0.0
        # draw true-positive segments (match by exact name or terminal ATCC token)
        for g in sorted(prof, key=sp_label):
            if is_true_pred(g, tg, true_tokens):
                v = prof[g] / tot
                if v > 0:
                    axP.barh(yy, v, left=left, height=0.78, color=pal[species_of(g)],
                             edgecolor="white", linewidth=0.4)
                    left += v
        fp = sum(v for g, v in prof.items() if not is_true_pred(g, tg, true_tokens)) / tot
        if fp > 0:
            axP.barh(yy, fp, left=left, height=0.78, color=FPBLACK, edgecolor="white", linewidth=0.4)
        # metric bars — but if the tool detected nothing (empty profile), P/R are undefined, not 0
        empty = (r["metrics"] is not None) and (sum(prof.values()) == 0)
        if empty:
            axM[0].text(3.5, yy, "below detection limit — no strains resolved",
                        va="center", ha="center", fontsize=6.8, style="italic", color=MUTED,
                        clip_on=False, zorder=20, transform=axM[0].get_yaxis_transform())
        elif r["metrics"]:
            for j, (mk, _) in enumerate(METRICS):
                val = r["metrics"][mk]
                axM[j].barh(yy, val, height=0.78, color=BARGREY, edgecolor="none")
                axM[j].text(min(val, 0.98), yy, f"{val:.2f}", va="center",
                            ha="right" if val > 0.18 else "left",
                            fontsize=5.6, color=INK, clip_on=True)

    axP.set_xlim(0, 1); axP.set_ylim(-0.8, H)
    axP.set_xticks([0, .25, .5, .75, 1.0]); axP.set_xlabel("Relative abundance", fontsize=9, color=INK, labelpad=8)
    axP.xaxis.tick_bottom()
    axP.set_yticks(ys); axP.set_yticklabels([r["cond"] for r in rowlist], fontsize=7.2, color=INK)
    axP.tick_params(axis="x", labelsize=7, colors=MUTED)
    for s in ("right", "left", "bottom"): axP.spines[s].set_visible(False)
    axP.spines["top"].set_color(MUTED)

    for j, (mk, mlabel) in enumerate(METRICS):
        ax = axM[j]; ax.set_xlim(0, 1)
        ax.set_facecolor("none")  # transparent so a spanning "below detection" label shows through
        ax.set_title(mlabel, fontsize=6.6, color=INK, pad=10)
        ax.set_xticks([]); ax.tick_params(axis="y", left=False, labelleft=False)
        for s in ("top", "right", "left", "bottom"): ax.spines[s].set_visible(False)
        ax.axvline(0, color="#eeeeee", lw=0.8)

    # group + mock brackets on the far left (figure fraction via axP transform)
    from collections import OrderedDict
    def blocks(key):
        out, start = [], 0
        for i in range(1, len(rowlist) + 1):
            if i == len(rowlist) or rowlist[i][key] != rowlist[start][key] or \
               (key == "group" and rowlist[i]["mock"] != rowlist[start]["mock"]):
                out.append((start, i - 1)); start = i
        return out
    for a, b in blocks("group"):
        yc = (ys[a] + ys[b]) / 2
        axP.annotate(rowlist[a]["group"], xy=(-0.10, yc), xycoords=("axes fraction", "data"),
                     ha="right", va="center", fontsize=7.0, color=INK)
        axP.plot([-0.085, -0.085], [ys[b] - 0.4, ys[a] + 0.4],
                 transform=axP.get_yaxis_transform(), color=GROUPLINE, lw=0.9, clip_on=False)
    for a, b in blocks("mock"):
        yc = (ys[a] + ys[b]) / 2
        axP.annotate(rowlist[a]["mock"], xy=(-0.38, yc), xycoords=("axes fraction", "data"),
                     ha="center", va="center", fontsize=11, fontweight="bold", color=INK, rotation=90)

    fig.suptitle(title, fontsize=11.5, fontweight="bold", color=INK, x=0.01, y=0.998, ha="left")
    handles = [Patch(facecolor=pal[s], label=s.replace("_", " ")) for s in
               sorted(pal, key=lambda x: x.replace("_", " ")) if s in legend_species]
    handles.append(Patch(facecolor=FPBLACK, label="False positive / decoy"))
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.92, 0.5),
               fontsize=6.6, frameon=False, handlelength=1.0, labelspacing=0.35)
    for ext in ("png", "pdf"):
        fig.savefig(f"{outbase}.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outbase}.png/.pdf ({len(rowlist)} rows)")

def truth_key_of(r):
    return truth_key(r["_sample"]) if "_sample" in r else r.get("_truthkey", r["pkey"])

# ---------------- assemble Fig 12 (WMS) ----------------
def add_truth(rowlist, mock, sample):
    rowlist.append({**rowspec(mock, "Ground Truth", "", truth_key(sample)),
                    "_truthkey": truth_key(sample)})

WMS_STRUCT = [
    ("MSA1002", [("0_100ng_2", "0% host"), ("90_100ng_1", "90% host"),
                 ("95_100ng_1", "95% host"), ("99_100ng_1", "99% host")]),
    ("MSA1003", [("0_100ng_1", "rep1"), ("0_100ng_2", "rep2"), ("0_100ng_3", "rep3")]),
    ("MSA1005", [("0_100ng_1", "rep1"), ("0_100ng_2", "rep2"), ("0_100ng_3", "rep3")]),
    ("MSA1007", [("0_100ng_1", "rep1"), ("0_100ng_2", "rep2"), ("0_100ng_3", "rep3")]),
]
def build_wms():
    rl = []
    for mock, samps in WMS_STRUCT:
        add_truth(rl, mock, f"WMS_{mock}_{samps[0][0]}")
        # inStrain reported on its documented DEREPLICATED reference (fair-play): ref98 for the
        # 20-species mocks, ref164-derep for the 28-species mocks.
        itool_variant = "derep98" if mock in ("MSA1002", "MSA1003") else "derep164"
        for tool, variant, tag in [("Strain2bScan", "164", "Strain2bScan"),
                                    ("StrainScan", "-", "StrainScan"),
                                    ("inStrain", itool_variant, "inStrain")]:
            for scode, clabel in samps:
                s = f"WMS_{mock}_{scode}"
                pkey = f"{s}|{tool}|{variant}"
                if pkey not in profiles:
                    continue
                r = rowspec(mock, f"WMS ({tag})", clabel, pkey, s, tool, variant)
                r["_truthkey"] = truth_key(s, variant)
                rl.append(r)
        # strainscan-port WMS results (layers only; default path is identical to Strain2bScan)
        for tool, variant, tag in [("Strain2bScan-port", "164-layers", "S2bS-port layers")]:
            for scode, clabel in samps:
                s = f"WMS_{mock}_{scode}"
                pkey = f"{s}|{tool}|{variant}"
                if pkey not in profiles:
                    continue
                r = rowspec(mock, f"WMS ({tag})", clabel, pkey, s, tool, variant)
                r["_truthkey"] = truth_key(s, variant)
                rl.append(r)
        if mock in ("MSA1002", "MSA1003"):
            for tool, variant, tag in [("Strain2bScan-port", "120-layers", "S2bS-port layers (120)")]:
                for scode, clabel in samps:
                    s = f"WMS_{mock}_{scode}"
                    pkey = f"{s}|{tool}|{variant}"
                    if pkey not in profiles:
                        continue
                    r = rowspec(mock, f"WMS ({tag})", clabel, pkey, s, tool, variant)
                    r["_truthkey"] = truth_key(s, variant)
                    rl.append(r)
    return rl

# ---------------- assemble Fig 6 (2bRAD) ----------------
BRAD_STRUCT = [
    ("MSA1002", [("0_0.001ng_1", "0.001 ng"), ("0_0.01ng_1", "0.01 ng"), ("0_0.1ng_1", "0.1 ng"),
                 ("0_1ng_1", "1 ng"), ("90_100ng_1", "90% host"), ("95_100ng_1", "95% host"),
                 ("99_100ng_1", "99% host"), ("99.9_100ng_1", "99.9% host")]),
    ("MSA1003", [("0_100ng_1", "rep1"), ("0_100ng_2", "rep2"), ("0_100ng_3", "rep3")]),
    ("MSA1005", [("0_100ng_1", "rep1"), ("0_100ng_2", "rep2"), ("0_100ng_3", "rep3")]),
    ("MSA1007", [("0_100ng_1", "rep1"), ("0_100ng_2", "rep2"), ("0_100ng_3", "rep3")]),
]
def build_brad():
    rl = []
    for mock, samps in BRAD_STRUCT:
        add_truth(rl, mock, f"BcgI_{mock}_{samps[0][0]}")
        for scode, clabel in samps:
            s = f"BcgI_{mock}_{scode}"
            pkey = f"{s}|Strain2bScan|164"
            if pkey not in profiles:
                continue
            r = rowspec(mock, "2bRAD (Strain2bScan)", clabel, pkey, s, "Strain2bScan", "164")
            r["_truthkey"] = truth_key(s)
            rl.append(r)
    return rl

# ---------------- assemble Fig S: DB-expansion cost (20- vs 28-species tree) ----------------
# Only the mocks that exist in BOTH trees (MSA1002/1003). WMS modality, matching Fig 12.
SUPP_STRUCT = [
    ("MSA1002", [("0_100ng_2", "0% host"), ("90_100ng_1", "90% host"),
                 ("95_100ng_1", "95% host"), ("99_100ng_1", "99% host")]),
    ("MSA1003", [("0_100ng_1", "rep1"), ("0_100ng_2", "rep2"), ("0_100ng_3", "rep3")]),
]
def build_supp():
    rl = []
    for mock, samps in SUPP_STRUCT:
        add_truth(rl, mock, f"WMS_{mock}_{samps[0][0]}")
        for variant, tag in [("164", "28-species tree"), ("120", "20-species tree")]:
            for scode, clabel in samps:
                s = f"WMS_{mock}_{scode}"
                pkey = f"{s}|Strain2bScan|{variant}"
                if pkey not in profiles:
                    continue
                r = rowspec(mock, tag, clabel, pkey, s, "Strain2bScan", variant)
                r["_truthkey"] = truth_key(s, variant)
                rl.append(r)
    return rl

if __name__ == "__main__":
    supp = build_supp()
    legend_sp = {species_of(g) for r in supp for g, v in profiles.get(r["_truthkey"], {}).items() if v > 0}
    draw_figure(supp, "Figure S — Database-expansion cost: 20- vs 28-species tree",
                f"{FIGDIR}/figS_tree_expansion", legend_sp)
    wms = build_wms()
    legend_sp = {species_of(g) for r in wms for g, v in profiles.get(r["_truthkey"], {}).items() if v > 0}
    draw_figure(wms, "Figure 12 — Strain-level profiling on shotgun (WMS) mock communities",
                f"{FIGDIR}/fig12_wms_toolcompare", legend_sp)
    brad = build_brad()
    legend_sp = {species_of(g) for r in brad for g, v in profiles.get(r["_truthkey"], {}).items() if v > 0}
    draw_figure(brad, "Figure 6 — Strain-level profiling on native 2bRAD mock communities",
                f"{FIGDIR}/fig6_2brad", legend_sp)
