#!/usr/bin/env python3
"""Estimate the effect of --trace-gap filtering on existing mock pred files.

This is a Python post-hoc preview of what the new Rust --trace-gap / --trace-floor
flags would produce. It should be replaced by actual Strain2bScan outputs once the
HPC mock rerun with --trace-gap 10 --trace-floor 1e-4 is available.

Inputs:
  data/fig6_fig12_profiles.json
  data/fig6_fig12_metrics.tsv

Outputs:
  results/trace_gap_estimate.tsv
  figures/supplementary/trace_gap_estimate.png
"""
import csv, json, math
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAPER = Path(__file__).resolve().parent.parent
DATA = PAPER / "data"
RES = PAPER / "results"
FIG = PAPER / "figures"

METRICS = list(csv.DictReader(open(DATA / "fig6_fig12_metrics.tsv"), delimiter="\t"))
PROFILES = json.load(open(DATA / "fig6_fig12_profiles.json"))

def filter_trace_gap(prof, min_ratio, floor):
    """Drop the trace tail below the largest abundance gap."""
    items = sorted(prof.items(), key=lambda x: -x[1])
    if len(items) < 2 or min_ratio <= 0:
        cut = len(items)
    else:
        cut = len(items)
        best = min_ratio
        for i in range(len(items) - 1):
            hi, lo = items[i][1], items[i + 1][1]
            if lo > 0 and hi / lo >= best:
                best = hi / lo
                cut = i + 1
    filtered = {g: ab for g, ab in items[:cut] if ab >= floor}
    total = sum(filtered.values()) or 1.0
    return {g: ab / total for g, ab in filtered.items()}

def prf(pred, truth, genomes, threshold=0.0):
    pos = {g for g in genomes if truth.get(g, 0) > 0}
    det = {g for g, v in pred.items() if v >= threshold}
    tp, fp, fn = len(det & pos), len(det - pos), len(pos - det)
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    return prec, rec, f1, tp, fp, fn

def main():
    rows = []
    # Only WMS Strain2bScan default for the trace-gap demonstration.
    for r in METRICS:
        if r["kind"] != "WMS" or r["tool"] != "Strain2bScan" or r["variant"] != "164":
            continue
        sample = r["sample"]
        mock = r["mock"]
        pkey = f"{sample}|Strain2bScan|164"
        tkey = f"{sample}|truth|164"
        prof = PROFILES.get(pkey, {})
        truth = PROFILES.get(tkey, {})
        if not prof or not truth:
            continue
        genomes = set(truth) | set(prof)

        p0, r0, f0, tp0, fp0, fn0 = prf(prof, truth, genomes, 0.0)
        filtered = filter_trace_gap(prof, 10.0, 1e-4)
        p1, r1, f1, tp1, fp1, fn1 = prf(filtered, truth, genomes, 0.0)

        rows.append({
            "sample": sample, "mock": mock,
            "default_TP": tp0, "default_FP": fp0, "default_FN": fn0,
            "default_precision": p0, "default_recall": r0, "default_f1": f0,
            "tracegap_TP": tp1, "tracegap_FP": fp1, "tracegap_FN": fn1,
            "tracegap_precision": p1, "tracegap_recall": r1, "tracegap_f1": f1,
        })

    RES.mkdir(parents=True, exist_ok=True)
    cols = ["sample", "mock",
            "default_TP", "default_FP", "default_FN", "default_precision", "default_recall", "default_f1",
            "tracegap_TP", "tracegap_FP", "tracegap_FN", "tracegap_precision", "tracegap_recall", "tracegap_f1"]
    with open(RES / "trace_gap_estimate.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(r[k], 4) if isinstance(r[k], float) else r[k]) for k in cols})
    print(f"wrote {RES / 'trace_gap_estimate.tsv'}")

    # Aggregate per mock
    mocks = sorted({r["mock"] for r in rows})
    agg = {}
    for m in mocks:
        subset = [r for r in rows if r["mock"] == m]
        def _agg(prefix):
            tp = sum(r[f"{prefix}_TP"] for r in subset)
            fp = sum(r[f"{prefix}_FP"] for r in subset)
            fn = sum(r[f"{prefix}_FN"] for r in subset)
            prec = tp / (tp + fp) if tp + fp else 0.0
            rec = tp / (tp + fn) if tp + fn else 0.0
            f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
            return prec, rec, f1, tp, fp, fn
        agg[m] = {"default": _agg("default"), "tracegap": _agg("tracegap")}

    # Bar chart: precision/recall before/after per mock
    x = np.arange(len(mocks))
    width = 0.2
    fig, ax = plt.subplots(figsize=(8, 5))
    dp = [agg[m]["default"][0] for m in mocks]
    dr = [agg[m]["default"][1] for m in mocks]
    tp = [agg[m]["tracegap"][0] for m in mocks]
    tr = [agg[m]["tracegap"][1] for m in mocks]
    ax.bar(x - 1.5*width, dp, width, label="Default precision", color="#4c78a8")
    ax.bar(x - 0.5*width, dr, width, label="Default recall", color="#72b7b2")
    ax.bar(x + 0.5*width, tp, width, label="Trace-gap precision", color="#f58518")
    ax.bar(x + 1.5*width, tr, width, label="Trace-gap recall", color="#54a24b")
    ax.set_ylabel("Precision / Recall")
    ax.set_xticks(x)
    ax.set_xticklabels(mocks)
    ax.set_ylim(0, 1.05)
    ax.legend(loc="lower right")
    ax.set_title("Estimated effect of --trace-gap 10 --trace-floor 1e-4\n(Python post-process preview)")
    fig.tight_layout()
    out = FIG / "supplementary"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "trace_gap_estimate.png", dpi=200)
    fig.savefig(out / "trace_gap_estimate.pdf")
    plt.close(fig)
    print(f"wrote {out / 'trace_gap_estimate.png'}")

    # Print aggregate table
    print("\n=== Aggregate (micro-averaged) ===")
    print(f"{'mock':<10} {'default P/R/F1':>25} {'trace-gap P/R/F1':>25} {'FP reduction':>15}")
    for m in mocks:
        d = agg[m]["default"]
        t = agg[m]["tracegap"]
        fp_red = f"{(d[4] - t[4])}/{d[4]}" if d[4] else "0/0"
        print(f"{m:<10} {d[0]:.3f}/{d[1]:.3f}/{d[2]:.3f} ({d[3]}/{d[4]}/{d[5]})   "
              f"{t[0]:.3f}/{t[1]:.3f}/{t[2]:.3f} ({t[3]}/{t[4]}/{t[5]})   {fp_red:>15}")


if __name__ == "__main__":
    main()
