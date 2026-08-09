#!/usr/bin/env python3
"""Compare Strain2bScan default output with the real --trace-gap rerun on mocks.

Inputs:
  work/mock_retest/Strain2bScan-raw-data/wms_analysis/{wms164,brad164}/
  work/mock_retest/Strain2bScan-raw-data/wms_analysis_tracegap/{wms164,brad164}/
  work/mock_retest/Strain2bScan-raw-data/MSA_combined164_all_flat.members.tsv
  work/mock_retest/Strain2bScan-raw-data/MSA_combined164_bcgi_cont.members.tsv
  work/mock_retest/Strain2bScan-raw-data/Ground_truth/*_ground_truth.txt

Outputs:
  results/trace_gap_estimate.tsv
  figures/supplementary/trace_gap_estimate.png
  figures/supplementary/trace_gap_estimate.pdf
"""
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Re-use the parsing / metric helpers from update_mock_figs.py
import update_mock_figs as umf

PAPER = Path(__file__).resolve().parent.parent
RAW = PAPER / "work" / "mock_retest" / "Strain2bScan-raw-data"
RES = PAPER / "results"
FIG = PAPER / "figures" / "supplementary"

MANIFEST = RAW / "wms_analysis_manifest.csv"


def read_pred(pred_path, members_key):
    """Wrapper around update_mock_figs.read_s2b with path construction."""
    return umf.read_s2b(str(pred_path), members_key)


def genome_set_for(kind):
    return umf.genome_set("164_bcgi" if kind == "2bRAD" else "164_all")


def members_key_for(kind):
    return "164_bcgi" if kind == "2bRAD" else "164_all"


def main():
    rows = []
    with open(MANIFEST) as fh:
        for r in csv.DictReader(fh):
            sample = r["sample"]
            kind = r["kind"]          # WMS or 2bRAD
            enzyme = r["enzyme"]
            outdir = r["outdir"]
            # sample is like WMS_MSA1002_0_100ng_1 or BcgI_MSA1002_0_0.1ng_1
            mock = sample.split("_")[1]  # MSA1002, MSA1003, ...

            default_dir = RAW / "wms_analysis" / outdir
            trace_dir = RAW / "wms_analysis_tracegap" / outdir
            default_pred = default_dir / f"{sample}.pred"
            trace_pred = trace_dir / f"{sample}.pred"

            if not default_pred.exists() or not trace_pred.exists():
                continue

            key = members_key_for(kind)
            genomes = genome_set_for(kind)
            try:
                default_prof = read_pred(default_pred, key)
                trace_prof = read_pred(trace_pred, key)
            except ValueError as e:
                print(f"[skip {sample}] {e}")
                continue

            truth = umf.atcc_genomes(mock, genomes)

            m0 = umf.metrics(default_prof, truth, genomes)
            m1 = umf.metrics(trace_prof, truth, genomes)
            p0t, r0t, f0t = umf.pr_at_threshold(default_prof, truth, genomes, 1e-4)
            p1t, r1t, f1t = umf.pr_at_threshold(trace_prof, truth, genomes, 1e-4)

            rows.append({
                "kind": kind, "mock": mock, "sample": sample,
                "default_TP": m0["TP"], "default_FP": m0["FP"], "default_FN": m0["FN"],
                "default_precision": m0["precision"], "default_recall": m0["recall"], "default_f1": m0["f1"],
                "default_aupr": m0["aupr"],
                "default_precision_at_1e_4": p0t, "default_recall_at_1e_4": r0t, "default_f1_at_1e_4": f0t,
                "tracegap_TP": m1["TP"], "tracegap_FP": m1["FP"], "tracegap_FN": m1["FN"],
                "tracegap_precision": m1["precision"], "tracegap_recall": m1["recall"], "tracegap_f1": m1["f1"],
                "tracegap_aupr": m1["aupr"],
                "tracegap_precision_at_1e_4": p1t, "tracegap_recall_at_1e_4": r1t, "tracegap_f1_at_1e_4": f1t,
            })

    RES.mkdir(parents=True, exist_ok=True)
    cols = [k for k in rows[0].keys()] if rows else []
    with open(RES / "trace_gap_estimate.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {RES / 'trace_gap_estimate.tsv'} ({len(rows)} rows)")

    # Aggregate per (kind, mock)
    agg_keys = sorted({(r["kind"], r["mock"]) for r in rows})
    agg = {}
    for kind, mock in agg_keys:
        subset = [r for r in rows if r["kind"] == kind and r["mock"] == mock]
        def _agg(prefix):
            tp = sum(r[f"{prefix}_TP"] for r in subset)
            fp = sum(r[f"{prefix}_FP"] for r in subset)
            fn = sum(r[f"{prefix}_FN"] for r in subset)
            prec = tp / (tp + fp) if tp + fp else 0.0
            rec = tp / (tp + fn) if tp + fn else 0.0
            f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
            return prec, rec, f1, tp, fp, fn
        agg[(kind, mock)] = {"default": _agg("default"), "tracegap": _agg("tracegap")}

    # Bar chart
    labels = [f"{k[0]}\n{k[1]}" for k in agg_keys]
    x = np.arange(len(labels))
    width = 0.2
    fig, ax = plt.subplots(figsize=(10, 5))
    dp = [agg[k]["default"][0] for k in agg_keys]
    dr = [agg[k]["default"][1] for k in agg_keys]
    tp = [agg[k]["tracegap"][0] for k in agg_keys]
    tr = [agg[k]["tracegap"][1] for k in agg_keys]
    ax.bar(x - 1.5*width, dp, width, label="Default precision", color="#4c78a8")
    ax.bar(x - 0.5*width, dr, width, label="Default recall", color="#72b7b2")
    ax.bar(x + 0.5*width, tp, width, label="Trace-gap precision", color="#f58518")
    ax.bar(x + 1.5*width, tr, width, label="Trace-gap recall", color="#54a24b")
    ax.set_ylabel("Precision / Recall")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.05)
    ax.legend(loc="lower right")
    ax.set_title("Effect of --trace-gap 10 --trace-floor 1e-4 on mock communities")
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "trace_gap_estimate.png", dpi=200)
    fig.savefig(FIG / "trace_gap_estimate.pdf")
    plt.close(fig)
    print(f"wrote {FIG / 'trace_gap_estimate.png'}")

    print("\n=== Aggregate (micro-averaged) ===")
    print(f"{'kind/mock':<18} {'default P/R/F1':>28} {'trace-gap P/R/F1':>28} {'FP reduction':>15}")
    for k in agg_keys:
        d = agg[k]["default"]
        t = agg[k]["tracegap"]
        fp_red = f"{d[4] - t[4]}/{d[4]}" if d[4] else "0/0"
        print(f"{k[0]}/{k[1]:<10} {d[0]:.3f}/{d[1]:.3f}/{d[2]:.3f} ({d[3]}/{d[4]}/{d[5]})   "
              f"{t[0]:.3f}/{t[1]:.3f}/{t[2]:.3f} ({t[3]}/{t[4]}/{t[5]})   {fp_red:>15}")


if __name__ == "__main__":
    main()
