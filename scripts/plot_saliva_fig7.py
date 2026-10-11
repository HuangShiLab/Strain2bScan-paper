#!/usr/bin/env python3
"""Build the five-panel saliva Figure 7 from committed source tables.

Panels:
  A species-level PCoA; B strain-level PCoA; C per-species strain subject R2;
  D within- versus between-subject Bray-Curtis; E leave-one-timepoint-out host ID.

The script reads PERMANOVA values rather than recomputing them, but independently recomputes
Bray-Curtis distance and the leave-one-timepoint-out classifier and checks both against
`results/saliva_temporal_ml.tsv`.
"""
from pathlib import Path
import csv
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PAPER = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PAPER / "scripts"))
from saliva_table import read_long, COMMUNITY_ABUNDANCE  # noqa: E402

RNG = np.random.default_rng(20261011)


def read_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def matrix(rows, level):
    if level == "strain":
        key = lambda r: f"{r['species']}|{r['cluster']}"
    else:
        key = lambda r: r["species"]
    samples = sorted({r["sample"] for r in rows})
    features = sorted({key(r) for r in rows})
    si = {s: i for i, s in enumerate(samples)}
    fi = {f: j for j, f in enumerate(features)}
    M = np.zeros((len(samples), len(features)))
    for r in rows:
        M[si[r["sample"]], fi[key(r)]] += float(r[COMMUNITY_ABUNDANCE])
    total = M.sum(axis=1, keepdims=True)
    total[total == 0] = 1
    return samples, M / total


def bray_curtis(M):
    n = len(M)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            denominator = (M[i] + M[j]).sum()
            value = abs(M[i] - M[j]).sum() / denominator if denominator else 0.0
            D[i, j] = D[j, i] = value
    return D


def pcoa(D):
    n = len(D)
    A = -0.5 * D * D
    J = np.eye(n) - np.ones((n, n)) / n
    B = J @ A @ J
    values, vectors = np.linalg.eigh(B)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[order]
    coords = vectors[:, :2] * np.sqrt(np.maximum(values[:2], 0))
    positive = values[values > 0]
    explained = values[:2] / positive.sum() if len(positive) else np.array([0.0, 0.0])
    return coords, explained


def loto_accuracy(M, subject, timepoint):
    D = bray_curtis(M)
    correct = 0
    for held in sorted(set(timepoint)):
        test = np.where(timepoint == held)[0]
        train = np.where(timepoint != held)[0]
        for i in test:
            nearest = train[np.argmin(D[i, train])]
            correct += subject[nearest] == subject[i]
    return correct / len(subject)


def main():
    rows = read_long(PAPER / "results/saliva_strain_long.tsv",
                     require=(COMMUNITY_ABUNDANCE, "within_abund"))
    samples, Ms = matrix(rows, "strain")
    _, Mp = matrix(rows, "species")
    subject = np.array([s.split("-")[1] for s in samples])
    timepoint = np.array([s.split("-")[0] for s in samples])

    permanova = read_tsv(PAPER / "results/saliva_permanova.tsv")
    perm = {(r["factor"], r["level"]): r for r in permanova if r["factor"] == "subject"}
    temporal = {r["metric"]: r for r in read_tsv(PAPER / "results/saliva_temporal_ml.tsv")}
    acc_s2b = loto_accuracy(Ms, subject, timepoint)
    acc_species = loto_accuracy(Mp, subject, timepoint)
    expected_s2b = float(temporal["leave_one_timepoint_out_subject_accuracy"]["strain"])
    expected_species = float(temporal["leave_one_timepoint_out_subject_accuracy"]["species"])
    if not math.isclose(acc_s2b, expected_s2b, abs_tol=5e-5) or not math.isclose(acc_species, expected_species, abs_tol=5e-5):
        raise SystemExit(f"LOTO mismatch: computed {acc_species}/{acc_s2b}, archived {expected_species}/{expected_s2b}")

    palette = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b",
               "#17becf", "#e377c2"]
    subjects = sorted(set(subject), key=int)
    colors = {s: palette[i % len(palette)] for i, s in enumerate(subjects)}
    coords_s, ev_s = pcoa(bray_curtis(Mp))
    coords_b, ev_b = pcoa(bray_curtis(Ms))

    per_species = sorted(read_tsv(PAPER / "results/saliva_perspecies_subject.tsv"),
                         key=lambda r: float(r["strain_R2_subject"]))
    D = bray_curtis(Ms)
    within = np.array([D[i, j] for i in range(len(samples)) for j in range(i + 1, len(samples))
                       if subject[i] == subject[j]])
    between = np.array([D[i, j] for i in range(len(samples)) for j in range(i + 1, len(samples))
                        if subject[i] != subject[j]])
    expected_within = float(temporal["within_subject_BC_mean"]["strain"])
    expected_between = float(temporal["between_subject_BC_mean"]["strain"])
    if not math.isclose(within.mean(), expected_within, abs_tol=5e-5) or not math.isclose(between.mean(), expected_between, abs_tol=5e-5):
        raise SystemExit("Bray-Curtis group means do not match archived saliva table")
    p_within = float(temporal["within_vs_between_mannwhitney_p"]["strain"])

    fig = plt.figure(figsize=(15.5, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 6)
    axes = {
        "A": fig.add_subplot(gs[0, 0:2]),
        "B": fig.add_subplot(gs[0, 2:4]),
        "C": fig.add_subplot(gs[0, 4:6]),
        "D": fig.add_subplot(gs[1, 0:3]),
        "E": fig.add_subplot(gs[1, 3:6]),
    }

    for label, coords, explained, level in (("A", coords_s, ev_s, "species"), ("B", coords_b, ev_b, "strain")):
        ax = axes[label]
        row = perm[("subject", level)]
        for s in subjects:
            idx = np.where(subject == s)[0]
            ax.scatter(coords[idx, 0], coords[idx, 1], s=48, color=colors[s],
                       label=f"{s}" if label == "A" else None, edgecolor="k", linewidth=0.35, alpha=0.9)
        ax.set_title(f"{label}  {level.capitalize()}-level PCoA\nsubject R²={float(row['R2']):.3f}, p={float(row['p_value']):.1e}")
        ax.set_xlabel(f"PCo1 ({explained[0]:.0%})")
        ax.set_ylabel(f"PCo2 ({explained[1]:.0%})")
        ax.grid(alpha=0.22)
    axes["A"].legend(title="Subject", fontsize=8, loc="best", ncol=2)

    ax = axes["C"]
    ypos = np.arange(len(per_species))
    r2 = [float(r["strain_R2_subject"]) for r in per_species]
    sig = [float(r["p_value"]) < 0.05 for r in per_species]
    ax.barh(ypos, r2, color=["#1f77b4" if x else "#bbbbbb" for x in sig], edgecolor="k", linewidth=0.3)
    ax.set_yticks(ypos)
    ax.set_yticklabels([r["species"].replace("_", " ") for r in per_species], fontsize=8)
    ax.axvline(float(perm[("subject", "strain")]["R2"]), color="k", linestyle="--", linewidth=1,
               label="whole-community strain R²")
    ax.set_xlabel("strain-level subject PERMANOVA R²")
    ax.set_title("C  Per-species individual signal")
    ax.grid(alpha=0.22, axis="x")
    ax.legend(loc="lower right", fontsize=8)

    ax = axes["D"]
    ax.boxplot([within, between], tick_labels=["Within subject", "Between subjects"], widths=0.55,
               patch_artist=True, boxprops={"facecolor": "#8ecae6", "edgecolor": "k"},
               medianprops={"color": "k"})
    for i, values in enumerate((within, between), 1):
        ax.scatter(np.full(len(values), i) + RNG.normal(0, 0.045, len(values)), values,
                   s=12, color="#1f77b4", alpha=0.40, linewidths=0)
    ax.set_ylabel("strain-profile Bray-Curtis distance")
    ax.set_title(f"D  Within- vs between-subject distance\np = {p_within:.1e}")
    ax.grid(alpha=0.22, axis="y")

    ax = axes["E"]
    bars = ax.bar(["Species", "Strain"], [acc_species, acc_s2b], width=0.48,
                  color=["#9467bd", "#1f77b4"], edgecolor="k")
    ax.axhline(1 / len(subjects), color="k", linestyle=":", linewidth=1.2, label="Chance = 12.5%")
    for bar, value in zip(bars, [acc_species, acc_s2b]):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.025, f"{value:.1%}",
                ha="center", fontweight="bold")
    ax.set_ylim(0, 1.10)
    ax.set_ylabel("Leave-one-timepoint-out accuracy")
    ax.set_title("E  Within-subject host classification\n8 subjects; one timepoint held out")
    ax.grid(alpha=0.22, axis="y")
    ax.legend(loc="upper left", fontsize=8)

    fig.suptitle("Native 2bRAD saliva profiles are individual-specific within one day (8 subjects × 4 timepoints)",
                 fontsize=14, fontweight="bold")
    out = PAPER / "figures/saliva_fig07"
    fig.savefig(out.with_suffix(".png"), dpi=180, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(f"wrote {out}.png/.pdf; LOTO species={acc_species:.4f}, strain={acc_s2b:.4f}")


if __name__ == "__main__":
    main()
