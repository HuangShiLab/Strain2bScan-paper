#!/usr/bin/env python3
"""Fig 7 (bottom row): per-species individual signal + ML host-ID.
* left: per-species strain-level subject discrimination (moved from the old top-right panel)
* right: leave-one-timepoint-out subject classification accuracy
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial"]
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import mannwhitneyu

PAPER = Path(__file__).resolve().parent.parent
WORK = PAPER / "results"
OUT_FIG = PAPER / "figures"

rows = [l.rstrip("\n").split("\t") for l in open(WORK / "saliva_strain_long.tsv")][1:]
samples = sorted(set(r[0] for r in rows))
si = {s: i for i, s in enumerate(samples)}
subject = {s: s.split("-")[1] for s in samples}
timepoint = {s: s.split("-")[0] for s in samples}


def build(level):
    key = (lambda r: f"{r[3]}|{r[4]}") if level == "strain" else (lambda r: r[3])
    feats = sorted(set(key(r) for r in rows))
    fi = {f: j for j, f in enumerate(feats)}
    M = np.zeros((len(samples), len(feats)))
    for r in rows:
        M[si[r[0]], fi[key(r)]] += float(r[7])
    rs = M.sum(1, keepdims=True)
    rs[rs == 0] = 1
    return M / rs


def bc(M):
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            den = (M[i] + M[j]).sum()
            D[i, j] = D[j, i] = np.abs(M[i] - M[j]).sum() / den if den else 0
    return D


Mst = build("strain")
Dst = bc(Mst)
subj = np.array([subject[s] for s in samples])
tp = np.array([timepoint[s] for s in samples])

# ---- temporal stability: within-subject vs between-subject strain BC ----
within, between = [], []
for i in range(len(samples)):
    for j in range(i + 1, len(samples)):
        (within if subj[i] == subj[j] else between).append(Dst[i, j])
within, between = np.array(within), np.array(between)
u, p_wb = mannwhitneyu(within, between, alternative="less")
print(f"strain BC within-subject mean={within.mean():.3f} vs between={between.mean():.3f} (p={p_wb:.1e})")


def loto_accuracy(M):
    D = bc(M)
    tps = sorted(set(tp))
    correct = tot = 0
    for held in tps:
        test = np.where(tp == held)[0]
        train = np.where(tp != held)[0]
        for i in test:
            d = D[i, train]
            nn = train[np.argmin(d)]
            correct += (subj[nn] == subj[i])
            tot += 1
    return correct / tot


acc_strain = loto_accuracy(Mst)
acc_species = loto_accuracy(build("species"))
print(f"leave-one-timepoint-out subject accuracy: strain={acc_strain:.3f}, species={acc_species:.3f}  "
      f"(chance={1/len(set(subj)):.3f})")

with open(WORK / "saliva_temporal_ml.tsv", "w") as w:
    w.write("metric\tstrain\tspecies\n")
    w.write(f"within_subject_BC_mean\t{within.mean():.4f}\t-\n")
    w.write(f"between_subject_BC_mean\t{between.mean():.4f}\t-\n")
    w.write(f"within_vs_between_mannwhitney_p\t{p_wb:.3e}\t-\n")
    w.write(f"leave_one_timepoint_out_subject_accuracy\t{acc_strain:.4f}\t{acc_species:.4f}\n")
    w.write(f"n_subjects\t{len(set(subj))}\t-\n")

# ---- read per-species discrimination and whole-community strain R^2 ----
perspec = []
with open(WORK / "saliva_perspecies_subject.tsv") as f:
    header = next(f)
    for line in f:
        sp, r2, p, pres, ncl = line.rstrip("\n").split("\t")
        perspec.append((sp, float(r2), float(p), int(pres), int(ncl)))

strain_community_r2 = None
with open(WORK / "saliva_permanova.tsv") as f:
    header = next(f)
    for line in f:
        fields = line.rstrip("\n").split("\t")
        if fields[0] == "subject" and fields[1] == "strain":
            strain_community_r2 = float(fields[4])
            break

# ---- figure: bottom row of Fig 7 ----
fig, (ax_species, ax_ml) = plt.subplots(
    1, 2, figsize=(10.5, 4.6),
    gridspec_kw={"width_ratios": [2.2, 1.2]}
)

# Left: per-species strain-level subject R^2 (horizontal bars)
names = [x[0].replace("_", " ") for x in perspec][::-1]
r2s = [x[1] for x in perspec][::-1]
sig = [x[2] < 0.05 for x in perspec][::-1]
ypos = range(len(names))
bar_colors = ["#2ca02c" if "Rothia mucilaginosa" in n else ("#1f77b4" if s else "#bbbbbb")
              for n, s in zip(names, sig)]
ax_species.barh(list(ypos), r2s, color=bar_colors, edgecolor="k", linewidth=0.3)
ax_species.set_yticks(list(ypos))
ax_species.set_yticklabels(names, fontsize=10.5)
ax_species.set_xlabel("strain-level PERMANOVA R² (subject)", fontsize=11)
ax_species.set_xlim(0, 1)
ax_species.tick_params(labelsize=10)
ax_species.axvline(strain_community_r2, color="k", ls="--", lw=1,
                   label=f"whole-community strain R²={strain_community_r2:.2f}")
ax_species.set_title("Per-species strain-level individual signal\n(green = R. mucilaginosa; grey = ns)",
                     fontsize=12)
ax_species.grid(alpha=0.25, axis="x")
ax_species.legend(fontsize=9, loc="lower right")

# Right: leave-one-timepoint-out accuracy (narrower bars, closer spacing)
bar_x = [0, 0.55]
ax_ml.bar(bar_x, [acc_species, acc_strain], color=["#9467bd", "#1f77b4"], width=0.30)
ax_ml.set_xticks(bar_x)
ax_ml.set_xticklabels(["species", "strain"], fontsize=10)
ax_ml.set_xlim(-0.35, 0.9)
ax_ml.set_ylim(0, 1.08)
ax_ml.set_ylabel("subject-ID accuracy", fontsize=11)
ax_ml.tick_params(labelsize=10)
ax_ml.set_title("Leave-one-timepoint-out\nhost-ID accuracy\n(train 3, predict 4th)",
                fontsize=11.5)
for i, v in zip(bar_x, [acc_species, acc_strain]):
    ax_ml.text(i, v + 0.02, f"{v:.0%}", ha="center", fontsize=12, weight="bold")
ax_ml.grid(alpha=0.25, axis="y")

fig.suptitle("Saliva 2bRAD: strain profiles are individual-specific and temporally stable (8 subjects × 4 timepoints)",
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT_FIG / "saliva_temporal_ml.png", dpi=150)
fig.savefig(OUT_FIG / "saliva_temporal_ml.pdf")
print(f"wrote {OUT_FIG / 'saliva_temporal_ml.png/pdf'}")
