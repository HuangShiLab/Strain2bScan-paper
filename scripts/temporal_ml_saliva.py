#!/usr/bin/env python3
"""Fig 10, panels D-E: (1) temporal stability — an individual's saliva strain profile is more similar
across the day (4 timepoints) than between individuals; (2) ML host-ID — leave-one-timepoint-out
subject classification from strain vs species features (strain2bfunc reported 100% at strain level)."""
import os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
rows = [l.rstrip("\n").split("\t") for l in open(f"{WORK}/results/saliva_strain_long.tsv")][1:]
samples = sorted(set(r[0] for r in rows)); si = {s: i for i, s in enumerate(samples)}
subject   = {s: s.split("-")[1] for s in samples}
timepoint = {s: s.split("-")[0] for s in samples}

def build(level):
    key = (lambda r: f"{r[3]}|{r[4]}") if level == "strain" else (lambda r: r[3])
    feats = sorted(set(key(r) for r in rows)); fi = {f: j for j, f in enumerate(feats)}
    M = np.zeros((len(samples), len(feats)))
    for r in rows:
        M[si[r[0]], fi[key(r)]] += float(r[7])
    rs = M.sum(1, keepdims=True); rs[rs == 0] = 1
    return M / rs

def bc(M):
    n = M.shape[0]; D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            den = (M[i] + M[j]).sum(); D[i, j] = D[j, i] = np.abs(M[i] - M[j]).sum() / den if den else 0
    return D

Mst = build("strain"); Dst = bc(Mst)
subj = np.array([subject[s] for s in samples]); tp = np.array([timepoint[s] for s in samples])

# ---- (1) temporal stability: within-subject vs between-subject strain BC ----
within, between = [], []
for i in range(len(samples)):
    for j in range(i + 1, len(samples)):
        (within if subj[i] == subj[j] else between).append(Dst[i, j])
within, between = np.array(within), np.array(between)
from scipy.stats import mannwhitneyu
u, p_wb = mannwhitneyu(within, between, alternative="less")
print(f"strain BC within-subject mean={within.mean():.3f} vs between={between.mean():.3f} (p={p_wb:.1e})")

# ---- (2) leave-one-timepoint-out subject classification (1-NN Bray-Curtis) ----
def loto_accuracy(M):
    D = bc(M); tps = sorted(set(tp)); correct = tot = 0
    for held in tps:
        test = np.where(tp == held)[0]; train = np.where(tp != held)[0]
        for i in test:
            d = D[i, train]; nn = train[np.argmin(d)]
            correct += (subj[nn] == subj[i]); tot += 1
    return correct / tot

acc_strain = loto_accuracy(Mst)
acc_species = loto_accuracy(build("species"))
print(f"leave-one-timepoint-out subject accuracy: strain={acc_strain:.3f}, species={acc_species:.3f}  (chance={1/len(set(subj)):.3f})")

with open(f"{WORK}/results/saliva_temporal_ml.tsv", "w") as w:
    w.write("metric\tstrain\tspecies\n")
    w.write(f"within_subject_BC_mean\t{within.mean():.4f}\t-\n")
    w.write(f"between_subject_BC_mean\t{between.mean():.4f}\t-\n")
    w.write(f"within_vs_between_mannwhitney_p\t{p_wb:.3e}\t-\n")
    w.write(f"leave_one_timepoint_out_subject_accuracy\t{acc_strain:.4f}\t{acc_species:.4f}\n")
    w.write(f"n_subjects\t{len(set(subj))}\t-\n")

# ---- figure ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.6))
bp = ax1.boxplot([within, between], tick_labels=["within-subject\n(across day)", "between-subject"],
                 patch_artist=True, widths=0.55, showfliers=False)
for patch, c in zip(bp["boxes"], ["#2ca02c", "#d62728"]): patch.set_facecolor(c); patch.set_alpha(0.55)
for i, data in enumerate([within, between], 1):
    ax1.scatter(np.random.default_rng(2).normal(i, 0.05, len(data)), data, s=10, color="k", alpha=0.3, zorder=3)
ax1.set_ylabel("strain-level Bray–Curtis distance")
ax1.set_title(f"Temporal stability: within-individual < between\n(p={p_wb:.1e}) — profiles stable across the day", fontsize=10.5)
ax1.grid(alpha=0.25, axis="y")
ax2.bar([0, 1], [acc_species, acc_strain], color=["#9467bd", "#1f77b4"], width=0.6)
ax2.axhline(1/len(set(subj)), color="k", ls="--", lw=1)
ax2.text(0.5, 1/len(set(subj)) + 0.03, f"chance = {1/len(set(subj)):.0%}", ha="center", fontsize=8.5, color="#333")
ax2.set_xticks([0, 1]); ax2.set_xticklabels(["species", "strain"]); ax2.set_ylim(0, 1.08)
ax2.set_ylabel("subject-ID accuracy"); ax2.set_title("Leave-one-timepoint-out host ID\n(train 3 timepoints → predict the 4th)", fontsize=10.5)
for i, v in enumerate([acc_species, acc_strain]): ax2.text(i, v + 0.02, f"{v:.0%}", ha="center", fontsize=10, weight="bold")
ax2.grid(alpha=0.25, axis="y")
fig.suptitle("Saliva 2bRAD: strain profiles are individual-specific and temporally stable (8 subjects × 4 timepoints)", fontsize=11.5)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(f"{WORK}/results/saliva_temporal_ml.png", dpi=150)
fig.savefig(f"{WORK}/results/saliva_temporal_ml.pdf")
print("wrote results/saliva_temporal_ml.tsv + .png/.pdf")
