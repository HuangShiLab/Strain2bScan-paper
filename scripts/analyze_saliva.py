#!/usr/bin/env python3
"""Fig 10 - saliva individual discrimination. Build strain-level and species-level abundance
matrices from Strain2bScan calls on the 32 diurnal saliva 2bRAD samples; Bray-Curtis + PERMANOVA
(adonis, subject as factor). Expect strain-level R^2 > species-level R^2."""
import os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
rng = np.random.default_rng(42)

rows = [l.rstrip("\n").split("\t") for l in open(f"{WORK}/results/saliva_strain_long.tsv")][1:]
samples = sorted(set(r[0] for r in rows))
# alias S<timepoint>-<subject>: prefix (S9/S11/S13/S17) = time of day; suffix (1..14) = subject.
subject   = {s: s.split("-")[1] for s in samples}
timepoint = {s: s.split("-")[0] for s in samples}
# feature keys
strain_feats  = sorted(set(f"{r[3]}|{r[4]}" for r in rows))
species_feats = sorted(set(r[3] for r in rows))
si = {s: i for i, s in enumerate(samples)}

def build(level):
    feats = strain_feats if level == "strain" else species_feats
    fi = {f: j for j, f in enumerate(feats)}
    M = np.zeros((len(samples), len(feats)))
    for r in rows:
        sup = float(r[7])
        key = f"{r[3]}|{r[4]}" if level == "strain" else r[3]
        M[si[r[0]], fi[key]] += sup
    # row-normalize to relative abundance
    rs = M.sum(1, keepdims=True); rs[rs == 0] = 1
    return M / rs

def braycurtis(M):
    n = M.shape[0]; D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.abs(M[i] - M[j]).sum(); den = (M[i] + M[j]).sum()
            D[i, j] = D[j, i] = num / den if den else 0.0
    return D

def permanova(D, labels, perms=4999):
    labels = np.asarray(labels); N = len(labels); groups = np.unique(labels)
    a = len(groups)
    def ss(D2, lab):
        tot = (D2[np.triu_indices(N, 1)]).sum() / N
        within = 0.0
        for g in np.unique(lab):
            idx = np.where(lab == g)[0]; ng = len(idx)
            if ng > 1:
                sub = D2[np.ix_(idx, idx)][np.triu_indices(ng, 1)]
                within += sub.sum() / ng
        return tot, within
    D2 = D ** 2
    tot, within = ss(D2, labels); between = tot - within
    F = (between / (a - 1)) / (within / (N - a))
    R2 = between / tot
    ge = 1
    for _ in range(perms):
        p = rng.permutation(labels)
        t2, w2 = ss(D2, p); b2 = t2 - w2
        Fp = (b2 / (a - 1)) / (w2 / (N - a))
        if Fp >= F: ge += 1
    return R2, F, ge / (perms + 1)

def loo_1nn(D, labels):
    """Leave-one-out 1-NN subject classification accuracy on a distance matrix."""
    labels = np.asarray(labels); n = len(labels); correct = 0
    for i in range(n):
        d = D[i].copy(); d[i] = np.inf
        j = int(np.argmin(d))
        correct += (labels[j] == labels[i])
    return correct / n

out = ["factor\tlevel\tn_samples\tn_features\tR2\tpseudo_F\tp_value\tLOO_1NN_acc"]
res = {}
labels = [subject[s] for s in samples]
tlabels = [timepoint[s] for s in samples]
for level in ("species", "strain"):
    M = build(level); D = braycurtis(M)
    R2, F, p = permanova(D, labels); acc = loo_1nn(D, labels)
    res[level] = (M, D, R2, F, p, acc)
    out.append(f"subject\t{level}\t{len(samples)}\t{M.shape[1]}\t{R2:.4f}\t{F:.3f}\t{p:.4f}\t{acc:.3f}")
    print(out[-1])
    # timepoint (time-of-day) as a control factor — expect near-zero
    tR2, tF, tp_ = permanova(D, tlabels); tacc = loo_1nn(D, tlabels)
    out.append(f"timepoint\t{level}\t{len(samples)}\t{M.shape[1]}\t{tR2:.4f}\t{tF:.3f}\t{tp_:.4f}\t{tacc:.3f}")
open(f"{WORK}/results/saliva_permanova.tsv", "w").write("\n".join(out) + "\n")

# ---- per-species strain-level subject discrimination (mirrors strain2bfunc's Rothia result) ----
def species_matrix(sp):
    feats = sorted(set(r[4] for r in rows if r[3] == sp)); fi = {f: j for j, f in enumerate(feats)}
    M = np.zeros((len(samples), len(feats)))
    for r in rows:
        if r[3] == sp: M[si[r[0]], fi[r[4]]] += float(r[7])
    rs = M.sum(1, keepdims=True); rs[rs == 0] = 1
    return M / rs
perspec = []
for sp in sorted(set(r[3] for r in rows)):
    M = species_matrix(sp); pres = int((M.sum(1) > 0).sum())
    if pres < 12: continue
    R2, F, p = permanova(braycurtis(M), labels)
    perspec.append((R2, p, pres, M.shape[1], sp))
perspec.sort(reverse=True)
with open(f"{WORK}/results/saliva_perspecies_subject.tsv", "w") as w:
    w.write("species\tstrain_R2_subject\tp_value\tn_samples_present\tn_clusters\n")
    for R2, p, pres, nf, sp in perspec:
        w.write(f"{sp}\t{R2:.4f}\t{p:.4f}\t{pres}\t{nf}\n")
print(f"per-species: top = {perspec[0][4]} R2={perspec[0][0]:.3f}; {sum(1 for x in perspec if x[1]<0.05)}/{len(perspec)} species significant")

# ---- PCoA figure (classical MDS on Bray-Curtis), both levels ----
def pcoa(D):
    n = D.shape[0]; A = -0.5 * D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = J @ A @ J
    w, V = np.linalg.eigh(B)
    order = np.argsort(w)[::-1]; w = w[order]; V = V[:, order]
    pos = w > 0
    coords = V[:, :2] * np.sqrt(np.abs(w[:2]))
    ev = w[pos] / w[pos].sum()
    return coords, ev

subs = sorted(set(subject.values()), key=int)
palette = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b", "#17becf", "#e377c2",
           "#bcbd22", "#7f7f7f"]
colors = {s: palette[i % len(palette)] for i, s in enumerate(subs)}
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
for ax, level in zip(axes[:2], ("species", "strain")):
    M, D, R2, F, p, acc = res[level]
    coords, ev = pcoa(D)
    for s in subs:
        idx = [i for i, sm in enumerate(samples) if subject[sm] == s]
        ax.scatter(coords[idx, 0], coords[idx, 1], s=70, color=colors[s], label=f"subject {s}",
                   edgecolor="k", linewidth=0.5, alpha=0.9)
    ax.set_title(f"{level.capitalize()}-level  (R²={R2:.2f}, p={p:.3f}, 1-NN acc={acc:.0%})", fontsize=11)
    ax.set_xlabel(f"PCo1 ({ev[0]*100:.0f}%)"); ax.set_ylabel(f"PCo2 ({ev[1]*100:.0f}%)")
    ax.grid(alpha=0.25)
axes[0].legend(fontsize=8, loc="best", ncol=2)
# Panel C: per-species strain-level subject R^2 (horizontal bars, Rothia highlighted)
axc = axes[2]
names = [x[4].replace("_", " ") for x in perspec][::-1]
r2s   = [x[0] for x in perspec][::-1]
sig   = [x[1] < 0.05 for x in perspec][::-1]
ypos = range(len(names))
bar_colors = ["#2ca02c" if "Rothia mucilaginosa" in n else ("#1f77b4" if s else "#bbbbbb")
              for n, s in zip(names, sig)]
axc.barh(list(ypos), r2s, color=bar_colors, edgecolor="k", linewidth=0.3)
axc.set_yticks(list(ypos)); axc.set_yticklabels(names, fontsize=7.5)
axc.set_xlabel("strain-level PERMANOVA R² (subject)"); axc.set_xlim(0, 1)
axc.axvline(res["species"][2], color="k", ls="--", lw=1, label=f"whole-community strain R²={res['strain'][2]:.2f}")
axc.set_title("Per-species strain-level individual signal\n(green = R. mucilaginosa; grey = ns)", fontsize=10)
axc.grid(alpha=0.25, axis="x")
fig.suptitle("Saliva 2bRAD (Strain2bScan): strain-level profiles discriminate individuals — 8 subjects × 4 timepoints",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f"{WORK}/results/saliva_individual_discrimination.png", dpi=150)
fig.savefig(f"{WORK}/results/saliva_individual_discrimination.pdf")
print(f"\nwrote saliva_permanova.tsv + saliva_individual_discrimination.png/pdf")
print(f"species R²={res['species'][2]:.3f}  strain R²={res['strain'][2]:.3f}")
