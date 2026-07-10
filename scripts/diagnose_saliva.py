#!/usr/bin/env python3
"""Diagnose the saliva null result: NaN/zero rows, subject-vs-timepoint grouping, per-species R^2."""
import numpy as np
WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
rng = np.random.default_rng(0)
rows = [l.rstrip("\n").split("\t") for l in open(f"{WORK}/results/saliva_strain_long.tsv")][1:]
samples = sorted(set(r[0] for r in rows))
subject = {r[0]: r[1] for r in rows}
timepoint = {r[0]: r[2] for r in rows}
si = {s: i for i, s in enumerate(samples)}

def matrix(feat_of, subset=None):
    feats = sorted(set(feat_of(r) for r in rows if subset is None or r[3] == subset))
    fi = {f: j for j, f in enumerate(feats)}
    M = np.zeros((len(samples), len(feats)))
    for r in rows:
        if subset is not None and r[3] != subset: continue
        M[si[r[0]], fi[feat_of(r)]] += float(r[7])
    rs = M.sum(1, keepdims=True); rs[rs == 0] = 1
    return M / rs

def bc(M):
    n = M.shape[0]; D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            den = (M[i] + M[j]).sum()
            D[i, j] = D[j, i] = np.abs(M[i] - M[j]).sum() / den if den else 0.0
    return D

def permanova(D, labels, perms=1999):
    labels = np.asarray(labels); N = len(labels); a = len(np.unique(labels))
    D2 = D ** 2
    def ss(lab):
        tot = D2[np.triu_indices(N, 1)].sum() / N; within = 0.0
        for g in np.unique(lab):
            idx = np.where(lab == g)[0]; ng = len(idx)
            if ng > 1: within += D2[np.ix_(idx, idx)][np.triu_indices(ng, 1)].sum() / ng
        return tot, within
    tot, within = ss(labels); between = tot - within
    F = (between / (a - 1)) / (within / (N - a)); R2 = between / tot
    ge = sum(1 for _ in range(perms) if (lambda t,w:(t-w)/(a-1)/(w/(N-a)))(*ss(rng.permutation(labels))) >= F)
    return R2, F, (ge + 1) / (perms + 1)

subj = [subject[s] for s in samples]; tp = [timepoint[s] for s in samples]
print("=== sanity ===")
Mst = matrix(lambda r: f"{r[3]}|{r[4]}")
print(f"strain matrix {Mst.shape}, zero-rows={int((Mst.sum(1)==0).sum())}, NaN={int(np.isnan(Mst).sum())}")
D = bc(Mst); print(f"BC NaN={int(np.isnan(D).sum())}, max={np.nanmax(D):.3f}")

print("\n=== whole-community PERMANOVA ===")
for name, M in [("species", matrix(lambda r: r[3])), ("strain", Mst)]:
    D = bc(M)
    for fac, lab in [("subject", subj), ("timepoint", tp)]:
        R2, F, p = permanova(D, lab)
        print(f"{name:8s} ~ {fac:9s}  R2={R2:.3f} F={F:.2f} p={p:.3f}")

print("\n=== per-species strain PERMANOVA ~ subject (which species carry individual signal?) ===")
species = sorted(set(r[3] for r in rows))
res = []
for sp in species:
    M = matrix(lambda r: r[4], subset=sp)
    present = int((M.sum(1) > 0).sum())
    if present < 12: continue
    D = bc(M); R2, F, p = permanova(D, subj, perms=999)
    res.append((R2, p, present, M.shape[1], sp))
for R2, p, pres, nf, sp in sorted(res, reverse=True):
    print(f"  {sp:32s} R2={R2:.3f} p={p:.3f}  (in {pres}/32 samples, {nf} clusters)")
