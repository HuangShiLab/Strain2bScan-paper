#!/usr/bin/env python3
"""Final report: compare Strain2bScan strain-level vs species-level on HPC saliva data.

Inputs (expected in work/saliva_hpc_results/):
  saliva_strain_long_2brad.tsv
  saliva_strain_long_wms.tsv
  Abundance_Stat.all.xls          (2bRAD species-level)
  merged_abundance_table_species.txt  (WMS species-level, MetaPhlAn)

Outputs:
  work/saliva_hpc_results/saliva_discrimination_summary.tsv
  figures/saliva_hpc_pcoa.png/pdf
"""
import os, sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
WORK = PAPER / "work" / "saliva_hpc_results"
FIGDIR = PAPER / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)


def alias_to_hpc(a):
    tp, subj = a.split("-")
    return f"s{int(tp[1:])}_{int(subj)}"


def read_strain(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def read_species_2brad(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        sample_cols = [h for h in header if h.startswith("s") and "_" in h and h[1].isdigit()]
        for line in f:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows, sample_cols


def read_species_wms(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        sample_cols = [h for h in header if h.startswith("S") and "-" in h]
        for line in f:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows, sample_cols


def braycurtis(M):
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.abs(M[i] - M[j]).sum()
            den = (M[i] + M[j]).sum()
            D[i, j] = D[j, i] = num / den if den else 0.0
    return D


def permanova(D, labels, perms=4999):
    labels = np.asarray(labels)
    N = len(labels)
    groups = np.unique(labels)
    a = len(groups)

    def ss(D2, lab):
        tot = D2[np.triu_indices(N, 1)].sum() / N
        within = 0.0
        for g in np.unique(lab):
            idx = np.where(lab == g)[0]
            ng = len(idx)
            if ng > 1:
                within += D2[np.ix_(idx, idx)][np.triu_indices(ng, 1)].sum() / ng
        return tot, within

    D2 = D ** 2
    tot, within = ss(D2, labels)
    between = tot - within
    F = (between / (a - 1)) / (within / (N - a))
    R2 = between / tot
    ge = 1
    for _ in range(perms):
        p = rng.permutation(labels)
        t2, w2 = ss(D2, p)
        b2 = t2 - w2
        Fp = (b2 / (a - 1)) / (w2 / (N - a))
        if Fp >= F:
            ge += 1
    return R2, F, ge / (perms + 1)


def loo_1nn(D, labels):
    labels = np.asarray(labels)
    n = len(labels)
    correct = 0
    for i in range(n):
        d = D[i].copy()
        d[i] = np.inf
        j = int(np.argmin(d))
        correct += (labels[j] == labels[i])
    return correct / n


def pcoa(D):
    n = D.shape[0]
    A = -0.5 * D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = J @ A @ J
    w, V = np.linalg.eigh(B)
    order = np.argsort(w)[::-1]
    w = w[order]
    V = V[:, order]
    pos = w > 0
    coords = V[:, :2] * np.sqrt(np.abs(w[:2]))
    ev = w[pos] / w[pos].sum()
    return coords, ev


def compare(source, strain_path, species_path, species_reader):
    source = source.lower()
    strain_rows = read_strain(strain_path)
    spec_rows, sample_cols = species_reader(species_path)

    if source == "2brad":
        spec_by_sp = {r["Species"]: r for r in spec_rows}
        strain_samples = sorted(set(r["sample"] for r in strain_rows))
        spec_samples = [f"S{int(c.split('_')[0][1:])}-{int(c.split('_')[1])}" for c in sample_cols]
    else:
        spec_by_sp = {r["sample-id"]: r for r in spec_rows}
        strain_samples = sorted(set(r["sample"] for r in strain_rows))
        spec_samples = sample_cols

    common = sorted(set(strain_samples) & set(spec_samples))
    print(f"{source}: {len(common)} common samples")
    if not common:
        return []

    # species-level matrix
    sp_names = sorted(spec_by_sp.keys())
    M_spec = np.zeros((len(common), len(sp_names)))
    for i, alias in enumerate(common):
        if source == "2brad":
            col = alias_to_hpc(alias)
        else:
            col = alias
        for j, sp in enumerate(sp_names):
            val = spec_by_sp[sp].get(col, "0")
            M_spec[i, j] = float(val)
    rs = M_spec.sum(1, keepdims=True)
    rs[rs == 0] = 1
    M_spec = M_spec / rs

    # strain-level species-aggregated matrix
    strain_sp = sorted(set(r["species"] for r in strain_rows if r["sample"] in common))
    M_strain_sp = np.zeros((len(common), len(strain_sp)))
    si = {s: i for i, s in enumerate(common)}
    fi = {s: i for i, s in enumerate(strain_sp)}
    for r in strain_rows:
        if r["sample"] not in si:
            continue
        M_strain_sp[si[r["sample"]], fi[r["species"]]] += float(r["sample_fraction"])
    rs = M_strain_sp.sum(1, keepdims=True)
    rs[rs == 0] = 1
    M_strain_sp = M_strain_sp / rs

    # strain-level full matrix
    strain_feats = sorted(set(f"{r['species']}|{r['cluster']}" for r in strain_rows if r["sample"] in common))
    M_strain = np.zeros((len(common), len(strain_feats)))
    fi2 = {f: i for i, f in enumerate(strain_feats)}
    for r in strain_rows:
        if r["sample"] not in si:
            continue
        M_strain[si[r["sample"]], fi2[f"{r['species']}|{r['cluster']}"]] += float(r["sample_fraction"])
    rs = M_strain.sum(1, keepdims=True)
    rs[rs == 0] = 1
    M_strain = M_strain / rs

    subjects = [s.split("-")[1] for s in common]
    timepts = [s.split("-")[0] for s in common]

    out = []
    for level, M in [("species_table", M_spec), ("strain_aggregated", M_strain_sp), ("strain_full", M_strain)]:
        D = braycurtis(M)
        for factor, labels in [("subject", subjects), ("timepoint", timepts)]:
            R2, F, p = permanova(D, labels)
            acc = loo_1nn(D, labels)
            out.append({
                "data_type": source,
                "level": level,
                "factor": factor,
                "n_samples": len(common),
                "n_features": M.shape[1],
                "R2": f"{R2:.4f}",
                "pseudo_F": f"{F:.3f}",
                "p_value": f"{p:.4f}",
                "LOO_1NN_acc": f"{acc:.3f}",
            })
    return out, common, M_strain


def main():
    results = []

    brad_out, brad_samples, brad_M = compare(
        "2bRAD",
        WORK / "saliva_strain_long_2brad.tsv",
        WORK / "Abundance_Stat.all.xls",
        read_species_2brad,
    )
    results.extend(brad_out)

    wms_strain_path = WORK / "saliva_strain_long_wms.tsv"
    if wms_strain_path.exists():
        wms_out, wms_samples, wms_M = compare(
            "WMS",
            wms_strain_path,
            WORK / "merged_abundance_table_species.txt",
            read_species_wms,
        )
        results.extend(wms_out)
    else:
        print(f"WMS strain results not found at {wms_strain_path}; skipping")
        wms_samples, wms_M = None, None

    out_tsv = WORK / "saliva_discrimination_summary.tsv"
    with open(out_tsv, "w") as w:
        w.write("data_type\tlevel\tfactor\tn_samples\tn_features\tR2\tpseudo_F\tp_value\tLOO_1NN_acc\n")
        for r in results:
            w.write("\t".join(str(r[c]) for c in [
                "data_type", "level", "factor", "n_samples", "n_features",
                "R2", "pseudo_F", "p_value", "LOO_1NN_acc"
            ]) + "\n")
    print(f"wrote {out_tsv}")

    # PCoA figure for both data types (strain full)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, (source, samples, M) in [(axes[0], ("2bRAD", brad_samples, brad_M)),
                                      (axes[1], ("WMS", wms_samples, wms_M))]:
        if samples is None:
            ax.set_visible(False)
            continue
        subjects = [s.split("-")[1] for s in samples]
        timepts = [s.split("-")[0] for s in samples]
        D = braycurtis(M)
        coords, ev = pcoa(D)
        subs = sorted(set(subjects), key=int)
        palette = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b",
                   "#17becf", "#e377c2", "#bcbd22", "#7f7f7f"]
        colors = {s: palette[i % len(palette)] for i, s in enumerate(subs)}
        for s in subs:
            idx = [i for i, sm in enumerate(samples) if subjects[i] == s]
            ax.scatter(coords[idx, 0], coords[idx, 1], s=70, color=colors[s],
                       label=f"subject {s}", edgecolor="k", linewidth=0.5, alpha=0.9)
        R2, F, p = permanova(D, subjects)
        ax.set_title(f"{source} strain-level (R²={R2:.2f}, p={p:.3f})")
        ax.set_xlabel(f"PCo1 ({ev[0]*100:.0f}%)")
        ax.set_ylabel(f"PCo2 ({ev[1]*100:.0f}%)")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc="best", ncol=2)
    fig.tight_layout()
    fig.savefig(FIGDIR / "saliva_hpc_pcoa.png", dpi=150)
    fig.savefig(FIGDIR / "saliva_hpc_pcoa.pdf")
    print(f"wrote {FIGDIR / 'saliva_hpc_pcoa.png/pdf'}")


if __name__ == "__main__":
    main()
