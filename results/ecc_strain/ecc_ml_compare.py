#!/usr/bin/env python3
"""Compare strain-level vs species-level features for ECC diagnosis."""
import csv, os
from collections import defaultdict
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score
from sklearn.preprocessing import StandardScaler

PRED_DIR = "/lustre1/g/aos_shihuang/Strain2bScan/work/ecc_strain/preds"
META = "/lustre1/g/aos_shihuang/Strain2b/data/saliva_data/ECC_saliva/2b/micro_result/metadata.txt"
SPECIES_TABLE = "/lustre1/g/aos_shihuang/Strain2b/data/saliva_data/ECC_saliva/2b/micro_result/feature_table_species.txt"

def read_meta(path):
    labels = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            labels[r["sample-id"]] = 1 if r["status"] == "C" else 0
    return labels

def read_preds(pred_dir):
    samples = defaultdict(list)
    for fname in os.listdir(pred_dir):
        if not fname.endswith(".pred.tsv"):
            continue
        sample = fname.replace(".pred.tsv", "")
        path = os.path.join(pred_dir, fname)
        with open(path) as f:
            header = None
            for ln in f:
                if not ln.strip():
                    continue
                if ln.startswith("#species"):
                    header = ln.lstrip("#").rstrip("\n").split("\t")
                    continue
                if ln.startswith("#"):
                    continue
                if header is None:
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) != len(header):
                    continue
                row = dict(zip(header, parts))
                species = row.get("species", "")
                cluster = row.get("cluster", "")
                sf = float(row.get("sample_fraction", "0") or 0)
                strain_id = f"{species}|{cluster}"
                samples[sample].append((strain_id, sf))
    return samples

def build_matrix(samples, min_prev=2):
    strain_prev = defaultdict(int)
    for s, calls in samples.items():
        for sid, _ in calls:
            strain_prev[sid] += 1
    keep = {sid for sid, n in strain_prev.items() if n >= min_prev}
    strains = sorted(keep)
    sample_list = sorted(samples.keys())
    si = {s: i for i, s in enumerate(sample_list)}
    fi = {sid: i for i, sid in enumerate(strains)}
    M = np.zeros((len(sample_list), len(strains)))
    for s, calls in samples.items():
        for sid, sf in calls:
            if sid in fi:
                M[si[s], fi[sid]] += sf
    rs = M.sum(1, keepdims=True)
    rs[rs == 0] = 1
    return sample_list, strains, M / rs

def read_species_matrix(path):
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    species = [r["sample-id"] for r in rows]
    samples = [c for c in rows[0].keys() if c != "sample-id"]
    M = np.zeros((len(samples), len(species)))
    for i, sp in enumerate(species):
        for j, smp in enumerate(samples):
            M[j, i] = float(rows[i][smp])
    rs = M.sum(1, keepdims=True)
    rs[rs == 0] = 1
    return samples, species, M / rs

def evaluate(X, y, name, cv):
    y = np.asarray(y)
    results = []
    for clf_name, clf in [("RF", RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1)),
                          ("LR", LogisticRegression(max_iter=1000, class_weight="balanced"))]:
        if clf_name == "LR":
            Xs = StandardScaler().fit_transform(X)
        else:
            Xs = X
        y_prob = cross_val_predict(clf, Xs, y, cv=cv, method="predict_proba")[:, 1]
        y_pred = (y_prob >= 0.5).astype(int)
        acc = accuracy_score(y, y_pred)
        auc = roc_auc_score(y, y_prob)
        f1 = f1_score(y, y_pred)
        results.append((name, clf_name, acc, auc, f1))
        print(f"  {name} {clf_name}: accuracy={acc:.3f} AUROC={auc:.3f} F1={f1:.3f}")
    return results

def main():
    labels = read_meta(META)
    print(f"Metadata: {sum(labels.values())} Caries / {len(labels)-sum(labels.values())} Healthy")

    print("\n=== Strain-level features ===")
    samples_strain = read_preds(PRED_DIR)
    sample_list, strains, M_strain = build_matrix(samples_strain, min_prev=2)
    print(f"Samples: {len(sample_list)}, strains: {len(strains)}")
    y = [labels[s] for s in sample_list]
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    res_strain = evaluate(M_strain, y, "strain", cv)

    print("\n=== Species-level features ===")
    sample_list_sp, species, M_sp = read_species_matrix(SPECIES_TABLE)
    common = [s for s in sample_list if s in sample_list_sp]
    idx_strain = [sample_list.index(s) for s in common]
    idx_sp = [sample_list_sp.index(s) for s in common]
    M_strain_common = M_strain[idx_strain]
    M_sp_common = M_sp[idx_sp]
    y_common = [labels[s] for s in common]
    print(f"Common samples: {len(common)}, species: {len(species)}")
    res_sp = evaluate(M_sp_common, y_common, "species", cv)
    res_strain_common = evaluate(M_strain_common, y_common, "strain_common", cv)

    print("\n=== Summary (5-fold CV) ===")
    print("feature\tclf\taccuracy\tAUROC\tF1")
    for r in res_sp + res_strain_common:
        print(f"{r[0]}\t{r[1]}\t{r[2]:.3f}\t\t{r[3]:.3f}\t{r[4]:.3f}")

if __name__ == "__main__":
    main()
