#!/usr/bin/env python3
"""Score completed StrainScan rerun outputs against the committed simulation truth.

Cluster IDs are tool-specific. Truth genomes are mapped through each StrainScan database's
`Cluster_Result/hclsMap_95_recls.txt`; final-report strain IDs are mapped through the same file.
"""
import csv, re, sys
from collections import defaultdict
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
WORK = PAPER / "work/strainscan_headtohead_rerun"

def read_tsv(path):
    with open(path) as fh:
        return [line.rstrip("\n").split("\t") for line in fh if not line.startswith("#")]

def cluster_map(db):
    genome_to_cluster = {}
    with open(db / "Cluster_Result/hclsMap_95_recls.txt") as fh:
        for line in fh:
            cluster, n, members = line.rstrip("\n").split("\t")
            for genome in members.split(","):
                if genome:
                    genome_to_cluster[genome] = f"C{cluster}"
    return genome_to_cluster

def report_clusters(path):
    if not path.exists():
        return set()
    clusters = set()
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("Cluster_ID"):
                clusters.add(row["Cluster_ID"])
    return clusters

def metrics(truth, pred):
    tp = len(truth & pred); fp = len(pred - truth); fn = len(truth - pred)
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else float("nan")
    return tp, fp, fn, precision, recall, f1

def depth_of(sample):
    m = re.search(r"_d([0-9.]+)$", sample)
    return m.group(1) if m else ""

def main():
    rows = []
    manifest = list(csv.DictReader(open(WORK / "run_manifest.tsv"), delimiter="\t"))
    maps = {}
    for row in manifest:
        sp = row["species"]
        if sp not in maps:
            db = WORK / "db" / sp
            maps[sp] = cluster_map(db) if (db / "Cluster_Result/hclsMap_95_recls.txt").exists() else {}
        g2c = maps[sp]
        truth_path = Path(row["truth"])
        if row["kind"] == "single":
            truth_rows = read_tsv(truth_path)
            truth = {g2c[r[1]] for r in truth_rows if r[1] in g2c}
            pred = report_clusters(WORK / "results/single" / row["sample"] / row["species"] / "final_report.txt")
        else:
            truth = set(); pred = set()
            for r in read_tsv(truth_path):
                sp = r[0]; genome = r[1]
                if sp not in maps:
                    db = WORK / "db" / sp
                    maps[sp] = cluster_map(db) if (db / "Cluster_Result/hclsMap_95_recls.txt").exists() else {}
                if genome in maps[sp]: truth.add(f"{sp}|{maps[sp][genome]}")
            species_dir = WORK / "results/multi" / row["sample"]
            for report in species_dir.glob("*/final_report.txt"):
                sp = report.parent.name
                for cluster in report_clusters(report):
                    pred.add(f"{sp}|{cluster}")
        tp, fp, fn, precision, recall, f1 = metrics(truth, pred)
        rows.append({"kind": row["kind"], "sample": row["sample"], "species": row["species"],
                     "depth": depth_of(row["sample"]), "n_truth": len(truth), "n_pred": len(pred),
                     "tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1})
    out = WORK / "scored_runs.tsv"
    with open(out, "w", newline="") as fh:
        fields = ["kind","sample","species","depth","n_truth","n_pred","tp","fp","fn","precision","recall","f1"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    print(f"wrote {out} ({len(rows)} scored runs)")
    complete = [r for r in rows if r["kind"] == "single"]
    print(f"single completed: {len(complete)}")

if __name__ == "__main__":
    main()
