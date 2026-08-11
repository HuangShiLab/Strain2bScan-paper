#!/usr/bin/env python3
"""Compute StrainScan precision/recall/F1 for the sim benchmark.

Reads StrainScan final_report.txt files from the HPC output tree and writes
figure_raw_data/sim_headtohead/strainscan_single_persample.tsv and
strainscan_multi_persample.tsv with the same columns as the previous hand-run
versions, so existing plotting scripts can consume them unchanged.

Metrics:
  - Single-species: truth cluster IDs are read from the .truth.tsv files.
    A predicted strain is a TP if its source genome maps to a truth cluster.
    TP = # truth clusters hit; FP = # predicted strains not mapping to a truth
    cluster; FN = # truth clusters missed.
  - Multi-species: truth species are read from the .truth.tsv files.
    Predicted strains are mapped to species via results/genome_to_species.tsv.
    TP = # true species hit; FP = # predicted species not in truth;
    FN = # true species missed.
"""
import csv
import re
from pathlib import Path
from collections import defaultdict

PAPER = Path(__file__).resolve().parent.parent
RES_ROOT = Path("/lustre1/g/aos_shihuang/LU/Strain2bScan-raw-data/sim_benchmark/strainscan_results")
LOCAL_FALLBACK = PAPER / "work" / "mock_retest" / "Strain2bScan-raw-data" / "sim_benchmark" / "results" / "strainscan"
FIGRAW = PAPER / "figure_raw_data" / "sim_headtohead"
TRUTH_SINGLE = PAPER / "figure_raw_data" / "sim_single_species"
TRUTH_MULTI = PAPER / "figure_raw_data" / "sim_multi_species"

SINGLE_SAMPLES = [
    ("Escherichia_coli", "Escherichia_coli__diff_k2_rep1_d5", "diff", 2, 1, 5),
    ("Cutibacterium_acnes", "Cutibacterium_acnes__diff_k2_rep1_d5", "diff", 2, 1, 5),
    ("Staphylococcus_epidermidis", "Staphylococcus_epidermidis__diff_k2_rep1_d5", "diff", 2, 1, 5),
    ("Prevotella_copri", "Prevotella_copri__diff_k2_rep1_d5", "diff", 2, 1, 5),
]
MULTI_DEPTHS = ["depth_low", "depth_med", "depth_high"]

# genome_to_species.tsv uses current taxonomic names; the simulation truth files
# use the names from the time the genomes were downloaded. Map synonyms so
# metrics are comparable.
SPECIES_SYNONYMS = {
    "Propionibacterium acnes": "Cutibacterium acnes",
    "Bacteroides dorei": "Phocaeicola dorei",
    "Segatella copri": "Prevotella copri",
    "Peptoclostridium difficile": "Clostridioides difficile",
    "Lactobacillus plantarum": "Lactiplantibacillus plantarum",
}


def normalize_species(name):
    name = name.strip()
    name = SPECIES_SYNONYMS.get(name, name)
    return name.replace(" ", "_")


def load_genome_to_species():
    out = {}
    path = PAPER / "results" / "genome_to_species.tsv"
    if path.exists():
        for row in csv.DictReader(open(path), delimiter="\t"):
            out[row["genome"]] = normalize_species(row["species"])
    return out


def read_final_report(path):
    rows = []
    if not path.exists():
        return rows
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def get_predicted_strains(res_root, sample_key, species):
    """Return set of strain names predicted for one sample/species pair."""
    out = set()
    report_dir = res_root / sample_key / species
    for report_path in report_dir.rglob("final_report.txt"):
        for row in read_final_report(report_path):
            name = row.get("Strain_Name", "").strip()
            if name:
                out.add(name)
    return out


def load_truth_single(sample_id):
    species = sample_id.split("__")[0]
    path = TRUTH_SINGLE / species / "truth" / f"{sample_id}.truth.tsv"
    clusters = set()
    genome_to_cluster = {}
    with open(path) as fh:
        next(fh)
        for ln in fh:
            p = ln.rstrip().split("\t")
            if len(p) >= 3:
                genome_to_cluster[p[0]] = p[1]
                clusters.add(p[1])
    return clusters, genome_to_cluster


def load_truth_multi(depth):
    path = TRUTH_MULTI / depth / "truth" / "sample01.truth.tsv"
    species = set()
    with open(path) as fh:
        next(fh)
        for ln in fh:
            p = ln.rstrip().split("\t")
            if len(p) >= 4:
                species.add(p[0])
    return species


def compute_single(res_root):
    rows = []
    for species, sample, strategy, k, rep, depth in SINGLE_SAMPLES:
        truth_clusters, genome_to_cluster = load_truth_single(sample)
        pred_strains = get_predicted_strains(res_root, f"single/{sample}", species)

        hit_clusters = set()
        fp = 0
        for strain in pred_strains:
            cluster = genome_to_cluster.get(strain)
            if cluster is not None:
                hit_clusters.add(cluster)
            else:
                # Try stripping version suffix
                base = re.split(r"[._]\d+$", strain)[0]
                cluster = genome_to_cluster.get(base)
                if cluster is not None:
                    hit_clusters.add(cluster)
                else:
                    fp += 1

        tp = len(hit_clusters)
        fn = len(truth_clusters - hit_clusters)
        n_truth = len(truth_clusters)
        n_pred = len(pred_strains)
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

        rows.append({
            "species": species,
            "sample": sample,
            "strategy": strategy,
            "k": k,
            "rep": rep,
            "depth": depth,
            "n_truth_cl": n_truth,
            "n_pred_cl": n_pred,
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "wall_s": "",
            "maxrss_mb": "",
        })
    return rows


def compute_multi(res_root, genome_to_species):
    rows = []
    species_set = set(genome_to_species.values())
    n_species_db = len(species_set)
    for depth in MULTI_DEPTHS:
        truth_species = load_truth_multi(depth)
        sample_key = f"multi/{depth}_sample01"
        pred_species = set()
        for sp_dir in (res_root / sample_key).iterdir():
            if not sp_dir.is_dir():
                continue
            species = sp_dir.name
            for report_path in sp_dir.rglob("final_report.txt"):
                for row in read_final_report(report_path):
                    strain = row.get("Strain_Name", "").strip()
                    if not strain:
                        continue
                    sp = genome_to_species.get(strain)
                    if not sp:
                        base = re.split(r"[._]\d+$", strain)[0]
                        sp = genome_to_species.get(base)
                    if sp:
                        pred_species.add(sp)

        tp = len(truth_species & pred_species)
        fp = len(pred_species - truth_species)
        fn = len(truth_species - pred_species)
        n_truth = len(truth_species)
        n_pred = len(pred_species)
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

        rows.append({
            "depth": depth,
            "sample": "sample01",
            "n_species_db": n_species_db,
            "n_truth": n_truth,
            "n_pred": n_pred,
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "total_wall_s": "",
            "max_rss_mb": "",
        })
    return rows


def main():
    res_root = RES_ROOT if RES_ROOT.exists() else LOCAL_FALLBACK
    if not res_root.exists():
        raise FileNotFoundError(f"StrainScan results not found at {RES_ROOT} or {LOCAL_FALLBACK}")

    genome_to_species = load_genome_to_species()

    single = compute_single(res_root)
    FIGRAW.mkdir(parents=True, exist_ok=True)
    with open(FIGRAW / "strainscan_single_persample.tsv", "w", newline="") as fh:
        fieldnames = ["species", "sample", "strategy", "k", "rep", "depth",
                      "n_truth_cl", "n_pred_cl", "TP", "FP", "FN",
                      "precision", "recall", "f1", "wall_s", "maxrss_mb"]
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in single:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {FIGRAW / 'strainscan_single_persample.tsv'} ({len(single)} rows)")

    multi = compute_multi(res_root, genome_to_species)
    with open(FIGRAW / "strainscan_multi_persample.tsv", "w", newline="") as fh:
        fieldnames = ["depth", "sample", "n_species_db", "n_truth", "n_pred",
                      "TP", "FP", "FN", "precision", "recall", "f1",
                      "total_wall_s", "max_rss_mb"]
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in multi:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {FIGRAW / 'strainscan_multi_persample.tsv'} ({len(multi)} rows)")


if __name__ == "__main__":
    main()
