#!/usr/bin/env python3
"""Parse StrainScan final_report.txt files produced by run_strainscan_sim_hpc.py
and write strain/species abundance tables for FigS5.

Inputs:
  /lustre1/.../sim_benchmark/strainscan_results/{single,multi}/<sample>/<species>/final_report.txt
  results/genome_to_species.tsv

Outputs:
  figure_raw_data/sim_headtohead/strainscan_single_abundance.tsv
  figure_raw_data/sim_headtohead/strainscan_multi_species_abundance.tsv

Abundance estimation:
  - For each final_report we take the 'Predicted_Depth (Ab*cls_depth)' column as
    an absolute depth estimate.
  - Single-species samples: depths are normalised to sum to 1 across strains.
  - Multi-species samples: strain depths are summed per species, then normalised
    to sum to 1 across species (matching the species-level display used for
    Strain2bScan and inStrain in FigS5).
"""
import csv
import re
from pathlib import Path
from collections import defaultdict

PAPER = Path(__file__).resolve().parent.parent
RES_ROOT = Path("/lustre1/g/aos_shihuang/LU/Strain2bScan-raw-data/sim_benchmark/strainscan_results")
FIGRAW = PAPER / "figure_raw_data" / "sim_headtohead"

# Fallback to a locally-synced copy if the HPC path is unavailable
LOCAL_FALLBACK = PAPER / "work" / "mock_retest" / "Strain2bScan-raw-data" / "sim_benchmark" / "results" / "strainscan"

SINGLE_SAMPLES = [
    ("Escherichia_coli", "Escherichia_coli__diff_k2_rep1_d5"),
    ("Cutibacterium_acnes", "Cutibacterium_acnes__diff_k2_rep1_d5"),
    ("Staphylococcus_epidermidis", "Staphylococcus_epidermidis__diff_k2_rep1_d5"),
    ("Prevotella_copri", "Prevotella_copri__diff_k2_rep1_d5"),
]
MULTI_DEPTHS = ["depth_low", "depth_med", "depth_high"]

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
    path = PAPER / "results" / "genome_to_species.tsv"
    out = {}
    if path.exists():
        for row in csv.DictReader(open(path), delimiter="\t"):
            out[row["genome"]] = normalize_species(row["species"])
    return out


def read_final_report(path):
    """Return list of dicts from a StrainScan final_report.txt.

    StrainScan writes two header variants:
      - single-cluster: Strain_ID, Strain_Name, Cluster_ID,
        Relative_Abundance_Inside_Cluster, Predicted_Depth, ...
      - multi-cluster: ID, Strain_Name, Cluster_ID, Relative_Abundance,
        Predicted_Depth (Enet), Predicted_Depth (Ab*cls_depth), ...
    """
    rows = []
    if not path.exists():
        return rows
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # Normalise depth column name.
            depth = row.get("Predicted_Depth (Ab*cls_depth)", "")
            if depth == "":
                depth = row.get("Predicted_Depth", "")
            if depth == "":
                depth = row.get("Predicted_Depth (Enet)", "")
            row["_depth"] = depth
            rows.append(row)
    return rows


def parse_depth(val):
    try:
        return float(val.strip()) if val else 0.0
    except (ValueError, AttributeError):
        return 0.0


def aggregate_single(res_root):
    out = []
    for species, sample in SINGLE_SAMPLES:
        sample_dir = res_root / "single" / sample / species
        depths = defaultdict(float)
        for report_path in sample_dir.rglob("final_report.txt"):
            for row in read_final_report(report_path):
                strain = row.get("Strain_Name", "").strip()
                depth = parse_depth(row.get("_depth", "0"))
                if strain and depth > 0:
                    depths[strain] += depth
        total = sum(depths.values())
        for strain, depth in sorted(depths.items(), key=lambda x: -x[1]):
            out.append({
                "sample": sample,
                "species": species,
                "strain": strain,
                "predicted_depth": depth,
                "relative_abundance": depth / total if total else 0.0,
            })
    return out


def aggregate_multi(res_root, genome_to_species):
    out = []
    for depth in MULTI_DEPTHS:
        sample = f"{depth}_sample01"
        sample_dir = res_root / "multi" / sample
        species_depths = defaultdict(float)
        for species_dir in sample_dir.iterdir():
            if not species_dir.is_dir():
                continue
            for report_path in species_dir.rglob("final_report.txt"):
                for row in read_final_report(report_path):
                    strain = row.get("Strain_Name", "").strip()
                    depth_val = parse_depth(row.get("_depth", "0"))
                    if not strain or depth_val <= 0:
                        continue
                    sp = genome_to_species.get(strain)
                    if not sp:
                        # strain name may contain version suffix; try stripping
                        base = re.split(r"[._]\d+$", strain)[0]
                        sp = genome_to_species.get(base)
                    if sp:
                        species_depths[sp] += depth_val
        total = sum(species_depths.values())
        for sp, depth_val in sorted(species_depths.items(), key=lambda x: -x[1]):
            out.append({
                "sample": sample,
                "species": sp,
                "predicted_depth": depth_val,
                "relative_abundance": depth_val / total if total else 0.0,
            })
    return out


def main():
    res_root = RES_ROOT if RES_ROOT.exists() else LOCAL_FALLBACK
    if not res_root.exists():
        raise FileNotFoundError(f"StrainScan results not found at {RES_ROOT} or {LOCAL_FALLBACK}")

    genome_to_species = load_genome_to_species()

    single = aggregate_single(res_root)
    FIGRAW.mkdir(parents=True, exist_ok=True)
    with open(FIGRAW / "strainscan_single_abundance.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "species", "strain", "predicted_depth", "relative_abundance"], delimiter="\t")
        w.writeheader()
        for r in single:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {FIGRAW / 'strainscan_single_abundance.tsv'} ({len(single)} rows)")

    multi = aggregate_multi(res_root, genome_to_species)
    with open(FIGRAW / "strainscan_multi_species_abundance.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "species", "predicted_depth", "relative_abundance"], delimiter="\t")
        w.writeheader()
        for r in multi:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {FIGRAW / 'strainscan_multi_species_abundance.tsv'} ({len(multi)} rows)")


if __name__ == "__main__":
    main()
