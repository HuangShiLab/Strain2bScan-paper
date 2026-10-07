#!/usr/bin/env python3
"""Score StrainGST calls in StrainGE's own 0.90-reference cluster space."""
from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/strainge_rerun"
RUNS = ROOT / "work/strainge_rerun/runs"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def cluster_map(path: Path) -> dict[str, str]:
    genome_to_rep = {}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if fields and fields[0]:
            for genome in fields:
                label = Path(genome).name
                if label.endswith(".hdf5"):
                    label = label[:-5]
                genome_to_rep[label] = fields[0]
    return genome_to_rep


def databases(scope: str) -> dict[str, dict[str, str]]:
    return {row["species"]: row for row in read_tsv(OUT / f"{scope}_database_build.tsv")}


def parse_straingst(path: Path) -> tuple[bool, dict[str, float], dict[str, float]]:
    """Return completed status and all/default-threshold predicted abundances."""
    abundances = {}
    thresholded = {}
    if not path.exists() or not path.stat().st_size:
        return False, abundances, thresholded
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 15 or not fields[0].isdigit():
            continue
        try:
            abundance = float(fields[11]) / 100.0
            score = float(fields[14])
        except (ValueError, IndexError):
            continue
        strain = fields[1]
        if score >= 0.02:
            abundances[strain] = max(abundances.get(strain, 0.0), abundance)
            if abundance >= 1e-4:
                thresholded[strain] = max(thresholded.get(strain, 0.0), abundance)
    return True, abundances, thresholded


def collapse(profile: dict[str, float], g2rep: dict[str, str]) -> dict[str, float]:
    out = defaultdict(float)
    for genome, abundance in profile.items():
        if genome in g2rep:
            out[g2rep[genome]] += abundance
    return dict(out)


def metrics(truth: set[str], pred: set[str]) -> tuple[int, int, int, float, float, float]:
    tp, fp, fn = len(truth & pred), len(pred - truth), len(truth - pred)
    precision = tp / (tp + fp) if tp + fp else math.nan
    recall = tp / (tp + fn) if tp + fn else math.nan
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return tp, fp, fn, precision, recall, f1


def profile_metrics(truth: dict[str, float], pred: dict[str, float]) -> tuple[float, float]:
    keys = set(truth) | set(pred)
    if not keys:
        return math.nan, math.nan
    tsum = sum(truth.values()) or 1.0
    psum = sum(pred.values()) or 1.0
    t = {k: truth.get(k, 0.0) / tsum for k in keys}
    p = {k: pred.get(k, 0.0) / psum for k in keys}
    bray = sum(abs(t[k] - p[k]) for k in keys) / sum(t[k] + p[k] for k in keys)
    l2_similarity = 1.0 - math.sqrt(sum((t[k] - p[k]) ** 2 for k in keys))
    return 1.0 - bray, l2_similarity


def depth_of(sample: str) -> str:
    match = re.search(r"_d([0-9.]+)$", sample)
    return match.group(1) if match else ""


def score_sim() -> None:
    dbs = databases("simulation")
    maps = {species: cluster_map(ROOT / row["clusters"]) for species, row in dbs.items()}
    manifest = read_tsv(ROOT / "results/strainscan_rerun/run_manifest.tsv")
    manifest = [r for r in manifest if (r["kind"] == "single" and "_rep1_" in r["sample"]) or r["kind"] == "multi"]
    rows = []
    for item in manifest:
        sample = item["sample"]
        outdir = RUNS / "simulation" / sample
        truth_profiles: dict[str, dict[str, float]] = {}
        if item["kind"] == "single":
            for line in (ROOT / item["truth"]).read_text().splitlines():
                if line.startswith("#") or not line:
                    continue
                fields = line.split("\t")
                truth_profiles.setdefault(item["species"], {})
                truth_profiles[item["species"]][fields[0]] = float(fields[2])
        else:
            for line in (ROOT / item["truth"]).read_text().splitlines():
                if line.startswith("#") or not line:
                    continue
                fields = line.split("\t")
                truth_profiles.setdefault(fields[0], {})
                truth_profiles[fields[0]][fields[1]] = float(fields[3])
        pred_profiles = {}
        completed = True
        species_names = [item["species"]] if item["kind"] == "single" else sorted(dbs)
        for species in species_names:
            done, pred, _ = parse_straingst(outdir / f"straingst_{species}.tsv")
            completed &= done
            pred_profiles[species] = collapse(pred, maps[species])
        truth_sets = {sp: set(collapse(truth_profiles.get(sp, {}), maps[sp])) for sp in species_names}
        pred_sets = {sp: set(pred_profiles.get(sp, {})) for sp in species_names}
        truth_all = set().union(*truth_sets.values())
        pred_all = set().union(*pred_sets.values())
        tp, fp, fn, precision, recall, f1 = metrics(truth_all, pred_all)
        truth_ab = {f"{sp}|{rep}": abundance for sp, profile in truth_profiles.items()
                    for rep, abundance in collapse(profile, maps[sp]).items()}
        pred_ab = {f"{sp}|{rep}": abundance for sp, profile in pred_profiles.items()
                   for rep, abundance in profile.items()}
        bray_similarity, l2_similarity = profile_metrics(truth_ab, pred_ab)
        rows.append({
            "benchmark": "simulation", "kind": item["kind"], "sample": sample,
            "species": item["species"], "depth": depth_of(sample) if item["kind"] == "single" else sample.split("_")[1],
            "status": "completed" if completed else "failed", "n_truth": len(truth_all), "n_pred": len(pred_all),
            "tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1,
            "bray_curtis_similarity": bray_similarity, "l2_similarity": l2_similarity,
        })
    write_tsv(OUT / "simulation_persample.tsv", list(rows[0]), rows)
    summarize(rows, ["species"], OUT / "simulation_by_species.tsv")
    summarize([r for r in rows if r["kind"] == "single"], ["depth"], OUT / "simulation_single_by_depth.tsv")
    summarize([r for r in rows if r["kind"] == "multi"], ["depth"], OUT / "simulation_multi_by_depth.tsv")


def atcc_reference_name(genome_files: list[str], atcc: str) -> str:
    number = re.sub(r"^ATCC[_ ]*", "", atcc).strip()
    matches = [name for name in genome_files if "ATCC_" in name and re.search(rf"{re.escape(number)}(_|$)", name)]
    if not matches:
        raise KeyError(f"cannot map {atcc} in StrainGE panel")
    return matches[0]


def score_mock() -> None:
    dbs = databases("mock")
    maps = {species: cluster_map(ROOT / row["clusters"]) for species, row in dbs.items()}
    all_names = sorted({genome for mapping in maps.values() for genome in mapping})
    rows = []
    for sample, mock, *_ in [
        ("MSA1002_99", "MSA1002"), ("MSA1003_0", "MSA1003"),
        ("MSA1005_0", "MSA1005"), ("MSA1007_0", "MSA1007"),
    ]:
        outdir = RUNS / "mock" / sample
        truth_file = ROOT / f"results/mock_benchmark_f26f234/Ground_truth/{mock}_ground_truth.txt"
        lines = [line.split("\t") for line in truth_file.read_text().splitlines() if line.strip()]
        truth_raw = {}
        if mock == "MSA1003":
            for fields in lines[1:]:
                truth_raw[fields[0]] = float(fields[1])
        else:
            index = lines[0].index("seq_abd")
            for fields in lines[1:]:
                truth_raw[fields[0]] = float(fields[index])
        truth_profiles = {}
        pred_profiles = {}
        completed = True
        for species, metadata in dbs.items():
            done, pred, thresholded = parse_straingst(outdir / f"straingst_{species}.tsv")
            completed &= done
            pred_profiles[species] = collapse(thresholded, maps[species])
            truth_profiles[species] = {}
        for atcc, abundance in truth_raw.items():
            name = atcc_reference_name(all_names, atcc)
            species = next(sp for sp, mapping in maps.items() if name in mapping)
            truth_profiles[species][name] = abundance
        truth_sets = {sp: set(truth_profiles.get(sp, {})) for sp in dbs}
        pred_sets = {sp: set(pred_profiles.get(sp, {})) for sp in dbs}
        truth_all = set().union(*truth_sets.values())
        pred_all = set().union(*pred_sets.values())
        tp, fp, fn, precision, recall, f1 = metrics(truth_all, pred_all)
        truth_ab = {f"{sp}|{rep}": abundance for sp, profile in truth_profiles.items()
                    for rep, abundance in collapse(profile, maps[sp]).items()}
        pred_ab = {f"{sp}|{rep}": abundance for sp, profile in pred_profiles.items()
                   for rep, abundance in profile.items()}
        bray_similarity, l2_similarity = profile_metrics(truth_ab, pred_ab)
        rows.append({
            "benchmark": "ATCC_DNA_mock", "sample": sample, "mock": mock, "status": "completed" if completed else "failed",
            "n_truth": len(truth_all), "n_pred": len(pred_all), "tp": tp, "fp": fp, "fn": fn,
            "precision": precision, "recall": recall, "f1": f1,
            "bray_curtis_similarity": bray_similarity, "l2_similarity": l2_similarity,
        })
    write_tsv(OUT / "mock_persample.tsv", list(rows[0]), rows)


def median(values: list[float]) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    middle = len(finite) // 2
    return finite[middle] if len(finite) % 2 else (finite[middle - 1] + finite[middle]) / 2


def summarize(rows: list[dict[str, object]], group_fields: list[str], path: Path) -> None:
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(str(row[field]) for field in group_fields)].append(row)
    out = []
    metric_fields = ["precision", "recall", "f1", "bray_curtis_similarity", "l2_similarity"]
    for key, values in sorted(groups.items()):
        item = dict(zip(group_fields, key))
        item.update({
            "n": len(values),
            **{field: median([float(value[field]) for value in values]) for field in metric_fields},
        })
        out.append(item)
    write_tsv(path, group_fields + ["n", *metric_fields], out)


def main() -> int:
    score_sim()
    score_mock()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
