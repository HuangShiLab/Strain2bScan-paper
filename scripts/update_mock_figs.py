#!/usr/bin/env python3
"""Merge freshly re-profiled Strain2bScan mock outputs with existing comparison-tool
results, regenerate Fig6/Fig12 metrics + profiles, and rebuild the manuscript figures.

Inputs (populated by rsync from HPC):
  work/mock_retest/Strain2bScan-raw-data/
    MSA_combined164_all_flat.members.tsv
    MSA1002_combined_all_flat.members.tsv
    MSA_combined164_bcgi_cont.members.tsv
    MSA1002_combined_bcgi_cont.members.tsv
    Ground_truth/{MSA1002,MSA1003,MSA1005,MSA1007}_ground_truth.txt
    wms_analysis/{brad164,brad120,wms164,shot120}/{sample}.pred

Outputs:
  data/fig6_fig12_metrics.tsv
  data/fig6_fig12_profiles.json
  figures/fig6_2brad.{png,pdf}
  figures/fig12_wms_toolcompare.{png,pdf}
  figures/numbered/Fig06_fig6_2brad.{png,pdf}
  figures/numbered/Fig12_fig12_wms_toolcompare.{png,pdf}
"""
import csv, glob, json, os, re, sys, math
import numpy as np
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
RAW = PAPER / "work" / "mock_retest" / "Strain2bScan-raw-data"
DATA = PAPER / "data"
FIG = PAPER / "figures"

# ------------------------------------------------------------------ members / clusters
def load_members(path):
    c2g = {}
    for ln in open(path):
        if ln.startswith("#") or not ln.strip():
            continue
        g, c = ln.rstrip("\n").split("\t")[:2]
        c2g.setdefault(c, []).append(g)
    return c2g

def rep_of(members):
    a = [g for g in members if "ATCC_" in g]
    return a[0] if a else sorted(members)[0]

MEMBERS = {
    "164_all":  load_members(RAW / "MSA_combined164_all_flat.members.tsv"),
    "120_all":  load_members(RAW / "MSA1002_combined_all_flat.members.tsv"),
    "164_bcgi": load_members(RAW / "MSA_combined164_bcgi_cont.members.tsv"),
    "120_bcgi": load_members(RAW / "MSA1002_combined_bcgi_cont.members.tsv"),
}
def genome_set(key):
    return sorted({g for gs in MEMBERS[key].values() for g in gs})
GENOMES164 = genome_set("164_all")

# ------------------------------------------------------------------ ground truth
def load_truth(mock):
    f = RAW / "Ground_truth" / f"{mock}_ground_truth.txt"
    rows = [l.rstrip("\n").split("\t") for l in open(f) if l.strip()]
    hdr = rows[0]
    out = {}
    if mock == "MSA1003":
        for r in rows[1:]:
            if len(r) >= 2 and r[0].strip():
                out[r[0].strip()] = float(r[1])
    else:
        i = hdr.index("seq_abd") if "seq_abd" in hdr else 2
        for r in rows[1:]:
            if len(r) > i and r[0].strip():
                out[r[0].strip()] = float(r[i])
    return out

def atcc_genomes(mock, all_genomes):
    t = load_truth(mock)
    out = {}
    for atcc, ab in t.items():
        num = re.sub(r"^ATCC[_ ]*", "", atcc).strip()
        hit = [g for g in all_genomes if "ATCC_" in g and re.search(rf"{re.escape(num)}(_|$)", g)]
        if hit:
            out[hit[0]] = out.get(hit[0], 0.0) + ab
    s = sum(out.values())
    return {k: v / s for k, v in out.items()} if s else out

# ------------------------------------------------------------------ parse new Strain2bScan pred (v0.2+ format)
def read_s2b(pred, members_key):
    c2g = MEMBERS[members_key]
    out = {}
    if not os.path.exists(pred):
        return None
    in_tsv = False
    for ln in open(pred):
        if ln.startswith("#cluster"):
            in_tsv = True
            continue
        if not in_tsv or not ln.strip():
            continue
        p = ln.rstrip("\n").split("\t")
        if len(p) < 2:
            continue
        c, ab = p[0], p[1]
        try:
            ab = float(ab)
        except ValueError:
            continue
        if c in c2g:
            r = rep_of(c2g[c])
            out[r] = out.get(r, 0.0) + ab
    return out

# ------------------------------------------------------------------ metrics
def metrics(pred, truth, genomes):
    pos = {g for g in genomes if truth.get(g, 0) > 0}
    det = {g for g, v in pred.items() if v > 0}
    tp, fp, fn = len(det & pos), len(det - pos), len(pos - det)
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    pts = []
    for t in [0.0] + sorted({v for v in pred.values() if v > 0}):
        d = {g for g, v in pred.items() if v >= t and v > 0}
        if not d:
            continue
        tp_ = len(d & pos)
        pts.append((tp_ / len(pos) if pos else 0.0, tp_ / len(d)))
    pts = sorted(set(pts))
    aupr = 0.0
    if pts:
        rs = [p[0] for p in pts]; ps = [p[1] for p in pts]
        if rs[0] > 0:
            rs = [0.0] + rs; ps = [ps[0]] + ps
        aupr = float(np.trapz(ps, rs))
    p = np.array([pred.get(g, 0.0) for g in genomes], float)
    t = np.array([truth.get(g, 0.0) for g in genomes], float)
    if p.sum() > 0:
        p = p / p.sum()
    l2 = float(np.sqrt(((p - t) ** 2).sum()))
    bc = float(np.abs(p - t).sum() / (p + t).sum()) if (p + t).sum() else 1.0
    return dict(TP=tp, FP=fp, FN=fn, precision=prec, recall=rec, f1=f1, aupr=aupr,
                bray_curtis=bc, l2=l2)

# ------------------------------------------------------------------ sample manifest
def wms(mock, tag):
    return dict(sample=f"WMS_{mock}_{tag}", mock=mock, kind="WMS")
def brad(mock, tag):
    return dict(sample=f"BcgI_{mock}_{tag}", mock=mock, kind="2bRAD")

WMS_SAMPLES = (
    [wms("MSA1002", t) for t in ["0_100ng_1", "0_100ng_2", "90_100ng_1", "95_100ng_1", "99_100ng_1"]] +
    [wms("MSA1003", t) for t in ["0_100ng_1", "0_100ng_2", "0_100ng_3"]] +
    [wms("MSA1005", t) for t in ["0_100ng_1", "0_100ng_2", "0_100ng_3"]] +
    [wms("MSA1007", t) for t in ["0_100ng_1", "0_100ng_2", "0_100ng_3"]]
)
BRAD_SAMPLES = (
    [brad("MSA1002", t) for t in ["0_1ng_1", "0_0.1ng_1", "0_0.01ng_1", "0_0.001ng_1",
                                   "90_100ng_1", "95_100ng_1", "99_100ng_1", "99.9_100ng_1"]] +
    [brad("MSA1003", t) for t in ["0_100ng_1", "0_100ng_2", "0_100ng_3"]] +
    [brad("MSA1005", t) for t in ["0_100ng_1", "0_100ng_2", "0_100ng_3"]] +
    [brad("MSA1007", t) for t in ["0_100ng_1", "0_100ng_2", "0_100ng_3"]]
)

def collect_new_s2b(pr):
    s, mock, kind = pr["sample"], pr["mock"], pr["kind"]
    res = {}
    if kind == "2bRAD":
        res[("Strain2bScan", "164")] = read_s2b(RAW / "wms_analysis" / "brad164" / f"{s}.pred", "164_bcgi")
        if mock in ("MSA1002", "MSA1003"):
            res[("Strain2bScan", "120")] = read_s2b(RAW / "wms_analysis" / "brad120" / f"{s}.pred", "120_bcgi")
    else:
        res[("Strain2bScan", "164")] = read_s2b(RAW / "wms_analysis" / "wms164" / f"{s}.pred", "164_all")
        if mock in ("MSA1002", "MSA1003"):
            res[("Strain2bScan", "120")] = read_s2b(RAW / "wms_analysis" / "shot120" / f"{s}.pred", "120_all")
    return {k: v for k, v in res.items() if v is not None}

# ------------------------------------------------------------------ main
def main():
    if not (RAW / "MSA_combined164_all_flat.members.tsv").exists():
        print(f"ERROR: HPC mirror not found at {RAW}")
        print("Run rsync from HPC first:")
        print("  rsync -av shihuang@hpc2021.hku.hk:/lustre1/g/aos_shihuang/LU/Strain2bScan-raw-data/"
              "\{MSA_combined164_all_flat.members.tsv,MSA1002_combined_all_flat.members.tsv,"
              "MSA_combined164_bcgi_cont.members.tsv,MSA1002_combined_bcgi_cont.members.tsv,"
              "Ground_truth,wms_analysis\} work/mock_retest/Strain2bScan-raw-data/")
        sys.exit(1)

    old_profiles = json.load(open(DATA / "fig6_fig12_profiles.json"))
    old_metrics = list(csv.DictReader(open(DATA / "fig6_fig12_metrics.tsv"), delimiter="\t"))

    profiles = {}
    rows = []
    ALL = WMS_SAMPLES + BRAD_SAMPLES
    for pr in ALL:
        mock = pr["mock"]
        for (tool, variant), prof in collect_new_s2b(pr).items():
            if tool != "Strain2bScan":
                genomes = GENOMES164
            elif pr["kind"] == "2bRAD":
                key = "120_bcgi" if variant == "120" else "164_bcgi"
                genomes = genome_set(key)
            else:
                key = "120_all" if variant == "120" else "164_all"
                genomes = genome_set(key)
            truth = atcc_genomes(mock, genomes)
            m = metrics(prof, truth, genomes)
            m.update(sample=pr["sample"], mock=mock, kind=pr["kind"], tool=tool, variant=variant)
            rows.append(m)
            profiles[f"{pr['sample']}|{tool}|{variant}"] = prof
            profiles[f"{pr['sample']}|truth|{variant}"] = truth

    # merge: keep old non-Strain2bScan entries and truth, add new Strain2bScan
    new_metric_keys = {(r["sample"], r["tool"], r["variant"]) for r in rows}
    merged_metrics = [r for r in old_metrics if (r["sample"], r["tool"], r["variant"]) not in new_metric_keys]
    merged_metrics += rows

    new_profile_keys = set(profiles.keys())
    merged_profiles = {k: v for k, v in old_profiles.items() if k not in new_profile_keys}
    merged_profiles.update(profiles)

    cols = ["kind", "mock", "sample", "tool", "variant", "TP", "FP", "FN",
            "precision", "recall", "f1", "aupr", "bray_curtis", "l2"]
    with open(DATA / "fig6_fig12_metrics.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader()
        for r in sorted(merged_metrics, key=lambda x: (x["kind"], x["mock"], x["sample"], x["tool"], x["variant"])):
            w.writerow({k: (round(r[k], 4) if isinstance(r[k], float) else r[k]) for k in cols})
    json.dump(merged_profiles, open(DATA / "fig6_fig12_profiles.json", "w"))

    print(f"wrote {DATA / 'fig6_fig12_metrics.tsv'}")
    print(f"wrote {DATA / 'fig6_fig12_profiles.json'}")

    # regenerate figures
    import subprocess
    for cmd in ["python3 scripts/plot_figs_h.py", "python3 scripts/number_figures.py"]:
        print(f"\n$ {cmd}")
        subprocess.run(cmd, shell=True, cwd=PAPER, check=True)

if __name__ == "__main__":
    main()
