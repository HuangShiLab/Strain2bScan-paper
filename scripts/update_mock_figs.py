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
PORT_RAW = PAPER / "work" / "mock_retest" / "Strain2bScan-port-results"
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
    # Guard against the stale per-species cluster naming. Those files silently
    # broke the pred -> members mapping; fail hard rather than guess.
    if c2g and any("__C" in c for c in c2g):
        raise ValueError(
            f"{path} contains per-species cluster names (e.g. 'Species__C0'). "
            "Regenerate the members file so cluster IDs match the DB bare IDs (C0, C1, ...)."
        )
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
    missing = []
    for atcc, ab in t.items():
        num = re.sub(r"^ATCC[_ ]*", "", atcc).strip()
        # anchored match: the ATCC number must be preceded by an underscore (or
        # string start) so ATCC_6940 does not collide with ATCC_16940.
        pat = rf"(?:^|_){re.escape(num)}(_|$)"
        hit = [g for g in all_genomes if "ATCC_" in g and re.search(pat, g)]
        if not hit:
            missing.append(atcc)
            continue
        out[hit[0]] = out.get(hit[0], 0.0) + ab
    if missing:
        raise ValueError(
            f"{mock}: truth ATCC(s) not found in panel genome set: {missing}. "
            "Check that the members file matches the mock community."
        )
    s = sum(out.values())
    return {k: v / s for k, v in out.items()} if s else out

# ------------------------------------------------------------------ parse new Strain2bScan pred (v0.2+ format)
def read_s2b(pred, members_key):
    c2g = MEMBERS[members_key]
    out = {}
    if not os.path.exists(pred):
        return None
    in_tsv = False
    unmatched = set()
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
        if c not in c2g:
            unmatched.add(c)
            continue
        r = rep_of(c2g[c])
        out[r] = out.get(r, 0.0) + ab
    if unmatched:
        raise ValueError(
            f"{pred}: {len(unmatched)} cluster(s) not present in members '{members_key}': "
            f"{sorted(unmatched)[:5]}{'...' if len(unmatched) > 5 else ''}"
        )
    return out

# ------------------------------------------------------------------ metrics
def average_precision(pred, truth):
    """Step-function AP over the precision-recall curve (sklearn convention)."""
    pos = {g for g, v in truth.items() if v > 0}
    if not pos:
        return 0.0
    scores = sorted({v for v in pred.values() if v > 0}, reverse=True)
    prev_recall = 0.0
    ap = 0.0
    for t in scores:
        d = {g for g, v in pred.items() if v >= t and v > 0}
        tp_ = len(d & pos)
        precision = tp_ / len(d) if d else 0.0
        recall = tp_ / len(pos)
        ap += max(0.0, recall - prev_recall) * precision
        prev_recall = recall
    return ap

def pr_at_threshold(pred, truth, genomes, threshold):
    """Precision/recall/F1 when presence requires abundance >= threshold."""
    pos = {g for g in genomes if truth.get(g, 0) > 0}
    det = {g for g, v in pred.items() if v >= threshold}
    tp, fp, fn = len(det & pos), len(det - pos), len(pos - det)
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    return prec, rec, f1

def metrics(pred, truth, genomes):
    pos = {g for g in genomes if truth.get(g, 0) > 0}
    det = {g for g, v in pred.items() if v > 0}
    tp, fp, fn = len(det & pos), len(det - pos), len(pos - det)
    detected_any = sum(pred.values()) > 0
    if detected_any:
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
        aupr = average_precision(pred, truth)
    else:
        # Empty profile: precision/recall/f1/aupr are undefined, not 0.
        prec = rec = f1 = aupr = float("nan")
    p = np.array([pred.get(g, 0.0) for g in genomes], float)
    t = np.array([truth.get(g, 0.0) for g in genomes], float)
    if p.sum() > 0:
        p = p / p.sum()
    l2 = float(np.sqrt(((p - t) ** 2).sum()))
    bc = float(np.abs(p - t).sum() / (p + t).sum()) if (p + t).sum() else 1.0
    pt, rt, f1t = pr_at_threshold(pred, truth, genomes, 1e-4)
    return dict(TP=tp, FP=fp, FN=fn, precision=prec, recall=rec, f1=f1, aupr=aupr,
                precision_at_1e_4=pt, recall_at_1e_4=rt, f1_at_1e_4=f1t,
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
        res[("Strain2bScan-port", "164-tracegap")] = read_s2b(RAW / "wms_analysis_tracegap" / "brad164" / f"{s}.pred", "164_bcgi")
        if mock in ("MSA1002", "MSA1003"):
            res[("Strain2bScan", "120")] = read_s2b(RAW / "wms_analysis" / "brad120" / f"{s}.pred", "120_bcgi")
            res[("Strain2bScan-port", "120-tracegap")] = read_s2b(RAW / "wms_analysis_tracegap" / "brad120" / f"{s}.pred", "120_bcgi")
    else:
        res[("Strain2bScan", "164")] = read_s2b(RAW / "wms_analysis" / "wms164" / f"{s}.pred", "164_all")
        res[("Strain2bScan-port", "164-tracegap")] = read_s2b(RAW / "wms_analysis_tracegap" / "wms164" / f"{s}.pred", "164_all")
        if mock in ("MSA1002", "MSA1003"):
            res[("Strain2bScan", "120")] = read_s2b(RAW / "wms_analysis" / "shot120" / f"{s}.pred", "120_all")
            res[("Strain2bScan-port", "120-tracegap")] = read_s2b(RAW / "wms_analysis_tracegap" / "shot120" / f"{s}.pred", "120_all")
        # strainscan-port layers only: the default path is identical to the
        # Strain2bScan default, so plotting it again would duplicate rows.
        res[("Strain2bScan-port", "164-layers")] = read_s2b(PORT_RAW / "mock" / "wms164_port_layers" / f"{s}.pred", "164_all")
        if mock in ("MSA1002", "MSA1003"):
            res[("Strain2bScan-port", "120-layers")] = read_s2b(PORT_RAW / "mock" / "shot120_port_layers" / f"{s}.pred", "120_all")
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
            if tool == "Strain2bScan":
                base = variant
            elif tool == "Strain2bScan-port":
                base = variant.split("-")[0]
            else:
                base = None

            if base is None:
                genomes = GENOMES164
            elif pr["kind"] == "2bRAD":
                key = "120_bcgi" if base == "120" else "164_bcgi"
                genomes = genome_set(key)
            else:
                key = "120_all" if base == "120" else "164_all"
                genomes = genome_set(key)
            truth = atcc_genomes(mock, genomes)
            m = metrics(prof, truth, genomes)
            m.update(sample=pr["sample"], mock=mock, kind=pr["kind"], tool=tool, variant=variant)
            rows.append(m)
            profiles[f"{pr['sample']}|{tool}|{variant}"] = prof
            profiles[f"{pr['sample']}|truth|{variant}"] = truth

    # merge: drop all old Strain2bScan / Strain2bScan-port rows (they are
    # regenerated from raw pred here), keep old comparison-tool rows and truth.
    s2bs_tools = {"Strain2bScan", "Strain2bScan-port"}
    merged_metrics = [r for r in old_metrics if r["tool"] not in s2bs_tools]
    merged_metrics += rows

    new_profile_keys = set(profiles.keys())
    merged_profiles = {k: v for k, v in old_profiles.items() if k not in new_profile_keys and not k.split("|")[1] in s2bs_tools}
    merged_profiles.update(profiles)

    # Recompute threshold-based metrics for every row from the merged profiles so
    # comparison tools (StrainScan, inStrain) get the new columns too.
    def genomes_for_row(r):
        if r["kind"] == "2bRAD":
            key = "120_bcgi" if r["variant"] == "120" else "164_bcgi"
            return genome_set(key)
        # WMS: comparison tools always use the 164 panel in this dataset.
        return GENOMES164

    for r in merged_metrics:
        pkey = f"{r['sample']}|{r['tool']}|{r['variant']}"
        tkey = f"{r['sample']}|truth|{r['variant']}"
        if pkey in merged_profiles and tkey in merged_profiles:
            genomes = genomes_for_row(r)
            pt, rt, f1t = pr_at_threshold(merged_profiles[pkey], merged_profiles[tkey], genomes, 1e-4)
            r["precision_at_1e_4"] = pt
            r["recall_at_1e_4"] = rt
            r["f1_at_1e_4"] = f1t

    def _fmt(v):
        if isinstance(v, float):
            if math.isnan(v):
                return "nan"
            return round(v, 4)
        return v

    cols = ["kind", "mock", "sample", "tool", "variant", "TP", "FP", "FN",
            "precision", "recall", "f1", "aupr",
            "precision_at_1e_4", "recall_at_1e_4", "f1_at_1e_4",
            "bray_curtis", "l2"]
    with open(DATA / "fig6_fig12_metrics.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader()
        for r in sorted(merged_metrics, key=lambda x: (x["kind"], x["mock"], x["sample"], x["tool"], x["variant"])):
            w.writerow({k: _fmt(r[k]) for k in cols})
    with open(DATA / "fig6_fig12_profiles.json", "w") as fh:
        json.dump(merged_profiles, fh)

    print(f"wrote {DATA / 'fig6_fig12_metrics.tsv'}")
    print(f"wrote {DATA / 'fig6_fig12_profiles.json'}")

    # regenerate figures
    import subprocess
    for cmd in ["python3 scripts/plot_figs_h.py", "python3 scripts/number_figures.py"]:
        print(f"\n$ {cmd}")
        subprocess.run(cmd, shell=True, cwd=PAPER, check=True)

if __name__ == "__main__":
    main()
