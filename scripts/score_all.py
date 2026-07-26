#!/usr/bin/env python3
"""Comprehensive scorer -> long-form metrics + per-sample profiles for Fig 6 (2bRAD) & Fig 12 (WMS).

Tools & references
------------------
  Strain2bScan  combined tree, two variants reported side by side (user's choice):
                  tree=164 (28-species unified) and tree=120 (20-species, MSA1002/1003 only)
  StrainScan    per-species Cluster_ID (28-species DB), shotgun only
  inStrain      per-genome; native ref120 (MSA1002/1003), native ref164 (MSA1005/1007),
                and a dereplicated ref98 fair-play control (MSA1002/1003)

Everything is scored in a common per-mock label space: positives = the ATCC target genomes
truly present in that mock (truth seq_abd > 0); a cluster/genome detection is a true positive
iff it maps to one of those ATCC targets.

Metrics (per lab precedent):
  precision, recall            2bRAD-M (Genome Biol 2021) & Ye et al. (Cell 2019)
  AUPR                         Ye et al.: integrate PR curve while sweeping abundance threshold 0->1
  L2, Bray-Curtis              2bRAD-M: distance of predicted vs ground-truth (seq_abd) profile
"""
import csv, glob, os, re, sys, json
import numpy as np

RAW = "/Users/macstudio/Downloads/Strain2bScan-raw-data"
SP = "/Users/macstudio/Downloads/Strain2bScan-raw-data/wms_analysis"

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

# WMS (shotgun, in-silico all-enzyme digest) is scored against the all-enzyme trees;
# native BcgI 2bRAD is scored against the BcgI-containment trees. Cluster IDs are NOT
# comparable between the two, so each modality must use its own members file.
MEMBERS = {
    "164_all":  load_members(f"{RAW}/MSA_combined164_all.members.tsv"),
    "120_all":  load_members(f"{RAW}/MSA1002_combined_all.members.tsv"),
    "164_bcgi": load_members(f"{RAW}/MSA_combined164_bcgi_cont.members.tsv"),
    "120_bcgi": load_members(f"{RAW}/MSA1002_combined_bcgi_cont.members.tsv"),
}
def genome_set(key):
    return sorted({g for gs in MEMBERS[key].values() for g in gs})
GENOMES164 = genome_set("164_all")

# ------------------------------------------------------------------ ground truth
def load_truth(mock):
    f = f"{RAW}/Ground_truth/{mock}_ground_truth.txt"
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
    """Map each truth ATCC id -> the DB genome name (ATCC target). Returns {genome: seq_abd}."""
    t = load_truth(mock)
    out = {}
    for atcc, ab in t.items():
        num = re.sub(r"^ATCC[_ ]*", "", atcc).strip()
        hit = [g for g in all_genomes if "ATCC_" in g and re.search(rf"{re.escape(num)}(_|$)", g)]
        if hit:
            out[hit[0]] = out.get(hit[0], 0.0) + ab
    s = sum(out.values())
    return {k: v / s for k, v in out.items()} if s else out

# ------------------------------------------------------------------ tool readers
def read_s2b(pred, members_key):
    c2g = MEMBERS[members_key]
    out = {}
    if not os.path.exists(pred):
        return None
    for ln in open(pred):
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.rstrip("\n").split("\t")
        c, ab = p[0], float(p[1])
        if c in c2g:
            r = rep_of(c2g[c])
            out[r] = out.get(r, 0.0) + ab
    return out

def read_ss_maps(path):
    maps, cur = {}, None
    for ln in open(path):
        if ln.startswith("@@"):
            cur = ln[2:].strip(); maps[cur] = {}
        elif ln.strip() and cur:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 3:
                maps[cur]["C" + p[0]] = p[2].split(",")
    return maps
SS_MAPS = read_ss_maps(f"{SP}/ss_maps/all_maps.txt")

def read_strainscan(sample, all_genomes):
    got = False
    out = {}
    for sp, cl2g in SS_MAPS.items():
        frs = glob.glob(f"{SP}/ss_out/{sample}__{sp}/final_report.txt")
        if frs:
            got = True
        for fr in frs:
            for ln in open(fr):
                if ln.startswith("Strain_ID") or not ln.strip():
                    continue
                p = ln.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                cid, depth = p[2].strip(), float(p[4])
                full = []
                for m in cl2g.get(cid, []):
                    full += [g for g in all_genomes if g.split("__", 1)[-1] == m]
                if not full:
                    continue
                r = rep_of(full)
                out[r] = out.get(r, 0.0) + depth
    if not got:
        return None
    s = sum(out.values())
    return {k: v / s for k, v in out.items()} if s else {}

def read_instrain(isdir, sample, breadth_min=0.5, popani_min=0.99999):
    f = f"{isdir}/{sample}.IS/output/{sample}.IS_genome_info.tsv"
    if not os.path.exists(f):
        return None
    out = {}
    for r in csv.DictReader(open(f), delimiter="\t"):
        try:
            br, cov, pa = float(r["breadth"]), float(r["coverage"]), float(r["popANI_reference"] or 0)
        except Exception:
            continue
        if br >= breadth_min and pa >= popani_min:
            out[r["genome"]] = cov
    s = sum(out.values())
    return {k: v / s for k, v in out.items()} if s else {}

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

def collect(pred):
    """Return {(sample,tool,tree): profile} for one sample dict; None profiles skipped."""
    s, mock, kind = pred["sample"], pred["mock"], pred["kind"]
    res = {}
    if kind == "2bRAD":
        # Strain2bScan only (StrainScan/inStrain cannot process reduced-representation 2bRAD)
        res[("Strain2bScan", "164")] = read_s2b(f"{SP}/brad164/{s}.pred", "164_bcgi")
        if mock in ("MSA1002", "MSA1003"):
            res[("Strain2bScan", "120")] = read_s2b(f"{SP}/brad120/{s}.pred", "120_bcgi")
    else:
        res[("Strain2bScan", "164")] = read_s2b(f"{SP}/wms164/{s}.pred", "164_all")
        if mock in ("MSA1002", "MSA1003"):
            res[("Strain2bScan", "120")] = read_s2b(f"{SP}/shot120/{s}.pred", "120_all")
        res[("StrainScan", "-")] = read_strainscan(s, GENOMES164)
        # inStrain is run on the dereplicated reference appropriate for each mock
        if mock in ("MSA1002", "MSA1003"):
            res[("Strain2bScan", "120")] = read_s2b(f"{SP}/shot120/{s}.pred", "120_all")
            res[("inStrain", "derep98")] = read_instrain(f"{SP}/instrain", s)
        else:
            res[("inStrain", "derep164")] = read_instrain(f"{SP}/instrain", s)
    return {k: v for k, v in res.items() if v is not None}

# ------------------------------------------------------------------ main
def main():
    rows, profiles = [], {}
    ALL = WMS_SAMPLES + BRAD_SAMPLES
    for pr in ALL:
        mock = pr["mock"]
        for (tool, variant), prof in collect(pr).items():
            # genome set is identical across all/bcgi trees of the same size (only clustering differs)
            genomes = genome_set("120_all") if (tool == "Strain2bScan" and variant == "120") else GENOMES164
            truth = atcc_genomes(mock, genomes)
            m = metrics(prof, truth, genomes)
            m.update(sample=pr["sample"], mock=mock, kind=pr["kind"], tool=tool, variant=variant)
            rows.append(m)
            profiles[f"{pr['sample']}|{tool}|{variant}"] = prof
            profiles[f"{pr['sample']}|truth|{variant}"] = truth
    cols = ["kind", "mock", "sample", "tool", "variant", "TP", "FP", "FN",
            "precision", "recall", "f1", "aupr", "bray_curtis", "l2"]
    out = f"{SP}/all_metrics.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader()
        for r in rows:
            w.writerow({k: (round(r[k], 4) if isinstance(r[k], float) else r[k]) for k in cols})
    json.dump(profiles, open(f"{SP}/all_profiles.json", "w"))
    # human summary
    for kind in ("WMS", "2bRAD"):
        print(f"\n===== {kind} =====")
        print(f"{'sample':26} {'tool':12} {'var':8} {'P':>5} {'R':>5} {'F1':>5} {'AUPR':>5} {'BC':>5} {'L2':>5}")
        for r in rows:
            if r["kind"] != kind:
                continue
            print(f"{r['sample']:26} {r['tool']:12} {r['variant']:8} "
                  f"{r['precision']:5.2f} {r['recall']:5.2f} {r['f1']:5.2f} "
                  f"{r['aupr']:5.2f} {r['bray_curtis']:5.2f} {r['l2']:5.2f}")
    print(f"\n-> {out}")

if __name__ == "__main__":
    main()
