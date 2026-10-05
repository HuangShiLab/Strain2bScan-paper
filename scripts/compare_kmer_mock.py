#!/usr/bin/env python3
"""Compare kmer-mode vs enzyme-mode Strain2bScan results on the mock communities.

Reads:
  work/mock_retest/Strain2bScan-raw-data/wms_analysis/wms164_port_kmer/*.pred
  work/mock_retest/Strain2bScan-raw-data/wms_analysis/brad164_port_kmer/*.pred
  work/mock_retest/Strain2bScan-raw-data/kmer_members/{k31,k15}/*.members.tsv
Enzyme-mode reference numbers are recomputed with the same code path
(update_mock_figs.read_s2b / metrics) so the comparison is apples-to-apples.

kmer preds use Species__Cx cluster names (per-species multi-profile DBs), so
members are loaded here as {species}__{cluster} -> [genomes] instead of going
through update_mock_figs.load_members (which rejects __C names by design).
"""
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PAPER / "scripts"))
import update_mock_figs as umf

RAW = PAPER / "work" / "mock_retest" / "Strain2bScan-raw-data"


def load_kmer_members(kind_dir):
    c2g = {}
    for f in sorted((RAW / "kmer_members" / kind_dir).glob("*.members.tsv")):
        species = f.name[: -len(".members.tsv")]
        for ln in open(f):
            if ln.startswith("#") or not ln.strip():
                continue
            g, c = ln.rstrip("\n").split("\t")[:2]
            c2g.setdefault(f"{species}__{c}", []).append(g)
    return c2g


def read_kmer_pred(pred, c2g):
    out = {}
    if not pred.exists():
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
        try:
            ab = float(p[1])
        except ValueError:
            continue
        # Ambiguous call "Species__C0|C1": the species is resolved but the
        # cluster is not. Score it against the union of the candidate clusters
        # (rep_of prefers the ATCC genome), same convention as a single call.
        name = p[0]
        if "|" in name:
            sp, _, rest = name.partition("__")
            parts = rest.split("|")
            gs = [g for part in parts for g in c2g.get(f"{sp}__{part}", [])]
            if not gs:
                unmatched.add(name)
                continue
            r = umf.rep_of(gs)
            out[r] = out.get(r, 0.0) + ab
            continue
        if name not in c2g:
            unmatched.add(name)
            continue
        r = umf.rep_of(c2g[name])
        out[r] = out.get(r, 0.0) + ab
    if unmatched:
        raise ValueError(f"{pred}: clusters not in members: {sorted(unmatched)[:5]}")
    return out


def fmt(m):
    if m is None:
        return "      -       -       -       -       -       -       -"
    return (f"{m['TP']:3d}/{m['FP']:<3d}/{m['FN']:<3d}  "
            f"{m['precision']:.3f}  {m['recall']:.3f}  {m['f1']:.3f}  "
            f"{m['aupr']:.3f}  {m['bray_curtis']:.3f}")


def main():
    kmer_wms = load_kmer_members("k31")
    kmer_brad = load_kmer_members("k15")
    genomes_kmer = sorted({g for gs in kmer_wms.values() for g in gs})
    # sanity: kmer panel must cover the same genomes as the enzyme panel
    assert genomes_kmer == umf.genome_set("164_all"), "kmer panel != enzyme 164 panel"

    # The kmer preds passed through flatten_multiprofile.py with a global-abundance
    # floor (0.001; 0.0 for the MSA1003 titration). Apply the SAME floor to the
    # enzyme-mode profiles so the comparison is threshold-matched.
    def floor_of(sample):
        return 0.0 if "MSA1003" in sample else 0.001

    def floored(prof, sample):
        if prof is None:
            return None
        t = floor_of(sample)
        return {g: a for g, a in prof.items() if a >= t}

    hdr = f"{'sample':28s} {'mode':22s} {'TP/FP/FN':12s} {'prec':6s} {'rec':6s} {'f1':6s} {'aupr':6s} {'BC':6s}"
    print(hdr)
    print("-" * len(hdr))

    for pr in umf.WMS_SAMPLES + umf.BRAD_SAMPLES:
        s, mock, kind = pr["sample"], pr["mock"], pr["kind"]
        truth = umf.atcc_genomes(mock, genomes_kmer)
        rows = []
        if kind == "WMS":
            km = read_kmer_pred(RAW / "wms_analysis" / "wms164_port_kmer" / f"{s}.pred", kmer_wms)
            enz = umf.read_s2b(RAW / "wms_analysis" / "wms164" / f"{s}.pred", "164_all")
            lay = umf.read_s2b(PAPER / "work/mock_retest/Strain2bScan-port-results/mock/wms164_port_layers" / f"{s}.pred", "164_all")
            rows = [("kmer k31s100", km), ("enzyme (default)", enz), ("enzyme (layers)", lay)]
        else:
            km = read_kmer_pred(RAW / "wms_analysis" / "brad164_port_kmer" / f"{s}.pred", kmer_brad)
            enz = umf.read_s2b(RAW / "wms_analysis" / "brad164" / f"{s}.pred", "164_bcgi")
            rows = [("kmer k15s100", km), ("enzyme (default)", enz)]
        for i, (label, prof) in enumerate(rows):
            if kind == "WMS":
                prof = floored(prof, s)
            m = umf.metrics(prof, truth, genomes_kmer) if prof is not None else None
            print(f"{s if i == 0 else '':28s} {label:22s} {fmt(m)}")
        print()


if __name__ == "__main__":
    main()
