#!/usr/bin/env python3
"""Cross-check false-positive species across the four mock communities.

Uses the existing profile-mode .pred files (which lack a species column) and
maps cluster IDs to species via the corresponding members.tsv files.
"""
import csv, glob, json, re, sys
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
RAW = PAPER / "work" / "mock_retest" / "Strain2bScan-raw-data"

MEMBERS = {
    "164_all":  RAW / "MSA_combined164_all_flat.members.tsv",
    "120_all":  RAW / "MSA1002_combined_all_flat.members.tsv",
    "164_bcgi": RAW / "MSA_combined164_bcgi_cont.members.tsv",
    "120_bcgi": RAW / "MSA1002_combined_bcgi_cont.members.tsv",
}

TRUTH_DIR = RAW / "Ground_truth"


def load_members(path):
    c2gs = {}
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            g, c = ln.rstrip("\n").split("\t")[:2]
            c2gs.setdefault(c, []).append(g)
    return c2gs


def cluster_to_species(members_path):
    c2gs = load_members(members_path)
    c2sp = {}
    for c, gs in c2gs.items():
        sps = {g.split("__")[0] for g in gs}
        # take the species that appears most often in this cluster
        c2sp[c] = max(sps, key=lambda s: sum(1 for g in gs if g.startswith(s + "__")))
    return c2sp


def atcc_to_species(atcc):
    # Search members files for a genome name containing the ATCC string.
    for members in MEMBERS.values():
        with open(members) as f:
            for ln in f:
                if ln.startswith("#") or not ln.strip():
                    continue
                g = ln.split("\t")[0]
                if atcc in g:
                    return g.split("__")[0]
    return None


def truth_species(mock):
    f = TRUTH_DIR / f"{mock}_ground_truth.txt"
    atccs = []
    with open(f) as fh:
        for ln in fh:
            first = ln.lstrip()
            if not first or first.startswith("base_counts") or first.startswith("genome"):
                continue
            atccs.append(ln.split("\t")[0].strip())
    species = set()
    for atcc in atccs:
        sp = atcc_to_species(atcc)
        if sp:
            species.add(sp)
    return species


def predicted_species_from_profile(pred_path, c2sp):
    species = set()
    with open(pred_path) as f:
        hdr = f.readline()
        if not hdr.startswith("#cluster"):
            f.seek(0)
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            c = ln.split("\t")[0].strip()
            if c in c2sp:
                species.add(c2sp[c])
    return species


def main():
    if not RAW.exists():
        print(f"ERROR: HPC mirror not found at {RAW}", file=sys.stderr)
        sys.exit(1)

    c2sp = {k: cluster_to_species(v) for k, v in MEMBERS.items()}

    # map sample prefix to members key
    key_for = {
        "brad164": "164_bcgi",
        "brad120": "120_bcgi",
        "wms164":  "164_all",
        "shot120": "120_all",
    }

    mock_fps = {}
    for mock in ["MSA1002", "MSA1003", "MSA1005", "MSA1007"]:
        truth = truth_species(mock)
        fps = set()
        for subdir, key in key_for.items():
            d = RAW / "wms_analysis" / subdir
            if not d.exists():
                continue
            for pred in d.glob(f"*{mock}*.pred"):
                pred_sp = predicted_species_from_profile(pred, c2sp[key])
                fps.update(pred_sp - truth)
        mock_fps[mock] = fps
        print(f"== {mock} FP ({len(fps)}) ==")
        for sp in sorted(fps):
            print(f"  {sp}")
        print()

    print("== cross-mock FP overlap ==")
    all_fp = set().union(*mock_fps.values())
    for sp in sorted(all_fp):
        present = [m for m in mock_fps if sp in mock_fps[m]]
        print(f"{sp:<35} in {','.join(present)}")

    print("\n== FP that are true in another mock (index-hopping candidates) ==")
    for m, fps in mock_fps.items():
        truth_others = set().union(*(truth_species(o) for o in mock_fps if o != m))
        hop = fps & truth_others
        print(f"{m}: {len(hop)}/{len(fps)} -> {sorted(hop)}")


if __name__ == "__main__":
    main()
