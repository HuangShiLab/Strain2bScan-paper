#!/usr/bin/env python3
"""Create a reviewer-facing run manifest for the 22 exploratory public WGS libraries.

SRA spots/bases and instrument models are from NCBI SRA metadata fetched on 2026-10-11.
Marker/call outcomes are copied from the committed `ADVANTAGE_REPORT.md` and
`PRJNA1517970_body_site_summary.tsv` without re-profiling.
"""
from __future__ import annotations

import csv
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
OUT = PAPER / "results/realworld_cohort_benchmark/public_run_manifest.tsv"

# Conditions follow the committed run report. The three P08 runs are ordered by their SRA
# submitter titles: SRR32121885=W1, SRR32121884=W2, SRR32121883=W3.
CONDITIONS = {
    "SRR32076385": ("PRJNA1191225", "Escherichia coli LHCA45 isolate WGS"),
    "SRR32076393": ("PRJNA1191225", "Escherichia coli LHCA56 isolate WGS"),
    "SRR32076407": ("PRJNA1191225", "Escherichia coli LHCA72 isolate WGS"),
    "SRR32076422": ("PRJNA1191225", "Bifidobacterium breve LHCA81 isolate WGS"),
    "SRR32076429": ("PRJNA1191225", "Bifidobacterium bifidum LHCA82 isolate WGS"),
    "SRR32076444": ("PRJNA1191225", "Bifidobacterium longum LHCA43 isolate WGS"),
    "SRR32121885": ("PRJNA1191223", "preterm infant P08 stool W1"),
    "SRR32121884": ("PRJNA1191223", "preterm infant P08 stool W2"),
    "SRR32121883": ("PRJNA1191223", "preterm infant P08 stool W3"),
    "SRR40410528": ("PRJNA1517970", "meconium"),
    "SRR40410567": ("PRJNA1517970", "vagina"),
    "SRR40410657": ("PRJNA1517970", "meconium"),
    "SRR40410734": ("PRJNA1517970", "vagina"),
    "SRR40410753": ("PRJNA1517970", "vagina"),
    "SRR40410776": ("PRJNA1517970", "extraction blank"),
    "SRR40410851": ("PRJNA1517970", "meconium"),
    "SRR6747958": ("PRJNA288562", "subject T23 vagina GD273"),
    "SRR6748028": ("PRJNA288562", "subject T23 distal gut GD84"),
    "SRR6748033": ("PRJNA288562", "subject T23 distal gut GD273"),
    "SRR6748088": ("PRJNA288562", "subject T23 saliva GD273"),
    "SRR6748097": ("PRJNA288562", "subject T23 saliva GD84"),
    "SRR6748207": ("PRJNA288562", "subject T23 vagina GD84"),
}

# Generic MSA outcomes from ADVANTAGE_REPORT.md. Body-site panel outcomes supersede generic
# outcomes only for PRJNA1517970 and are retained in separate columns.
GENERIC = {
    "SRR32076385": (219245, 0), "SRR32076393": (242121, 0), "SRR32076407": (183313, 0),
    "SRR32076422": (142707, 0), "SRR32076429": (129681, 0), "SRR32076444": (65183, 0),
    "SRR32121883": (830986, 6), "SRR32121884": (694111, 2), "SRR32121885": (748209, 1),
    "SRR40410528": (116410, 0), "SRR40410567": (9693, 0), "SRR40410657": (31171, 0),
    "SRR40410734": (24067, 0), "SRR40410753": (14802, 0), "SRR40410776": (464, 0),
    "SRR40410851": (23718, 0), "SRR6747958": (226983, 1), "SRR6748028": (1918932, 7),
    "SRR6748033": (1595769, 2), "SRR6748088": (876200, 11), "SRR6748097": (533573, 13),
    "SRR6748207": (45805, 0),
}

# SRA metadata fetched with NCBI EFetch on 2026-10-11.
SRA = {
    "SRR32076385": (530692, 157695345, 148.58, "Illumina HiSeq 2500"),
    "SRR32076393": (694301, 206230492, 148.52, "Illumina HiSeq 2500"),
    "SRR32076407": (405951, 120852478, 148.85, "Illumina HiSeq 2500"),
    "SRR32076422": (292272, 87006759, 148.85, "Illumina HiSeq 2500"),
    "SRR32076429": (265039, 78810882, 148.68, "Illumina HiSeq 2500"),
    "SRR32076444": (117024, 34917748, 149.19, "Illumina HiSeq 2500"),
    "SRR32121885": (4862900, 1215725000, 125.0, "Illumina HiSeq 2500"),
    "SRR32121884": (4684342, 1171085500, 125.0, "Illumina HiSeq 2500"),
    "SRR32121883": (4403456, 1100864000, 125.0, "Illumina HiSeq 2500"),
    "SRR40410528": (130092, 35779278, 137.52, "Illumina NovaSeq X"),
    "SRR40410567": (29852, 5373843, 90.01, "Illumina NovaSeq X"),
    "SRR40410657": (63704, 12967250, 101.78, "Illumina NovaSeq X"),
    "SRR40410734": (123942, 24263782, 97.88, "Illumina NovaSeq X"),
    "SRR40410753": (71724, 13985625, 97.5, "Illumina NovaSeq X"),
    "SRR40410776": (969, 181397, 93.6, "Illumina NovaSeq X"),
    "SRR40410851": (37274, 8265014, 110.87, "Illumina NovaSeq X"),
    "SRR6747958": (636817, 186911987, 146.75, "Illumina HiSeq 2500"),
    "SRR6748028": (2804411, 824219750, 146.95, "Illumina HiSeq 2500"),
    "SRR6748033": (3803560, 1120072205, 147.24, "Illumina HiSeq 2500"),
    "SRR6748088": (771556, 226981334, 147.09, "Illumina HiSeq 2500"),
    "SRR6748097": (377872, 111169454, 147.1, "Illumina HiSeq 2500"),
    "SRR6748207": (243062, 71385301, 146.85, "Illumina HiSeq 2500"),
}


def main():
    bodysite = {r["run"]: r for r in
                csv.DictReader((PAPER / "results/realworld_cohort_benchmark/PRJNA1517970_body_site_summary.tsv").open(), delimiter="\t")}
    fields = ["run", "project", "condition_or_site", "sra_spots", "sra_bases",
              "mean_read_length_bp", "instrument", "msa_distinct_markers",
              "msa_strain_calls", "bodysite_panel_calls"]
    rows = []
    for run in sorted(CONDITIONS, key=lambda r: (CONDITIONS[r][0], CONDITIONS[r][1], r)):
        project, condition = CONDITIONS[run]
        spots, bases, mean_len, instrument = SRA[run]
        markers, calls = GENERIC[run]
        body_calls = bodysite.get(run, {}).get("body_site_panel_calls", "")
        rows.append({
            "run": run, "project": project, "condition_or_site": condition,
            "sra_spots": spots, "sra_bases": bases, "mean_read_length_bp": mean_len,
            "instrument": instrument, "msa_distinct_markers": markers,
            "msa_strain_calls": calls, "bodysite_panel_calls": body_calls or "--",
        })
    if len(rows) != 22:
        raise SystemExit(f"expected 22 runs, got {len(rows)}")
    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {OUT} ({len(rows)} runs)")


if __name__ == "__main__":
    main()
