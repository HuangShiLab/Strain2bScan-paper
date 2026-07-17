#!/usr/bin/env python3
"""Fig 12 panels C/D: strain-ID accuracy and wall time of Strain2bScan vs StrainScan vs inStrain on the
MSA-1002-even shotgun sample. Reads results/mock_tool_comparison.tsv. (The full 4-panel Fig 12 is
assembled together with panels A/B in scripts/plot_mock_strainid.py; this documents the C/D data.)"""
import csv
rows=[r for r in csv.DictReader(open("results/mock_tool_comparison.tsv"),delimiter="\t")
      if r["sample"]=="WMS_MSA1002_0_100ng_2"]
for r in rows:
    print(r["tool"], "acc", r["strains_resolved"], "time_s", r["wall_time_s"], "env", r["env"])
