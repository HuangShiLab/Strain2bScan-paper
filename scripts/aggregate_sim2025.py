#!/usr/bin/env python3
"""Aggregate the complete 2,025-sample single-species simulation benchmark."""
import csv, collections, statistics
from pathlib import Path
SRC = Path("figure_raw_data/sim_headtohead/strain2bscan_single_native_persample.tsv")
OUT = Path("results/sim2025_summary")
OUT.mkdir(parents=True, exist_ok=True)
rows = list(csv.DictReader(SRC.open(), delimiter="\t"))
fields = ["group","n","mean_precision","mean_recall","mean_f1","median_precision","median_recall","median_f1","fraction_perfect_f1"]
def write(path, groups):
    with path.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t"); w.writeheader()
        for key,rs in groups:
            vals=[tuple(float(r[m]) for m in ("precision","recall","f1")) for r in rs]
            means=[sum(v[i] for v in vals)/len(vals) for i in range(3)]
            meds=[statistics.median(v[i] for v in vals) for i in range(3)]
            w.writerow({"group":key,"n":len(vals),
                        **{f"mean_{m}":means[i] for i,m in enumerate(("precision","recall","f1"))},
                        **{f"median_{m}":meds[i] for i,m in enumerate(("precision","recall","f1"))},
                        "fraction_perfect_f1":sum(v[2]==1 for v in vals)/len(vals)})
for name,key in [("by_species.tsv","species"),("by_depth.tsv","depth"),("by_strategy.tsv","strategy")]:
    groups=collections.defaultdict(list)
    for r in rows: groups[r[key]].append(r)
    write(OUT/name,sorted(groups.items()))
print(f"wrote {OUT} for {len(rows)} samples")
