#!/usr/bin/env python3
"""Fig 5 (multi-species): enzyme count vs strain-ID accuracy across 14 sim-pool species.
Reads results/enzyme_sweep_multi.tsv; writes figures/enzyme_sweep.{png,pdf}."""
import csv, statistics as st, os
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, numpy as np
S2B="#1f77b4"
rows=list(csv.DictReader(open("results/enzyme_sweep_multi.tsv"),delimiter="\t"))
ks=sorted({int(r["n_enzymes"]) for r in rows})
bysp=defaultdict(dict)
for r in rows: bysp[r["species"]][int(r["n_enzymes"])]=r
def med(k,f): return st.median([float(bysp[s][k][f]) for s in bysp if k in bysp[s]])
fig,ax=plt.subplots(1,2,figsize=(11,4.4))
# A precision + recall vs enzymes (median + per-species faint recall)
for s in bysp:
    xs=[k for k in ks if k in bysp[s]]; ax[0].plot(xs,[float(bysp[s][k]["recall"]) for k in xs],color=S2B,alpha=0.12,lw=1)
ax[0].plot(ks,[med(k,"recall") for k in ks],"o-",color=S2B,lw=2.5,ms=8,label="recall (median)")
ax[0].plot(ks,[med(k,"precision") for k in ks],"s--",color="#2ca02c",lw=2,ms=7,label="precision (median)")
ax[0].set_xscale("log",base=2); ax[0].set_xticks(ks); ax[0].set_xticklabels(ks)
ax[0].set_xlabel("number of 2bRAD enzymes"); ax[0].set_ylabel("accuracy"); ax[0].set_ylim(0,1.05)
ax[0].set_title("A  Precision/recall vs enzyme set (14 species)"); ax[0].legend(loc="lower right"); ax[0].grid(alpha=.3)
for k in ks: ax[0].annotate(f"{med(k,'recall'):.2f}",(k,med(k,'recall')),textcoords="offset points",xytext=(0,7),ha="center",fontsize=8)
# B strain-specific markers vs enzymes
ax[1].plot(ks,[med(k,"strain_specific") for k in ks],"o-",color=S2B,lw=2.5,ms=8)
for s in bysp:
    xs=[k for k in ks if k in bysp[s]]; ax[1].plot(xs,[float(bysp[s][k]["strain_specific"]) for k in xs],color=S2B,alpha=0.12,lw=1)
ax[1].set_xscale("log",base=2); ax[1].set_yscale("log"); ax[1].set_xticks(ks); ax[1].set_xticklabels(ks)
ax[1].set_xlabel("number of 2bRAD enzymes"); ax[1].set_ylabel("strain-specific markers (median, log)")
ax[1].set_title("B  Marker yield vs enzyme set"); ax[1].grid(alpha=.3)
fig.suptitle("2bRAD enzyme set is a resolution/cost knob — precision 1.0, recall rises with enzymes (14 species)",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.95])
os.makedirs("figures",exist_ok=True)
fig.savefig("figures/enzyme_sweep.png",dpi=150,bbox_inches="tight"); fig.savefig("figures/enzyme_sweep.pdf",bbox_inches="tight")
print("wrote figures/enzyme_sweep.{png,pdf}")
