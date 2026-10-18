#!/usr/bin/env python3
"""Fig 3 replacement: Strain2bScan accuracy over all 2,025 single-species simulations."""
import csv, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

rows=list(csv.DictReader(open("figure_raw_data/sim_headtohead/strain2bscan_single_native_persample.tsv"),delimiter="\t"))
def aggregate(key, order=None):
    g=collections.defaultdict(list)
    for r in rows:g[r[key]].append(r)
    keys=order or sorted(g)
    out=[]
    for k in keys:
        rs=g[k]; vals=[tuple(float(r[m]) for m in ("precision","recall","f1")) for r in rs]
        out.append({"group":k,"n":len(rs),
                    "mean":[sum(v[i] for v in vals)/len(vals) for i in range(3)],
                    "median":[sorted(v[i] for v in vals)[len(vals)//2] for i in range(3)],
                    "perfect":sum(v[2]==1 for v in vals)/len(vals)})
    return out
by_sp=aggregate("species")
DEPTHS=["0.5","1","3","5","10"]; by_dep=aggregate("depth",DEPTHS); by_str=aggregate("strategy",["diff","same"])
fig,axes=plt.subplots(1,3,figsize=(15.5,4.8))
x=np.arange(len(by_sp)); width=.25
for j,(m,label) in enumerate(zip(("precision","recall","f1"),("Precision","Recall","F1"))):
    axes[0].bar(x+(j-1)*width,[d["mean"][j] for d in by_sp],width,label=label)
axes[0].set_xticks(x); axes[0].set_xticklabels([d["group"].replace("_"," ").split(" ")[0][0]+". "+d["group"].replace("_"," ").split(" ")[1] for d in by_sp],rotation=45,ha="right",fontsize=8,style="italic")
axes[0].set_ylim(0,1.08); axes[0].set_ylabel("mean score"); axes[0].set_title(f"Per-species accuracy (n={len(rows)})")
axes[0].legend(ncol=3,fontsize=8); axes[0].grid(axis="y",alpha=.25)
for j,(m,label) in enumerate(zip(("precision","recall","f1"),("Precision","Recall","F1"))):
    axes[1].plot(range(5),[d["mean"][j] for d in by_dep],"o-",label=f"mean {label.lower()}")
    axes[1].plot(range(5),[d["median"][j] for d in by_dep],"s--",alpha=.6,label=f"median {label.lower()}")
axes[1].set_xticks(range(5)); axes[1].set_xticklabels([d["group"]+"×" for d in by_dep]); axes[1].set_xlabel("per-strain depth")
axes[1].set_ylabel("score"); axes[1].set_ylim(.4,1.06); axes[1].set_title("Accuracy across depth"); axes[1].legend(fontsize=7,ncol=2); axes[1].grid(alpha=.25)
strategy_names={"diff":"different 0.95 clusters","same":"same 0.95 cluster"}
for i,d in enumerate(by_str):
    bars=axes[2].bar([i-.18,i,i+.18],[d["mean"][0],d["mean"][1],d["mean"][2]],.24,color=["#1f77b4","#2ca02c","#ff7f0e"])
    axes[2].text(i,max(b.get_height() for b in bars)+.03,f"n={d['n']}",ha="center",fontsize=9)
axes[2].set_xticks(range(2)); axes[2].set_xticklabels([strategy_names[d["group"]] for d in by_str]); axes[2].set_ylim(0,1.08)
axes[2].set_ylabel("mean score"); axes[2].set_title("Accuracy by strain-pair design"); axes[2].grid(axis="y",alpha=.25)
handles=[plt.Rectangle((0,0),1,1,color=c) for c in ["#1f77b4","#2ca02c","#ff7f0e"]]
axes[2].legend(handles,["Precision","Recall","F1"],fontsize=8)
fig.tight_layout()
fig.savefig("figures/source_panels/sim2025_summary.png",dpi=180,bbox_inches="tight")
fig.savefig("figures/source_panels/sim2025_summary.pdf",bbox_inches="tight")
print("wrote figures/source_panels/sim2025_summary.png/.pdf")
