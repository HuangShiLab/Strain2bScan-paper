#!/usr/bin/env python3
"""Figure 2: genome-distance preservation summary and representative rank-rank panels."""
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAPER=os.path.expanduser("~/Downloads/Strain2bScan-paper")
PD=Path(PAPER)/"results/pairdist"
TSV=f"{PAPER}/results/mash_2brad_vs_16s_recomputed.tsv"
OUT=f"{PAPER}/figures/source_panels/mash_2brad_vs_16s"

rows=[l.rstrip("\n").split("\t") for l in open(TSV)
      if l.strip() and not l.startswith("#") and not l.startswith("species")]
stat={r[0]:{"n":int(r[1]),"rho2":float(r[3]),"lo2":float(r[4]),"hi2":float(r[5]),
             "rho16":float(r[6]),"lo16":float(r[7]),"hi16":float(r[8])} for r in rows}

def rankpct(values):
    values=np.asarray(values,float); n=len(values)
    order=np.argsort(values,kind="mergesort"); ranks=np.empty(n)
    i=0
    while i<n:
        j=i
        while j+1<n and values[order[j+1]]==values[order[i]]: j+=1
        ranks[order[i:j+1]]=(i+j)/2.0
        i=j+1
    return ranks/(n-1) if n>1 else ranks*0

def short_species(species):
    genus,species_epithet=species.split("_",1)
    return f"{genus[0]}. {species_epithet.replace('_',' ')}"

def genus_species(species):
    genus,species_epithet=species.split("_",1)
    return f"{genus} {species_epithet.replace('_',' ')}"

plt.rcParams.update({
    "font.family":"Arial","font.size":14,"axes.titlesize":16,"axes.labelsize":13.5,
    "xtick.labelsize":12,"ytick.labelsize":12.5,"legend.fontsize":12,
    "axes.linewidth":1.0,"pdf.fonttype":42,"ps.fonttype":42,
    "figure.dpi":250,"savefig.dpi":250,"savefig.bbox":"tight",
})

fig=plt.figure(figsize=(15.0,12.0))
gs=fig.add_gridspec(2,1,height_ratios=[1.08,0.94],hspace=0.28)

# Panel A: all-species summary.
ax=fig.add_subplot(gs[0])
data=sorted(stat.items(),key=lambda item:item[1]["rho2"]-item[1]["rho16"])
y=np.arange(len(data)); height=0.38
rho2=np.array([item[1]["rho2"] for item in data]); lo2=np.array([item[1]["lo2"] for item in data]); hi2=np.array([item[1]["hi2"] for item in data])
rho16=np.array([item[1]["rho16"] for item in data]); lo16=np.array([item[1]["lo16"] for item in data]); hi16=np.array([item[1]["hi16"] for item in data])
ax.barh(y+height/2,rho2,height,color="#1F77B4",label="2bRAD",edgecolor="black",linewidth=.5,
        xerr=np.vstack([rho2-lo2,hi2-rho2]),error_kw=dict(ecolor="#0B3D66",elinewidth=1.3,capsize=3))
ax.barh(y-height/2,rho16,height,color="#D62728",label="16S rRNA",edgecolor="black",linewidth=.5,
        xerr=np.vstack([rho16-lo16,hi16-rho16]),error_kw=dict(ecolor="#7A1416",elinewidth=1.3,capsize=3))
ax.set_yticks(y)
ax.set_yticklabels([genus_species(species) for species,_ in data],style="italic")
ax.set_xlabel("Spearman ρ between strain and whole-genome distance")
ax.set_xlim(-.35,1.10); ax.set_ylim(-.75,len(data)-.25); ax.axvline(0,color=".55",linewidth=.9)
median2=np.median(rho2); median16=np.median(rho16)
ax.axvline(median2,color="#1F77B4",ls="--",linewidth=1.6)
ax.axvline(median16,color="#D62728",ls="--",linewidth=1.6)
ax.text(.97,.95,f"2bRAD median = {median2:.2f}\n16S median = {median16:.2f}",
        transform=ax.transAxes,ha="right",va="top",fontsize=13.5,
        bbox=dict(boxstyle="round,pad=.38",facecolor="white",edgecolor=".75"))
ax.legend(loc="lower right",frameon=True,framealpha=.95)
ax.grid(alpha=.22,axis="x")
ax.text(-.065,1.025,"A",transform=ax.transAxes,fontsize=21,fontweight="bold",ha="right",va="bottom")
ax.set_title("Complete/near-complete genome panels; error bars show 95% CIs",loc="left",pad=10)

# Panel B: lowest-, median- and highest-performing 16S examples.
files=sorted(PD/f"{species}.tsv" for species in stat)
files=[path for path in files if path.exists()]
order=sorted(files,key=lambda path:stat[path.stem]["rho16"])
selected=[order[0],order[len(order)//2],order[-1]]
inner=gs[1].subgridspec(1,3,wspace=.24)
for k,path in enumerate(selected):
    species=path.stem; values=stat[species]
    ax=fig.add_subplot(inner[k])
    pairs=[line.rstrip("\n").split("\t") for line in open(path)][1:]
    wgs=rankpct([float(row[2]) for row in pairs])
    b2=rankpct([float(row[3]) for row in pairs])
    s16=rankpct([float(row[4]) for row in pairs])
    ax.scatter(wgs,b2,s=20,color="#1F77B4",alpha=.36,edgecolors="none",label="2bRAD")
    ax.scatter(wgs,s16,s=20,color="#D62728",alpha=.36,edgecolors="none",label="16S rRNA")
    ax.plot([0,1],[0,1],color=".45",ls="--",linewidth=1.1)
    ax.set_title(f"{short_species(species)}\nρ {values['rho2']:.2f} / {values['rho16']:.2f}\nn pairs={len(pairs)}",
                 fontsize=13.5,pad=9)
    ax.set_xlim(-.04,1.04); ax.set_ylim(-.04,1.04)
    ax.set_xticks([0,.5,1]); ax.set_yticks([0,.5,1])
    ax.tick_params(labelsize=11.5); ax.grid(alpha=.18)
    if k==0: ax.set_ylabel("Marker-distance rank")
    if k==1: ax.set_xlabel("Whole-genome-distance rank")
    if k==0:
        ax.legend(loc="upper left",frameon=True,framealpha=.92,handletextpad=.2,borderpad=.4)
    if k==2:
        ax.text(.03,.95,"B",transform=ax.transAxes,fontsize=21,fontweight="bold",ha="left",va="bottom")

fig.savefig(f"{OUT}.png",bbox_inches="tight",facecolor="white")
fig.savefig(f"{OUT}.pdf",bbox_inches="tight",facecolor="white")
print(f"2bRAD median {median2:.3f}; 16S median {median16:.3f}; representative 16S ranks: {[path.stem for path in selected]}")
