#!/usr/bin/env python3
"""Figure 2: 2bRAD versus 16S genome-distance preservation (summary panel)."""
import os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

PAPER=os.path.expanduser("~/Downloads/Strain2bScan-paper")
TSV=f"{PAPER}/results/mash_2brad_vs_16s_recomputed.tsv"
rows=[l.rstrip("\n").split("\t") for l in open(TSV) if l.strip() and not l.startswith("#") and not l.startswith("species")]
stat={r[0]:(int(r[1]),float(r[3]),float(r[4]),float(r[5]),float(r[6]),float(r[7]),float(r[8])) for r in rows}
d=sorted([(sp,*vals) for sp,vals in stat.items()],key=lambda x:x[2]-x[5])
def split(sp): g,s=sp.split("_",1); return g,s.replace("_"," ")
fig,ax=plt.subplots(figsize=(9.2,10.4))
y=np.arange(len(d)); h=.38
labels=[f"{split(x[0])[0]} {split(x[0])[1]}\nn = {x[1]}" for x in d]
s2=np.array([x[2] for x in d]); e2=np.array([[x[2]-x[3] for x in d],[x[4]-x[2] for x in d]])
s16=np.array([x[5] for x in d]); e16=np.array([[x[5]-x[6] for x in d],[x[7]-x[5] for x in d]])
ax.barh(y+h/2,s2,h,color="#1f77b4",label="2bRAD",edgecolor="k",lw=.4,xerr=e2,error_kw=dict(ecolor="#0b3d66",elinewidth=1.3,capsize=3))
ax.barh(y-h/2,s16,h,color="#d62728",label="16S rRNA",edgecolor="k",lw=.4,xerr=e16,error_kw=dict(ecolor="#7a1416",elinewidth=1.3,capsize=3))
ax.set_yticks(y); ax.set_yticklabels(labels,fontsize=12,style="italic")
ax.set_xlabel("Spearman ρ between strain and whole-genome distance",fontsize=14)
ax.set_xlim(-.25,1.08); ax.set_ylim(-.75,len(d)-.25); ax.axvline(0,color="#888",lw=.9)
med2=np.median(s2); med16=np.median(s16)
ax.axvline(med2,color="#1f77b4",ls="--",lw=1.5); ax.axvline(med16,color="#d62728",ls="--",lw=1.5)
ax.text(.97,.96,f"2bRAD median = {med2:.2f}\n16S median = {med16:.2f}",transform=ax.transAxes,ha="right",va="top",fontsize=13,bbox=dict(boxstyle="round,pad=.4",facecolor="white",edgecolor=".75"))
ax.legend(loc="lower right",fontsize=12,framealpha=.95); ax.grid(alpha=.25,axis="x")
ax.tick_params(axis="x",labelsize=12)
ax.set_title("2bRAD preserves whole-genome strain distances;\n16S does not",loc="left",fontsize=17,fontweight="bold",pad=12)
fig.tight_layout()
fig.savefig(f"{PAPER}/figures/source_panels/mash_2brad_vs_16s.png",dpi=250)
fig.savefig(f"{PAPER}/figures/source_panels/mash_2brad_vs_16s.pdf")
print(f"2bRAD median {med2:.3f} vs 16S {med16:.3f}")
