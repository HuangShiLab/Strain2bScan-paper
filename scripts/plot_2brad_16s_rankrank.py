#!/usr/bin/env python3
"""Supplement: per-species rank-rank comparison of 2bRAD and 16S distances."""
import os,glob,numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PAPER=os.path.expanduser("~/Downloads/Strain2bScan-paper"); PD=f"{PAPER}/results/pairdist"
TSV=f"{PAPER}/results/mash_2brad_vs_16s_recomputed.tsv"
rows=[l.rstrip("\n").split("\t") for l in open(TSV) if l.strip() and not l.startswith("#") and not l.startswith("species")]
stat={r[0]:(int(r[1]),float(r[3]),float(r[6])) for r in rows}
def rankpct(v):
 v=np.asarray(v,float); n=len(v); order=np.argsort(v,kind="mergesort"); rk=np.empty(n); i=0
 while i<n:
  j=i
  while j+1<n and v[order[j+1]]==v[order[i]]: j+=1
  rk[order[i:j+1]]=(i+j)/2; i=j+1
 return rk/(n-1) if n>1 else rk*.0
fig,axs=plt.subplots(3,5,figsize=(17,11)); axs=axs.ravel()
files=sorted(glob.glob(f"{PD}/*.tsv")); files.sort(key=lambda p:stat.get(os.path.basename(p)[:-4],(0,0,1))[2])
for k,p in enumerate(files):
 sp=os.path.basename(p)[:-4]; pr=[l.rstrip("\n").split("\t") for l in open(p)][1:]
 ax=axs[k]
 if pr:
  wg=rankpct([float(r[2]) for r in pr]); b2=rankpct([float(r[3]) for r in pr]); s16=rankpct([float(r[4]) for r in pr])
  ax.scatter(wg,b2,s=12,color="#1f77b4",alpha=.35,label="2bRAD")
  ax.scatter(wg,s16,s=12,color="#d62728",alpha=.35,label="16S")
  ax.plot([0,1],[0,1],color="#777",lw=1,ls="--")
 n,r2,r16=stat.get(sp,(0,0,0)); g,s=sp.split("_",1)
 ax.set_title(f"$\\it{{{g}.\\ {s.replace('_',' ')}}}$\\nρ 2bRAD {r2:.2f} / 16S {r16:.2f}",fontsize=14)
 ax.set_xlim(-.04,1.04); ax.set_ylim(-.04,1.04); ax.set_xticks([0,.5,1]); ax.set_yticks([0,.5,1]); ax.tick_params(labelsize=11); ax.grid(alpha=.18)
 if k%5==0: ax.set_ylabel("2bRAD / 16S rank percentile",fontsize=13)
 if k>=10: ax.set_xlabel("Whole-genome rank percentile",fontsize=13)
handles=[plt.Line2D([0],[0],marker="o",color="w",markerfacecolor="#1f77b4",label="2bRAD",markersize=10),plt.Line2D([0],[0],marker="o",color="w",markerfacecolor="#d62728",label="16S rRNA",markersize=10)]
fig.legend(handles=handles,loc="upper right",ncol=2,fontsize=14,frameon=True)
fig.tight_layout(rect=[0,0,.99,.96])
fig.savefig(f"{PAPER}/figures/source_panels/rankrank_2brad_16s.png",dpi=220)
fig.savefig(f"{PAPER}/figures/source_panels/rankrank_2brad_16s.pdf")
