#!/usr/bin/env python3
"""Simplified Strain2bScan workflow overview."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

BLUE_F,BLUE_E="#e7f0f9","#1f77b4"; DB_F="#cfe0f2"
GREEN_F,GREEN_E="#e4f1e8","#2e8b57"
AMB_F,AMB_E="#fdf0dc","#c8791a"; TXT="#12354f"; ARR="#5a6b7a"

def box(ax,x,y,w,h,text,fc=BLUE_F,ec=BLUE_E,fs=12,weight="normal"):
 ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.15,rounding_size=1.2",fc=fc,ec=ec,lw=1.6,zorder=3))
 ax.text(x,y,text,ha="center",va="center",fontsize=fs,color=TXT,weight=weight,zorder=5,linespacing=1.25)

def arrow(ax,x1,y1,x2,y2,label=None):
 ax.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle="-|>",color=ARR,lw=2,shrinkA=3,shrinkB=3),zorder=2)
 if label: ax.text((x1+x2)/2,(y1+y2)/2,label,ha="center",va="center",fontsize=10,color=ARR,style="italic",bbox=dict(boxstyle="round,pad=.18",fc="white",ec="none",alpha=.95),zorder=6)

fig,axs=plt.subplots(1,2,figsize=(13.5,7.0),gridspec_kw={"width_ratios":[1.45,1]})
ax=axs[0]; ax.set_xlim(0,100); ax.set_ylim(-4,106); ax.axis("off")
ax.set_title("A  Database construction and profiling",loc="left",fontsize=17,fontweight="bold",color=BLUE_E,pad=12)

box(ax,50,88,50,8,"Reference genomes per species",weight="bold")
arrow(ax,50,84,50,79,"type-IIB digestion")
box(ax,50,73,50,9,"Single-copy 2bRAD tags\n25–33 bp; canonical 64-bit markers")
arrow(ax,50,68.5,50,63.5)
box(ax,50,58,50,8,"Within-species 0.95 Jaccard clusters")
arrow(ax,50,54,50,49,"cluster-specific markers")
box(ax,50,44,50,8,"Cluster × marker database",fc=DB_F,weight="bold")
arrow(ax,50,40,50,35,"one digest / sample")
box(ax,50,30,50,9,"Shotgun reads → in-silico digest\nNative BcgI 2bRAD reads",fc=GREEN_F,ec=GREEN_E)
arrow(ax,50,25.5,50,20.5)
box(ax,50,16,50,8,"Species gate and cross-species marker restriction",fc=AMB_F,ec=AMB_E)
arrow(ax,50,12,50,7)
box(ax,50,3.5,50,7.5,"Cluster detection, depth–breadth filter, abundance",fc=AMB_F,ec=AMB_E)

ax=axs[1]
ax.set_title("B  Depth–breadth quality control",loc="left",fontsize=17,fontweight="bold",color=AMB_E,pad=12)
lam=np.linspace(.02,9,500)
ax.plot(lam,1-np.exp(-lam),color=BLUE_E,lw=3)
ax.scatter([.44],[.392],s=95,color=GREEN_E,zorder=5,edgecolor="white",lw=1.5)
ax.annotate("Genuine low-depth strain\n0.44×, breadth 0.39",(.44,.392),(1.2,.17),fontsize=12,color=GREEN_E,
            arrowprops=dict(arrowstyle="-",color=GREEN_E,lw=1.2))
ax.scatter([7.68],[.350],s=110,marker="X",color="#c0392b",zorder=5,edgecolor="white",lw=1.5)
ax.annotate("Partial-marker shadow\n7.7×, breadth 0.35",(7.68,.350),(4.2,.42),fontsize=12,color="#c0392b",
            arrowprops=dict(arrowstyle="-",color="#c0392b",lw=1.2))
ax.plot([7.68,7.68],[.350,1-np.exp(-7.68)],color="#c0392b",lw=1.4,ls="--")
ax.text(7.35,.68,"rejected",color="#c0392b",rotation=90,ha="right",va="center",fontsize=12)
ax.set_xlim(0,9.2); ax.set_ylim(0,1.08)
ax.set_xlabel("Estimated depth λ (reads / tag)",fontsize=13)
ax.set_ylabel("Observed breadth",fontsize=13)
ax.tick_params(labelsize=11.5)
for side in ("top","right"): ax.spines[side].set_visible(False)

fig.savefig("figures/source_panels/overview.png",dpi=250,facecolor="white")
fig.savefig("figures/source_panels/overview.pdf",facecolor="white")
print("wrote figures/source_panels/overview.png + .pdf")
