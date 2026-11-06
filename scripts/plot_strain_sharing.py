#!/usr/bin/env python3
"""Plot Strain2bScan strain-sharing rates for saliva and monkey cohorts."""
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/"results/strain_sharing_20261105"
mpl.rcParams.update({"font.family":"Arial","font.size":7.5,"axes.titlesize":8.5,
                     "axes.labelsize":8.0,"xtick.labelsize":7.0,"ytick.labelsize":7.0,
                     "pdf.fonttype":42,"figure.dpi":300})
sal_long=pd.read_csv(IN/"saliva_pairwise_strain_sharing.tsv",sep="\t")
sal_sum=pd.read_csv(IN/"saliva_strain_sharing_group_summary.tsv",sep="\t")
mon_sum=pd.read_csv(IN/"monkey_strain_sharing_group_summary.tsv",sep="\t")
age=pd.read_csv(IN/"monkey_between_age_corresponding_site_summary.tsv",sep="\t")
fig,axs=plt.subplots(2,2,figsize=(11.4,7.8)); axs=axs.ravel()

# A saliva dominant sharing
order=["same_subject_between_timepoints","same_timepoint_between_subjects","different_subject_between_timepoints"]
labels=["Same subject\nbetween timepoints","Different subject\nsame timepoint","Different subject\ndifferent timepoints"]
data=[sal_long.loc[sal_long.comparison_group==x,"strain_sharing_rate_dominant_cluster"] for x in order]
ax=axs[0]; bp=ax.boxplot(data,tick_labels=labels,patch_artist=True,widths=.55,showfliers=False,medianprops=dict(color="black",lw=1.0))
for b in bp["boxes"]: b.set_facecolor("#C6DBEF"); b.set_edgecolor("#2171B5"); b.set_linewidth(.7)
for i,d in enumerate(data): ax.scatter(np.random.default_rng(11+i).normal(i+1,.06,len(d)),d,s=7,color="#08519C",alpha=.45,linewidths=0)
ax.set_ylabel("Dominant-cluster sharing rate"); ax.set_title("A  Saliva: within-subject versus between-subject",loc="left",fontweight="bold")

# B saliva any sharing
data=[sal_long.loc[sal_long.comparison_group==x,"strain_sharing_rate_any_cluster"] for x in order]
ax=axs[1]; bp=ax.boxplot(data,tick_labels=labels,patch_artist=True,widths=.55,showfliers=False,medianprops=dict(color="black",lw=1.0))
for b in bp["boxes"]: b.set_facecolor("#FC8D59"); b.set_edgecolor("#D7301F"); b.set_linewidth(.7)
for i,d in enumerate(data): ax.scatter(np.random.default_rng(21+i).normal(i+1,.06,len(d)),d,s=7,color="#7F0000",alpha=.4,linewidths=0)
ax.set_ylabel("Any-shared-cluster rate"); ax.set_title("B  Saliva: coexistence-aware sharing",loc="left",fontweight="bold")

# C monkey selected comparisons, dominant
ax=axs[2]
sel=["within_monkey_between_position_same_region","within_monkey_mucosal_luminal_same_region",
     "between_monkey_same_age_same_region_compartment","between_age_corresponding_region_compartment",
     "between_monkey_between_region_same_compartment"]
labels=["Within monkey\nbetween positions","Within monkey\nmucosal/luminal",
        "Between monkeys\nsame age/site","Old/young\ncorresponding site","Between monkeys\nbetween regions"]
plot=mon_sum.set_index("comparison_group").loc[sel]
x=np.arange(len(plot))
ax.bar(x,plot.strain_sharing_rate_dominant_cluster_mean,color="#4393C3",width=.58,
       yerr=plot.strain_sharing_rate_dominant_cluster_mean.plot_error if False else None,
       capsize=2,edgecolor="white",linewidth=.4)
for xx,v in zip(x,plot.strain_sharing_rate_dominant_cluster_mean): ax.text(xx,v+.012,f"{v:.3f}",ha="center",fontsize=6.2)
ax.set_xticks(x); ax.set_xticklabels(labels,fontsize=6.2); ax.set_ylabel("Mean dominant-cluster sharing rate")
ax.set_ylim(0,.48); ax.set_title("C  Monkey GI spatial and age contrasts",loc="left",fontweight="bold")

# D old/young corresponding-site heatmap: dominant then any
ax=axs[3]
vis=age[age.n_pairs_analyzable>0].copy()
mat=vis[["strain_sharing_rate_dominant_cluster_mean","strain_sharing_rate_any_cluster_mean"]].to_numpy()
im=ax.imshow(mat,cmap="RdYlBu_r",vmin=0,vmax=1,aspect="auto")
ax.set_xticks([0,1]); ax.set_xticklabels(["Dominant cluster","Any shared cluster"])
ax.set_yticks(np.arange(len(vis)))
ax.set_yticklabels([f"{r.region} {r.compartment}\n(n={int(r.n_pairs)})" for _,r in vis.iterrows()],fontsize=6.2)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        ax.text(j,i,f"{mat[i,j]:.2f}",ha="center",va="center",fontsize=6.4,
                color="black" if .25<mat[i,j]<.75 else "white")
cb=fig.colorbar(im,ax=ax,pad=.02,aspect=18); cb.set_label("Old/young sharing rate",fontsize=6.5)
cb.ax.tick_params(labelsize=6); cb.outline.set_linewidth(.4)
ax.set_title("D  Old/young corresponding GI sites",loc="left",fontweight="bold")

fig.suptitle("Strain2bScan strain-sharing rates across timepoints, monkeys and GI sites",x=.012,ha="left",fontsize=10.5,fontweight="bold")
fig.tight_layout(rect=[0,.01,.99,.96])
for ext in ["png","pdf"]:
    fig.savefig(IN/f"strain_sharing_summary.{ext}",bbox_inches="tight",facecolor="white")
print(IN/"strain_sharing_summary.png")
