#!/usr/bin/env python3
"""Build Figure 11: rhesus macaque GI strain-resolved cohort analysis."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import textwrap
from matplotlib.patches import Patch

ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/"results/monkey_spatial_148/hpc_outputs/analysis"
SHARING=ROOT/"results/strain_sharing_20261105"
CURATED=ROOT/"results/monkey_spatial_148"
OUT=ROOT/"figures/numbered"
OUT.mkdir(parents=True,exist_ok=True)

mpl.rcParams.update({
    "font.family":"Arial","font.size":16.0,"axes.titlesize":18.5,"axes.labelsize":15.5,
    "xtick.labelsize":13.5,"ytick.labelsize":13.5,"legend.fontsize":12.5,
    "axes.linewidth":0.9,"pdf.fonttype":42,"ps.fonttype":42,
    "figure.dpi":300,"savefig.dpi":300,"savefig.bbox":"tight",
})
meta=pd.read_csv(IN/"full_148_20261027/full_148_sample_metadata.tsv",sep="\t")
species=pd.read_csv(IN/"full_148_20261027/full148_species_GI_control_summary.tsv",sep="\t")
clusters=pd.read_csv(IN/"full_148_20261027/full148_2brad_cluster_GI_control_region_summary.tsv",sep="\t")
paired_tests=pd.read_csv(IN/"full_148_20261027/full148_mucosal_luminal_paired_tests.tsv",sep="\t")
paired_tests=paired_tests[paired_tests.p_bonferroni<0.05].copy()
concordance=pd.read_csv(IN/"strain_profiling_20261023/paired_luminal_modality_concordance.tsv",sep="\t")
ba_coupling=pd.read_csv(CURATED/"bileacid_coupling_key_results.tsv",sep="\t")
func_coupling=pd.read_csv(CURATED/"function_coupling_key_results.tsv",sep="\t")
sharing=pd.read_csv(SHARING/"monkey_strain_sharing_group_summary.tsv",sep="\t")
age_sharing=pd.read_csv(SHARING/"monkey_between_age_corresponding_site_summary.tsv",sep="\t")

def clean_species(x):
    return str(x).replace("_sp0"," sp0").replace("_sp9"," sp9").replace("_A"," A").replace("_B"," B")

def cluster_label(row,maxlen=38):
    if "species" in row:
        species=str(row["species"]); cluster=str(row.get("cluster",""))
    else:
        combined=str(row.get("strain_cluster","")); species,cluster=combined.rsplit("|",1)
    species=species.replace("_sp"," sp").replace("_A"," A").replace("_B"," B")
    label=f"{species} {cluster}"
    return label if len(label)<=maxlen else label[:maxlen-1]+"…"

def coupling_label(row,maxlen=58):
    feature=row.get("feature",row.get("pathway",""))
    short={
      "conjugated_or_highmass":"Conjugated BA",
      "cholate_like_Q1_407":"Cholate-like BA",
      "dihydroxy_Q1_391":"Dihydroxy BA",
      "LCA_like_Q1_375":"LCA-like BA",
      "CENTFERM-PWY: pyruvate fermentation to butanoate":"Pyruvate → butanoate",
      "PWY-5676: acetyl-CoA fermentation to butanoate II":"Acetyl-CoA → butanoate II",
      "Oribacterium_sp900768305|C2_PYRUVTETHSYN-PWY: pyruvate fermentation to butanoate":"Pyruvate → butanoate",
      "P108-PWY: pyruvate fermentation to propanoate I":"Pyruvate → propanoate",
      "PWY-5677: succinate fermentation to butanoate":"Succinate → butanoate",
      "P441-PWY: superpathway of N-acetylneuraminate degradation":"N-acetylneuraminate degradation",
      "GLCMANNANAUT-PWY: superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation":"N-acetyl-sugar degradation",
    }.get(feature,feature)
    display=cluster_label(row,maxlen=40)
    for long_name,short_name in [
        ("Cryptobacteroides sp000433355","Cryptobact."),("Cryptobacteroides sp900546445","Cryptobact."),
        ("Treponema_D_succinifaciens","T._succinifaciens"),("Treponema_D_berlinense","T._berlinense"),
        ("Succinivibrio sp000431835","Succinivibrio"),("UBA4372 sp900543815","UBA4372"),
        ("UBA636 sp900546285","UBA636"),("Oribacterium sp900768305","Oribacterium"),
    ]: display=display.replace(long_name,short_name)
    label=textwrap.fill(f"{display} — {short}",width=maxlen,break_long_words=False)
    return label if len(label)<=maxlen else label[:maxlen-1]+"…"

def coupling_ax(ax,data,title):
    plot=data.sort_values("rho").copy()
    plot["label"]=plot.apply(coupling_label,axis=1)
    y=np.arange(len(plot)); colors=["#B2182B" if v<0 else "#2166AC" for v in plot.rho]
    ax.barh(y,plot.rho,color=colors,height=.62)
    ax.axvline(0,color="0.35",lw=.8)
    ax.set_yticks(y); ax.set_yticklabels(plot.label,fontsize=13.5)
    ax.set_xlabel("Spearman ρ with strain-cluster composition")
    ax.set_xlim(-1.05,1.05); ax.set_title(title,loc="left",fontweight="bold")

fig=plt.figure(figsize=(21.0,18.5),layout="constrained")
gs=fig.add_gridspec(3,3,height_ratios=[.82,.92,.92])

# A
ax=fig.add_subplot(gs[0,0])
counts=meta[meta.sample_group=="GI_spatial"].groupby(["monkey_id","gi_region","compartment"]).size().unstack(fill_value=0)
order=[(m,r) for m in ["O1","O2","Y1","Y2"] for r in ["E","S","SI","LI","K"]]
counts=counts.reindex(pd.MultiIndex.from_tuples(order,names=["monkey_id","gi_region"]),fill_value=0)
x=np.arange(len(counts)); width=.78
ax.bar(x,counts.get("luminal",0),width,label="Luminal",color="#4C78A8",edgecolor="white",lw=.4)
ax.bar(x,counts.get("mucosal",0),width,bottom=counts.get("luminal",0),label="Mucosal",color="#F58518",edgecolor="white",lw=.4)
for xx in [4.5,9.5,14.5]: ax.axvline(xx,color="0.75",lw=.7,ls=":")
ax.set_xticks(x); ax.set_xticklabels([r for _,r in counts.index],fontsize=11.0,rotation=45,ha="right")
for center,label in [(2,"O1"),(7,"O2"),(12,"Y1"),(17,"Y2")]:
    ax.text(center,-.18,label,transform=ax.get_xaxis_transform(),ha="center",va="top",fontsize=12.5,fontweight="bold")
ax.set_ylabel("GI libraries"); ax.set_ylim(0,max(counts.sum(axis=1))*1.30)
ax.set_title("A  Cohort design",loc="left",fontweight="bold")
ax.legend(frameon=True,framealpha=.95,edgecolor=".75",ncol=2,loc="upper center",bbox_to_anchor=(.5,-.42))

# B
ax=fig.add_subplot(gs[0,1])
top=species.sort_values("prevalence_GI",ascending=False).head(10).iloc[::-1]
y=np.arange(len(top)); h=.18
for j,(col,label,color) in enumerate([("prevalence_GI","GI","#2F4B7C"),("prevalence_luminal","Luminal","#4C78A8"),("prevalence_mucosal","Mucosal","#F58518"),("prevalence_body_control","Body control","#B0B0B0")]):
    ax.barh(y+(j-1.5)*h,top[col],height=h,label=label,color=color)
ax.set_yticks(y); ax.set_yticklabels([clean_species(x) for x in top.species],fontsize=11.5)
ax.set_xlabel("Species prevalence"); ax.set_xlim(0,1); ax.xaxis.set_major_formatter(lambda v,_:f"{v:.0%}")
ax.set_title("B  Species recovery",loc="left",fontweight="bold")
ax.legend(frameon=False,ncol=2,loc="upper center",bbox_to_anchor=(.5,-.30))

# C
ax=fig.add_subplot(gs[0,2])
li=clusters[(clusters.peak_region=="LI")&(clusters.region_specificity>=.8)&(clusters.max_region_prevalence>=.3)].sort_values(["max_region_prevalence","region_specificity"],ascending=False).head(10).iloc[::-1]
y=np.arange(len(li))
for i,(_,r) in enumerate(li.iterrows()): ax.plot([r.prevalence_mucosal,r.prevalence_luminal],[i,i],color="#D0D0D0",lw=1)
ax.scatter(li.prevalence_mucosal,y,s=36,color="#F58518",label="Mucosal")
ax.scatter(li.prevalence_luminal,y,s=40,color="#4C78A8",label="Luminal")
ax.scatter(li.prevalence_LI,y,s=48,facecolor="none",edgecolor="#B2182B",lw=1.2,label="Peak LI")
panel_labels=[f"{clean_species(r.species)} {r.cluster}" for _,r in li.iterrows()]
ax.set_yticks(y); ax.set_yticklabels(panel_labels,fontsize=11.5)
ax.set_xlabel("Cluster prevalence"); ax.set_xlim(0,.85); ax.xaxis.set_major_formatter(lambda v,_:f"{v:.0%}")
ax.set_title("C  LI-enriched strain clusters (89 total)",loc="left",fontweight="bold")
ax.legend(frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(.5,-.30),columnspacing=.8,handletextpad=.2)

# D
ax=fig.add_subplot(gs[1,0])
plot=paired_tests.sort_values("p_bonferroni",ascending=False).copy(); plot["minus_log_q"]=-np.log10(plot.p_bonferroni)
y=np.arange(len(plot)); ax.barh(y,plot.minus_log_q,color="#2166AC",height=.62)
ax.set_yticks(y); ax.set_yticklabels([f"{clean_species(r.species)} {r.cluster}" for _,r in plot.iterrows()],fontsize=11.5)
ax.set_xlabel(r"$-\log_{10}$(Bonferroni $q$)"); ax.set_xlim(0,max(plot.minus_log_q)*1.36)
for i,(_,r) in enumerate(plot.iterrows()):
    ax.text(r.minus_log_q+.08,i,f"{int(r.n_luminal_positive)} L / {int(r.n_mucosal_positive)} M",va="center",fontsize=8.2,color=".2")
ax.set_title("D  Mucosal/luminal tests",loc="left",fontweight="bold")

# E
ax=fig.add_subplot(gs[1,1])
sc=ax.scatter(concordance.n_shotgun_calls,concordance.n_2brad_calls,c=concordance.jaccard_calls,cmap="viridis_r",s=40,edgecolor=".55",lw=.5,norm=Normalize(0,1))
lim=max(concordance.n_shotgun_calls.max(),concordance.n_2brad_calls.max())*1.08
ax.plot([0,lim],[0,lim],color=".45",lw=.7,ls="--"); ax.set_xlim(0,lim); ax.set_ylim(0,lim)
ax.set_xlabel("Shotgun strain-cluster calls"); ax.set_ylabel("2bRAD strain-cluster calls")
cb=fig.colorbar(sc,ax=ax,pad=.015,aspect=22); cb.set_label("Strain-call Jaccard",fontsize=9.5); cb.ax.tick_params(labelsize=8.5); cb.outline.set_linewidth(.6)
ax.set_title("E  Paired luminal calls",loc="left",fontweight="bold")

# F BA coupling
ax=fig.add_subplot(gs[1,2]); coupling_ax(ax,ba_coupling,"F  Bile-acid associations")

# G functional coupling
ax=fig.add_subplot(gs[2,0]); coupling_ax(ax,func_coupling,"G  Pathway associations")

# H sharing
ax=fig.add_subplot(gs[2,1])
sel=["within_monkey_between_position_same_region","within_monkey_mucosal_luminal_same_region","between_monkey_same_age_same_region_compartment","between_age_corresponding_region_compartment","between_monkey_between_region_same_compartment"]
plot=sharing.set_index("comparison_group").loc[sel]; x=np.arange(len(plot))
ax.bar(x,plot.strain_sharing_rate_dominant_cluster_mean,color="#4393C3",width=.58,edgecolor="white",lw=.4)
for xx,v in zip(x,plot.strain_sharing_rate_dominant_cluster_mean): ax.text(xx,v+.012,f"{v:.2f}",ha="center",va="bottom",fontsize=9,color=".15")
ax.set_xticks(x); ax.set_xticklabels(["1","2","3","4","5"],fontsize=10); ax.set_ylabel("Mean dominant-cluster sharing")
ax.set_ylim(0,.46); ax.set_title("H  Dominant-cluster sharing contrasts",loc="left",fontweight="bold")
labels=["Within monkey","Within mucosal/luminal","Same age/site","Old/young site","Between regions"]
handles=[Patch(facecolor="#4393C3",edgecolor="white",label=f"{i}  {label}") for i,label in enumerate(labels,1)]
ax.legend(handles=handles,frameon=False,ncol=2,fontsize=11.5,loc="upper center",bbox_to_anchor=(.5,-.24),columnspacing=.8,handlelength=1,handletextpad=.4)

# I age sharing
ax=fig.add_subplot(gs[2,2])
vis=age_sharing[age_sharing.n_pairs_analyzable>0].copy()
matrix=vis[["strain_sharing_rate_dominant_cluster_mean","strain_sharing_rate_any_cluster_mean"]].to_numpy()
im=ax.imshow(matrix,cmap="RdYlBu_r",vmin=0,vmax=1,aspect="auto")
ax.set_xticks([0,1]); ax.set_xticklabels(["Dominant\ncluster","Any shared\ncluster"],fontsize=9)
ax.set_yticks(np.arange(len(vis))); ax.set_yticklabels([f"{r.region} {r.compartment}\n(n={int(r.n_pairs)})" for _,r in vis.iterrows()],fontsize=11.5)
for i in range(matrix.shape[0]):
 for j in range(matrix.shape[1]):
  value=matrix[i,j]; ax.text(j,i,f"{value:.2f}",ha="center",va="center",fontsize=9,color="black" if .25<value<.75 else "white")
cb=fig.colorbar(im,ax=ax,pad=.02,aspect=18); cb.set_label("Old/young sharing rate",fontsize=9.5); cb.ax.tick_params(labelsize=8.5); cb.outline.set_linewidth(.6)
ax.set_title("I  Old/young corresponding GI sites",loc="left",fontweight="bold")

fig.suptitle("Rhesus macaque GI cohort: 128 GI samples and 18-species/303-cluster monkey MAG panel",x=.012,ha="left",fontsize=22,fontweight="bold")
for ext in ["png","pdf"]:
    fig.savefig(OUT/f"Fig11_monkey_spatial_cohort.{ext}",bbox_inches="tight",facecolor="white")
print(OUT/"Fig11_monkey_spatial_cohort.png")
