#!/usr/bin/env python3
"""Two figures: Part I (native BcgI 2bRAD) and Part II (in-silico shotgun) strain-ID + abundance on the
ATCC MSA-1002 (even) / MSA-1003 (staggered) DNA mocks, combined 120-genome tree, --min-abundance 0."""
import csv, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
SP="/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
RAW="/Users/macstudio/Downloads/Strain2bScan-raw-data"; OUT=f"{RAW}/MSA_demo"
S2B,GREY="#1f77b4","#8c8c8c"
size={r["species"]:int(r["length_bp"]) for r in csv.DictReader(open(f"{OUT}/genome_qc_120.tsv"),delimiter="\t") if r["role"].startswith("ATCC")}

def atcc_map(tm):  # cluster -> species  (ATCC clusters only)
    return {l.split("\t")[1].strip():l.split("\t")[0] for l in open(tm) if not l.startswith("species")}

def depths(pred, cl2sp, mincov=0.2):  # species -> depth (ATCC clusters passing coverage)
    d={}
    for l in open(pred):
        if not l.startswith("C"): continue
        r=l.rstrip().split("\t"); cl=r[0]
        if cl in cl2sp and float(r[2])>=mincov: d[cl2sp[cl]]=float(r[4])
    return d
def mass_ra(d):
    m={s:d[s]*size[s] for s in d}; t=sum(m.values()) or 1; return {s:100*m[s]/t for s in m}

# ---------------- PART I : native BcgI 2bRAD ----------------
cl2sp_b=atcc_map(f"{OUT}/truth_map_bcgi_cont.tsv")
hosts=["90","95","99","99.9"]; host_det=[19,19,19,11]
ng=["0.001","0.01","0.1","1"]; ng_det=[0,16,19,19]
evenI=mass_ra(depths(f"{SP}/msa_demo2/cont_BcgI_MSA1002_90_100ng_1.pred",cl2sp_b))
stagI=mass_ra(depths(f"{SP}/msa_demo2/cont_BcgI_MSA1003_0_100ng_1.pred",cl2sp_b))

fig,ax=plt.subplots(1,3,figsize=(15,4.4))
# A host robustness
ax[0].bar(range(4),host_det,color=S2B,width=0.62)
for i,v in enumerate(host_det): ax[0].text(i,v+0.3,f"{v}/20",ha="center",fontsize=9)
ax[0].set_xticks(range(4)); ax[0].set_xticklabels([f"{h}%" for h in hosts]); ax[0].set_ylim(0,21)
ax[0].set_xlabel("human DNA fraction"); ax[0].set_ylabel("strains resolved (/20)")
ax[0].set_title("A  Host robustness (native 2bRAD, 100 ng)"); ax[0].grid(alpha=.3,axis="y")
# B titration
ax[1].plot([float(x) for x in ng],ng_det,"o-",color=S2B,lw=2,ms=8)
for x,v in zip([float(x) for x in ng],ng_det): ax[1].annotate(f"{v}/20",(x,v),textcoords="offset points",xytext=(0,8),ha="center",fontsize=9)
ax[1].set_xscale("log"); ax[1].set_xticks([float(x) for x in ng]); ax[1].set_xticklabels(ng)
ax[1].set_xlabel("2bRAD DNA input (ng)"); ax[1].set_ylabel("strains resolved (/20)"); ax[1].set_ylim(-0.5,21)
ax[1].set_title("B  Low-biomass titration"); ax[1].grid(alpha=.3)
# C abundance even vs staggered
spx=sorted(size, key=lambda s:-stagI.get(s,0))
y=np.arange(len(spx))
ax[2].barh(y-0.2,[evenI.get(s,0) for s in spx],0.4,color=GREY,label="MSA1002 (even)")
ax[2].barh(y+0.2,[stagI.get(s,1e-3) for s in spx],0.4,color=S2B,label="MSA1003 (staggered)")
ax[2].set_xscale("log"); ax[2].set_yticks(y); ax[2].set_yticklabels([s.split("_")[0][0]+". "+s.split("_")[1][:9] for s in spx],fontsize=6,style="italic")
ax[2].invert_yaxis(); ax[2].set_xlabel("species abundance (mass %, from depth, log)"); ax[2].set_title("C  Abundance: even vs staggered")
ax[2].legend(fontsize=7,loc="lower right"); ax[2].grid(alpha=.3,axis="x")
fig.suptitle("Part I — native BcgI 2bRAD: strain identification + abundance on ATCC MSA-1002/1003 DNA mock (combined 120-genome tree)",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(f"{OUT}/partI_native2brad.png",dpi=150,bbox_inches="tight"); fig.savefig(f"{OUT}/partI_native2brad.pdf",bbox_inches="tight")
print("wrote partI_native2brad.png")

# ---------------- PART II : in-silico shotgun ----------------
cl2sp_a=atcc_map(f"{OUT}/truth_map_all.tsv")
evenII=mass_ra(depths(f"{SP}/msa_prof/combined_WMS_MSA1002_0_100ng_2.pred",cl2sp_a))
stagII=mass_ra(depths(f"{SP}/msa_prof/combined_WMS_MSA1003_0_100ng_1.pred",cl2sp_a))
fig2,ax2=plt.subplots(1,2,figsize=(11,4.6))
# A detection
ax2[0].bar([0,1],[20,20],color=S2B,width=0.5)
for i,l in enumerate(["MSA1002\n(even)","MSA1003\n(staggered)"]):
    ax2[0].text(i,20.3,"20/20",ha="center",fontsize=11)
ax2[0].set_xticks([0,1]); ax2[0].set_xticklabels(["MSA1002 (even)","MSA1003 (staggered)"]); ax2[0].set_ylim(0,22)
ax2[0].set_ylabel("strains resolved (/20)"); ax2[0].set_title("A  Strain identification (0 false positives)"); ax2[0].grid(alpha=.3,axis="y")
# B abundance
spx=sorted(size, key=lambda s:-stagII.get(s,0)); y=np.arange(len(spx))
ax2[1].barh(y-0.2,[evenII.get(s,0) for s in spx],0.4,color=GREY,label="MSA1002 (even)")
ax2[1].barh(y+0.2,[stagII.get(s,1e-3) for s in spx],0.4,color=S2B,label="MSA1003 (staggered)")
ax2[1].set_xscale("log"); ax2[1].set_yticks(y); ax2[1].set_yticklabels([s.split("_")[0][0]+". "+s.split("_")[1][:9] for s in spx],fontsize=6,style="italic")
ax2[1].invert_yaxis(); ax2[1].set_xlabel("species abundance (mass %, from depth, log)"); ax2[1].set_title("B  Abundance: even vs staggered")
ax2[1].legend(fontsize=7,loc="lower right"); ax2[1].grid(alpha=.3,axis="x")
fig2.suptitle("Part II — in-silico shotgun (WMS): strain identification + abundance on ATCC MSA-1002/1003 DNA mock (combined 120-genome tree)",fontsize=11)
fig2.tight_layout(rect=[0,0,1,0.95]); fig2.savefig(f"{OUT}/partII_shotgun.png",dpi=150,bbox_inches="tight"); fig2.savefig(f"{OUT}/partII_shotgun.pdf",bbox_inches="tight")
print("wrote partII_shotgun.png")
