#!/usr/bin/env python3
"""Build Figure 8 computational-performance panels from the 225-sample simulation benchmark."""
from pathlib import Path
import csv, statistics
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import numpy as np

PAPER=Path(__file__).resolve().parents[1]
DEPTHS=[0.5,1,3,5,10]; COLORS={"Strain2bScan":"#1f77b4","StrainScan":"#6e6e6e","StrainGST":"#2ca02c","inStrain":"#cc79a7"}

def rows(path): return list(csv.DictReader(open(PAPER/path, newline=""), delimiter="\t"))
def median_by_depth(records, wall_field="wall_s", rss_field="maxrss_mb"):
 out={}
 for d in DEPTHS:
  rr=[r for r in records if float(r.get("depth",d))==d]
  if rr:
   out[d]=(statistics.median(float(r[wall_field]) for r in rr),statistics.median(float(r[rss_field]) for r in rr)/1024/1024)
 return out

same=rows("results/sim_headtohead_profile_sameenv.tsv")
s2b=[]; ss=[]
for r in same:
 d=float(r["depth"]); n=int(r["n"])
 s2b += [(d,)*n for _ in []]
 # archived table already contains depth medians; use those values directly
 s2b.append((d,float(r["s2b_emu_wall_s"]),float(r["s2b_emu_rss_mb"])))
 ss.append((d,float(r["ss_wall_s"]),float(r["ss_rss_mb"])))
sg=[(float(r["depth"]),float(r["total_wall_s"]),float(r["maxrss_bytes"])) for r in rows("results/strainge_rerun/simulation_run_efficiency.tsv") if r["kind"]=="single"]
try:
 isn=rows("results/instrain_sim225/profiling_runtime.tsv")
 isn_records=[(float(r["depth"]),float(r["wall_s"]),float(r["maxrss_mb"])*1024*1024) for r in isn]
except FileNotFoundError:
 isn_records=[]
sg_single=[(float(r["depth"]),float(r["total_wall_s"]),float(r["maxrss_bytes"])) for r in rows("results/strainge_rerun/simulation_run_efficiency.tsv") if r["kind"]=="single"]

def med(records,key):
 out=[]
 for d in DEPTHS:
  vals=[r[key] for r in records if r[0]==d]
  if vals: out.append((d,statistics.median(vals)))
 return out

fig,axs=plt.subplots(2,2,figsize=(14.0,10.0)); axs=axs.ravel()
# A wall time
ax=axs[0]
sg_med=med(sg_single,1); isn_med=med(isn_records,1)
sg_mem=med(sg_single,2); isn_mem=med(isn_records,2)
series={"Strain2bScan":[(d,v) for d,v,_ in s2b],"StrainScan":[(d,v) for d,v,_ in ss],
        "StrainGST":sg_med,"inStrain":isn_med}
styles={"Strain2bScan":"o-","StrainScan":"s--","StrainGST":"^-.","inStrain":"v:"}
for tool,xy in series.items():
 if xy: ax.plot([x[0] for x in xy],[x[1] for x in xy],styles[tool],color=COLORS[tool],lw=2.4,ms=8,label=tool)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xticks(DEPTHS); ax.set_xticklabels(DEPTHS)
ax.set_xlabel("per-strain depth (×)"); ax.set_ylabel("median wall time (s, log)")
ax.set_title("A  Per-sample wall time across 225-sample design",loc="left",fontweight="bold")
ax.grid(True,which="both",alpha=.25); ax.legend(fontsize=10)

# B RSS
ax=axs[1]
series={"Strain2bScan":[(d,m/1024) for d,_,m in s2b],"StrainScan":[(d,m/1024) for d,_,m in ss],
        "StrainGST":[(d,m/1024/1024) for d,m in sg_mem],"inStrain":[(d,m/1024/1024) for d,m in isn_mem]}
for tool,xy in series.items():
 if xy: ax.plot([x[0] for x in xy],[x[1] for x in xy],styles[tool],color=COLORS[tool],lw=2.4,ms=8,label=tool)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xticks(DEPTHS); ax.set_xticklabels(DEPTHS)
ax.set_xlabel("per-strain depth (×)"); ax.set_ylabel("median peak RSS (GB, log)")
ax.set_title("B  Per-sample peak memory",loc="left",fontweight="bold")
ax.grid(True,which="both",alpha=.25); ax.legend(fontsize=10)

# C database build by species
ax=axs[2]
build_s2b=rows("results/sim_headtohead_build_cost.tsv")
build_ss=build_s2b
build_sg={r["species"]:float(r["build_wall_s"]) for r in rows("results/strainge_rerun/simulation_database_build.tsv")}
species=[r["species"] for r in build_s2b]
y=np.arange(len(species)); h=.26
ax.barh(y-h,[float(r["s2b_build_s"]) for r in build_s2b],h,color=COLORS["Strain2bScan"],label="Strain2bScan")
ax.barh(y,[float(r["ss_build_s"]) if r["ss_build_s"]!="NA" else np.nan for r in build_ss],h,color=COLORS["StrainScan"],label="StrainScan")
ax.barh(y+h,[build_sg[s] for s in species],h,color=COLORS["StrainGST"],label="StrainGST")
ax.set_yticks(y); ax.set_yticklabels([s.replace("_"," ") for s in species],fontsize=9,style="italic")
ax.set_xscale("log"); ax.set_xlabel("database-build wall time (s, log)")
ax.set_title("C  Database build across 15 species",loc="left",fontweight="bold")
ax.grid(axis="x",alpha=.25); ax.legend(fontsize=9,loc="lower right")

# D multi-species
ax=axs[3]
multi=rows("results/sim_headtohead_multi.tsv")
ge_multi=[r for r in rows("results/strainge_rerun/simulation_run_efficiency.tsv") if r["kind"]=="multi"]
labels=[r["depth"].replace("depth_","")+"×" for r in multi]; x=np.arange(len(multi)); w=.24
ax.bar(x-w,[float(r["s2b_wall_s"]) for r in multi],w,color=COLORS["Strain2bScan"],label="Strain2bScan")
ax.bar(x,[float(r["ss_total_wall_s"]) for r in multi],w,color=COLORS["StrainScan"],label="StrainScan")
ax.bar(x+w,[statistics.median(float(r["total_wall_s"]) for r in ge_multi if r["depth"]==r2["depth"].replace("depth_","")) for r2 in multi],w,color=COLORS["StrainGST"],label="StrainGST")
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_xlabel("community depth"); ax.set_ylabel("wall time / community sample (s, log)")
ax.set_title("D  Multi-species throughput",loc="left",fontweight="bold")
ax.grid(axis="y",alpha=.25); ax.legend(fontsize=9)

fig.suptitle("Computational performance on the 15-species simulation benchmark",x=.012,ha="left",fontsize=17,fontweight="bold")
fig.tight_layout(rect=[0,0,.99,.96])
for ext in ["png","pdf"]: fig.savefig(PAPER/f"figures/source_panels/computational_performance.{ext}",bbox_inches="tight",facecolor="white")
print(PAPER/"figures/source_panels/computational_performance.png")
