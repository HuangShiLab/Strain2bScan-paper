#!/usr/bin/env python3
"""Permanova and leave-one-timepoint-out analysis of the full-depth WMS saliva rebuild."""
from __future__ import annotations
import csv, statistics
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"results/wms32_rebuild/wms32_strain_calls_all.tsv"
OUT=ROOT/"results/wms32_rebuild"
rng=np.random.default_rng(42)

def sample_label(raw): return raw.removeprefix("WMS_").replace("_","-")

rows=list(csv.DictReader(SRC.open(),delimiter="\t"))
samples=sorted({sample_label(r["#sample"]) for r in rows})
subject_map={s:s.split("-")[1] for s in samples}
timepoint_map={s:s.split("-")[0].removeprefix("S") for s in samples}
subject=np.array([subject_map[s] for s in samples])
timepoint=np.array([timepoint_map[s] for s in samples])
si={s:i for i,s in enumerate(samples)}
species=sorted({r["species"] for r in rows})
features=sorted({f'{r["species"]}|{r["cluster"]}' for r in rows})
sp_i={x:i for i,x in enumerate(species)}; ft_i={x:i for i,x in enumerate(features)}
Ms=np.zeros((len(samples),len(features))); Mp=np.zeros((len(samples),len(species)))
for r in rows:
    i=si[sample_label(r["#sample"])]; frac=float(r["sample_fraction"])
    Ms[i,ft_i[f'{r["species"]}|{r["cluster"]}']]+=frac; Mp[i,sp_i[r["species"]]]+=frac
for M in (Ms,Mp):
    total=M.sum(1,keepdims=True); total[total==0]=1; M/=total

def bray(M):
    n=len(M); D=np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            den=(M[i]+M[j]).sum(); D[i,j]=D[j,i]=abs(M[i]-M[j]).sum()/den if den else 0
    return D

def permanova(D,labels,perms=4999):
    labels=np.asarray(labels); N=len(labels); a=len(np.unique(labels))
    def ss(D2,l):
        total=(D2[np.triu_indices(N,1)].sum()/N); within=0.0
        for g in np.unique(l):
            idx=np.where(l==g)[0]; ng=len(idx)
            if ng>1: within += D2[np.ix_(idx,idx)][np.triu_indices(ng,1)].sum()/ng
        return total,within
    D2=D*D; total,w=ss(D2,labels); F=((total-w)/(a-1))/(w/(N-a)); ge=1
    for _ in range(perms):
        l=rng.permutation(labels); t2,w2=ss(D2,l); Fp=((t2-w2)/(a-1))/(w2/(N-a)); ge+=Fp>=F
    return (total-w)/total,F,ge/(perms+1)

def loto(D,labels,groups):
    labels=np.asarray(labels); groups=np.asarray(groups); correct=0
    for held in sorted(set(groups)):
        test=np.where(groups==held)[0]; train=np.where(groups!=held)[0]
        for i in test: correct += labels[train[np.argmin(D[i,train])]]==labels[i]
    return correct/len(labels)

def write(path,fields,records):
    with path.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(records)

results=[]
for level,M in (("species",Mp),("strain",Ms)):
    D=bray(M)
    r2s,fp,pp=permanova(D,subject); r2t,ft,pt=permanova(D,timepoint)
    results.append({"factor":"subject","level":level,"n_samples":len(samples),"n_features":M.shape[1],
                    "R2":r2s,"pseudo_F":fp,"p_value":pp,"LOO_1NN_acc":loto(D,subject,timepoint)})
    results.append({"factor":"timepoint","level":level,"n_samples":len(samples),"n_features":M.shape[1],
                    "R2":r2t,"pseudo_F":ft,"p_value":pt,"LOO_1NN_acc":loto(D,timepoint,subject)})
write(OUT/"wms32_temporal_ml.tsv",["factor","level","n_samples","n_features","R2","pseudo_F","p_value","LOO_1NN_acc"],results)

# Per-species within-species strain structure.
Ds={}
perspecies=[]
for sp in species:
    idx=[i for i,r in enumerate(rows) if r["species"]==sp]
    sp_samples=sorted({sample_label(rows[i]["#sample"]) for i in idx})
    if len(sp_samples)<4: continue
    clusters=sorted({rows[i]["cluster"] for i in idx}); ci={c:j for j,c in enumerate(clusters)}
    local={s:k for k,s in enumerate(sp_samples)}; M=np.zeros((len(sp_samples),len(clusters)))
    for i in idx:
        r=rows[i]; M[local[sample_label(r["#sample"])],ci[r["cluster"]]]=float(r["abundance"])
    totals=M.sum(1,keepdims=True); totals[totals==0]=1; M/=totals
    D=np.zeros((len(sp_samples),)*2)
    for i in range(len(sp_samples)):
        for j in range(i+1,len(sp_samples)):
            den=(M[i]+M[j]).sum(); D[i,j]=D[j,i]=abs(M[i]-M[j]).sum()/den if den else 0
    Ds[sp]=D
    r2,F,p=permanova(D,[subject_map[s] for s in sp_samples],999)
    perspecies.append({"species":sp,"strain_R2_subject":r2,"p_value":p,"n_samples_present":len(sp_samples),"n_clusters":len(clusters)})
perspecies.sort(key=lambda x:x["strain_R2_subject"],reverse=True)
write(OUT/"wms32_perspecies_subject.tsv",["species","strain_R2_subject","p_value","n_samples_present","n_clusters"],perspecies)
print(f"wrote {OUT/'wms32_temporal_ml.tsv'} and {OUT/'wms32_perspecies_subject.tsv'}")
