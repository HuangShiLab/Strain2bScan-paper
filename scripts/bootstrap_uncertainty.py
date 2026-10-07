#!/usr/bin/env python3
"""Subject-aware bootstrap uncertainty for saliva host-ID and paired benchmark differences."""
import csv
from collections import defaultdict
import sys
from pathlib import Path
import numpy as np

PAPER = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PAPER / "scripts"))
from saliva_table import read_long, COMMUNITY_ABUNDANCE  # noqa: E402

RNG = np.random.default_rng(20261006)

def wilson(k, n, z=1.959963985):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z*z/n
    centre = (p + z*z/(2*n)) / den
    half = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return centre-half, centre+half

def saliva_bootstrap(nrep=20000):
    rows = read_long(PAPER / "results/saliva_strain_long.tsv", require=(COMMUNITY_ABUNDANCE,))
    samples = sorted({r["sample"] for r in rows})
    si = {s:i for i,s in enumerate(samples)}
    subject = np.array([s.split("-")[1] for s in samples])
    timepoint = np.array([s.split("-")[0] for s in samples])
    def matrix(level):
        key = ((lambda r: f"{r['species']}|{r['cluster']}") if level == "strain" else (lambda r: r["species"]))
        feats = sorted({key(r) for r in rows}); fi={f:j for j,f in enumerate(feats)}
        M=np.zeros((len(samples),len(feats)))
        for r in rows: M[si[r["sample"]],fi[key(r)]] += float(r[COMMUNITY_ABUNDANCE])
        total=M.sum(1,keepdims=True); total[total==0]=1
        return M/total
    def bc(M):
        n=len(M); D=np.zeros((n,n))
        for i in range(n):
            for j in range(i+1,n):
                den=(M[i]+M[j]).sum(); D[i,j]=D[j,i]=abs(M[i]-M[j]).sum()/den if den else 0
        return D
    def loto(M, subject, timepoint):
        D=bc(M); tps=sorted(set(timepoint)); correct=0
        for held in tps:
            test=np.where(timepoint==held)[0]; train=np.where(timepoint!=held)[0]
            labels = subject if subject is not None else subject
            for i in test:
                nn=train[np.argmin(D[i,train])]
                correct += subject[nn] == subject[i]
        return correct/len(samples)
    Ms=matrix("strain"); Mp=matrix("species")
    point=(loto(Ms,subject,timepoint),loto(Mp,subject,timepoint))
    # Subject-level cluster bootstrap: resample 8 subjects, retaining all 4 timepoints.
    subjects=np.unique(subject); stats=[]
    for _ in range(nrep):
        chosen=RNG.choice(subjects,len(subjects),replace=True)
        idx=np.concatenate([np.where(subject==s)[0] for s in chosen])
        sub=subject[idx]
        tp_sub=timepoint[idx]
        # Relabel duplicated subjects to preserve the classification identity after resampling.
        relabel={old:f"{old}_{rep}" for rep,old in enumerate(chosen)}
        sub=np.array([relabel[x] for x in sub])
        stats.append((loto(Ms[idx],sub,tp_sub),loto(Mp[idx],sub,tp_sub)))
    arr=np.array(stats); out={"n_subjects":len(subjects),"n_libraries":len(samples),"n_bootstrap":nrep}
    for j,level in enumerate(("strain","species")):
        lo,hi=np.quantile(arr[:,j],[.025,.975]); out[f"{level}_point"]=point[j]; out[f"{level}_ci"]=(lo,hi)
    out["difference_point"]=point[0]-point[1]; out["difference_ci"]=(np.quantile(arr[:,0]-arr[:,1],.025),np.quantile(arr[:,0]-arr[:,1],.975))
    return out

def benchmark_bootstrap(nrep=20000):
    def load(path):
        with open(path) as fh: return list(csv.DictReader(fh,delimiter="\t"))
    a=load(PAPER/"figure_raw_data/sim_headtohead/strain2bscan_single_native_persample.tsv")
    b=load(PAPER/"work/strainscan_headtohead_rerun/scored_runs.tsv")
    bmap={}
    for r in b:
        if r["kind"] != "single" or r.get("status", "completed") != "completed":
            continue
        final=PAPER/"work/strainscan_headtohead_rerun/results/single"/r["sample"]/"final_report.txt"
        finite_metrics=all(r[k] not in ("", "nan", "NaN") for k in ("recall", "f1"))
        if final.exists() and final.stat().st_size and finite_metrics:
            bmap[r["sample"]]=r
    pairs=[]
    for ra in a:
        rb=bmap.get(ra["sample"])
        if rb:
            pairs.append((ra["species"],float(ra["recall"]),float(rb["recall"]),float(ra["f1"]),float(rb["f1"])))
    by=defaultdict(list)
    for x in pairs: by[x[0]].append(x)
    species=sorted(by); point=(np.mean([x[1]-x[2] for x in pairs]),np.mean([x[3]-x[4] for x in pairs]))
    stats=[]
    for _ in range(nrep):
        chosen=RNG.choice(species,len(species),replace=True)
        sel=[x for s in chosen for x in by[s]]
        stats.append((np.mean([x[1]-x[2] for x in sel]),np.mean([x[3]-x[4] for x in sel])))
    arr=np.array(stats)
    return {"n_paired":len(pairs),"recall_difference_point":point[0],
            "recall_difference_ci":tuple(np.quantile(arr[:,0],[.025,.975])),
            "f1_difference_point":point[1],"f1_difference_ci":tuple(np.quantile(arr[:,1],[.025,.975]))}

if __name__ == "__main__":
    saliva=saliva_bootstrap(); bench=benchmark_bootstrap()
    with open(PAPER/"results/uncertainty_summary.tsv","w",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(["analysis","metric","point","ci_low","ci_high"])
        w.writerow(["saliva_leave_one_timepoint_out","strain_accuracy",saliva["strain_point"],*saliva["strain_ci"]])
        w.writerow(["saliva_leave_one_timepoint_out","species_accuracy",saliva["species_point"],*saliva["species_ci"]])
        w.writerow(["saliva_leave_one_timepoint_out","strain_minus_species_accuracy",saliva["difference_point"],*saliva["difference_ci"]])
        w.writerow(["benchmark_paired_single_species","recall_difference",bench["recall_difference_point"],*bench["recall_difference_ci"]])
        w.writerow(["benchmark_paired_single_species","f1_difference",bench["f1_difference_point"],*bench["f1_difference_ci"]])
    print("saliva",saliva); print("benchmark",bench)
