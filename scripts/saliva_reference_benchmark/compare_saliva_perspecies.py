#!/usr/bin/env python3
import sys
import numpy as np
from pathlib import Path

def read_batch(path):
    with open(path) as f:
        header=f.readline().lstrip('#').rstrip('\n').split('\t')
        rows=[dict(zip(header,line.rstrip('\n').split('\t'))) for line in f if line.strip()]
    return rows

def permanova(D,labels,perms=4999,seed=42):
    rng=np.random.default_rng(seed); labels=np.asarray(labels); N=len(labels)
    groups=np.unique(labels); a=len(groups)
    def ss(D2,lab):
        total=D2[np.triu_indices(N,1)].sum()/N; within=0.0
        for g in np.unique(lab):
            ix=np.where(lab==g)[0]; ng=len(ix)
            if ng>1: within += D2[np.ix_(ix,ix)][np.triu_indices(ng,1)].sum()/ng
        return total,within
    D2=D**2; total,within=ss(D2,labels); between=total-within
    F=(between/(a-1))/(within/(N-a)); R2=between/total; ge=1
    for _ in range(perms):
        lab=rng.permutation(labels); t,w=ss(D2,lab); b=t-w
        Fp=(b/(a-1))/(w/(N-a)); ge += Fp>=F
    return R2,F,ge/(perms+1)

def bray(M):
    n=len(M); D=np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            den=(M[i]+M[j]).sum()
            D[i,j]=D[j,i]=np.abs(M[i]-M[j]).sum()/den if den else 0.0
    return D

def bh(p):
    n=len(p); order=np.argsort(p); ranked=np.asarray(p)[order]
    q=ranked*n/(np.arange(n)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]
    out=np.empty_like(q); out[order]=np.minimum(q,1)
    return out

def analyze(name,rows,out):
    samples=sorted(set(r['sample'] for r in rows))
    labels=np.asarray([s.split('-')[1] for s in samples])
    si={s:i for i,s in enumerate(samples)}
    species=sorted(set(r['species'] for r in rows))
    records=[]
    for sp in species:
        srows=[r for r in rows if r['species']==sp]
        clusters=sorted(set(r['cluster'] for r in srows))
        fi={c:i for i,c in enumerate(clusters)}
        M=np.zeros((len(samples),len(clusters)))
        for r in srows:
            field='abundance' if 'abundance' in r else 'sample_fraction'
            M[si[r['sample']],fi[r['cluster']]] += float(r[field])
        total=M.sum(1,keepdims=True); total[total==0]=1; M=M/total
        pres=int((M.sum(1)>0).sum())
        if pres<10 or len(clusters)<2: continue
        R2,F,p=permanova(bray(M),labels)
        records.append((sp,pres,len(clusters),R2,F,p))
    q=bh([x[5] for x in records])
    for (sp,pres,ncl,R2,F,p),qq in zip(records,q):
        out.write(f"{name}\t{sp}\t{pres}\t{ncl}\t{R2:.6f}\t{F:.4f}\t{p:.6f}\t{qq:.6f}\t{int(qq<.05)}\n")

def main(paths):
    with open('saliva_reference_perspecies_subject.tsv','w') as o:
        o.write('panel\tspecies\tn_samples_present\tn_clusters\tstrain_R2_subject\tpseudo_F\tp_value\tBH_q\tsignificant\n')
        for name,path in paths:
            analyze(name,read_batch(path),o)
    print(Path('saliva_reference_perspecies_subject.tsv').resolve())
if __name__=='__main__': main([(x.split('=',1)[0],x.split('=',1)[1]) for x in sys.argv[1:]])
