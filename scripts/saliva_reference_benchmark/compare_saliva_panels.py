#!/usr/bin/env python3
import csv, sys
import numpy as np

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

def loo(D,labels):
    labels=np.asarray(labels); c=0
    for i in range(len(labels)):
        d=D[i].copy(); d[i]=np.inf; j=int(np.argmin(d)); c += labels[j]==labels[i]
    return c/len(labels)

def matrix(rows,level,field='sample_fraction'):
    samples=sorted(set(r['sample'] for r in rows))
    key=lambda r: f"{r['species']}|{r['cluster']}" if level=='strain' else r['species']
    feats=sorted(set(key(r) for r in rows))
    si={s:i for i,s in enumerate(samples)}; fi={f:i for i,f in enumerate(feats)}
    M=np.zeros((len(samples),len(feats)))
    for r in rows:
        M[si[r['sample']],fi[key(r)]] += float(r[field])
    total=M.sum(1,keepdims=True); total[total==0]=1
    return samples,M/total

def analyze(name,rows,out):
    samples,M_species=matrix(rows,'species'); _,M_strain=matrix(rows,'strain')
    labels=[s.split('-')[1] for s in samples]
    for level,M in [('species',M_species),('strain',M_strain)]:
        D=bray(M); R2,F,p=permanova(D,labels); acc=loo(D,labels)
        out.write(f"{name}\t{level}\t{M.shape[1]}\t{R2:.6f}\t{F:.3f}\t{p:.6f}\t{acc:.6f}\n")

def main(paths):
    with open('saliva_reference_panel_comparison.tsv','w') as o:
        o.write('panel\tlevel\tn_features\tR2\tpseudo_F\tp_value\tLOO_accuracy\n')
        for name,path in paths:
            analyze(name,read_batch(path),o)
    print(Path('saliva_reference_panel_comparison.tsv').resolve())
if __name__=='__main__': main([(x.split('=',1)[0],x.split('=',1)[1]) for x in sys.argv[1:]])
