from pathlib import Path
import re
import numpy as np,pandas as pd
from scipy.stats import spearmanr
W=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/work/monkey_panel_20261027')
W=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/work/monkey_panel_20261026')
W=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/work/monkey_panel_20261023')
F=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/analysis/functional_coupling_20261027')
OUT=F

def parse_sample(s):
    m=re.match(r'^(O1|O2|Y1|Y2)_(E|S|SI|LI|K)(?:_(\d+)(?:cm)?)?_(L|M)_MG_1$',s)
    if not m:return pd.Series({'monkey_id':None,'gi_region':None,'position_cm':np.nan,'compartment':None})
    monkey,region,pos,comp=m.groups()
    return pd.Series({'monkey_id':monkey,'gi_region':region,'position_cm':float(pos or np.nan),'compartment':{'L':'luminal','M':'mucosal'}[comp]})

# Shotgun strain profiles; abundance is within-species relative abundance.
d=pd.read_csv(W/'profile_shotgun/merged.tsv',sep='\t').rename(columns={'#sample':'sample_id'})
meta=pd.concat([d,d.sample_id.apply(parse_sample)],axis=1).groupby('sample_id').first()
# Within-species cluster composition matrix; zeros mean undetected in that sample.
S=d.pivot_table(index='sample_id',columns=['species','cluster'],values='abundance',aggfunc='max').fillna(0)
S.columns=[f'{a}|{b}' for a,b in S.columns]
# BA composition
BA=pd.read_csv(F/'bileacid_category_log2_matrix.tsv',sep='\t',index_col=0)
BA=2**BA
BA=BA.div(BA.sum(axis=1).replace(0,np.nan),axis=0)
BA=(BA+1e-6).apply(np.log2)
BA.index=BA.index.map(lambda x:x+'_L_MG_1')
# Function proportion of total unstratified HUMAnN pathway signal
allfun=pd.read_csv(F/'humann_pathabundance_cpm_unstratified.tsv',sep='\t',index_col=0)
sel=pd.read_csv(F/'humann_selected_function_matrix.tsv',sep='\t',index_col=0)
FUN=sel.div(allfun.loc[sel.index,allfun.columns].sum(axis=1).replace(0,np.nan),axis=0)*(1e6)
FUN=(FUN+1e-3).apply(np.log2)
# LI samples
li=meta[(meta.gi_region=='LI')&(meta.compartment=='luminal')].index
common_ba=li.intersection(S.index).intersection(BA.index)
common_fun=li.intersection(S.index).intersection(FUN.index)
print('LI samples',len(li),'strain+BA',len(common_ba),'strain+FUN',len(common_fun))

def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]
    out=np.empty(n); out[o]=np.minimum(q,1); return out

def residualize(df,cols,covs):
    Z=pd.get_dummies(df[covs[0]],drop_first=True).astype(float)
    Z['position']=df[covs[1]]; Z['intercept']=1.0
    Z=Z.loc[:,Z.std()!=0]
    out={}
    for c in cols:
        y=df[c].astype(float)
        beta=np.linalg.lstsq(Z.values,y.values,rcond=None)[0]
        out[c]=y.values-Z.values@beta
    return out

def assoc(matrixS,X,strain_cols,feature_cols,label,min_detect=8):
    covdf=meta.loc[matrixS.index,['monkey_id','position_cm']]
    s_in=pd.concat([matrixS,covdf],axis=1); x_in=pd.concat([X,covdf],axis=1)
    sr=residualize(s_in,strain_cols,covs=['monkey_id','position_cm'])
    fr=residualize(x_in,feature_cols,covs=['monkey_id','position_cm'])
    rows=[]
    for sc in strain_cols:
        if (matrixS[sc]>0).sum()<min_detect: continue
        for fc in feature_cols:
            r,p=spearmanr(sr[sc],fr[fc])
            rows.append({'strain_cluster':sc,'feature':fc,'n':len(matrixS),'n_detected':int((matrixS[sc]>0).sum()),'rho':r,'p':p})
    t=pd.DataFrame(rows); t['q']=bh(t.p); t.to_csv(OUT/f'{label}.tsv',sep='\t',index=False); return t

Sba=S.loc[common_ba]
det=(Sba>0).sum(); cols=det[det>=8].index.tolist()
print('clusters tested',len(cols))
ba_cols=list(BA.columns); fun_cols=list(FUN.columns)
tba=assoc(Sba.loc[common_ba,cols],BA.loc[common_ba,ba_cols],cols,ba_cols,'coupling_compositional_bileacid_li')
tfun=assoc(S.loc[common_fun,cols],FUN.loc[common_fun,fun_cols],cols,fun_cols,'coupling_compositional_functions_li')
print('\nTOP COMPOSITIONAL BILE ACID ASSOCIATIONS')
print(tba.sort_values('p').head(50).to_string(index=False))
print('\nTOP COMPOSITIONAL FUNCTION ASSOCIATIONS')
print(tfun.sort_values('p').head(80).to_string(index=False))
for label,t in [('BA',tba),('function',tfun)]:
    print(label,'tests',len(t),'q<0.10',sum(t.q<.10),'q<0.05',sum(t.q<.05),'p<0.001',sum(t.p<.001),'positive q<.05',sum((t.q<.05)&(t.rho>0)),'negative q<.05',sum((t.q<.05)&(t.rho<0)))
