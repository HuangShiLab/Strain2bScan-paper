from pathlib import Path
import re
import numpy as np,pandas as pd
from scipy.stats import spearmanr
W=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/work/monkey_panel_20261023')
F=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/analysis/functional_coupling_20261027')
OUT=F; OUT.mkdir(parents=True,exist_ok=True)

def parse_sample(s):
    m=re.match(r'^(O1|O2|Y1|Y2)_(E|S|SI|LI|K)(?:_(\d+)(?:cm)?)?_(L|M)_MG_1$',s)
    if not m:return pd.Series({'monkey_id':None,'gi_region':None,'position_cm':np.nan,'compartment':None})
    monkey,region,pos,comp=m.groups()
    return pd.Series({'monkey_id':monkey,'gi_region':region,'position_cm':float(pos or np.nan),'compartment':{'L':'luminal','M':'mucosal'}[comp]})

# strain shotgun matrix
d=pd.read_csv(W/'profile_shotgun/merged.tsv',sep='\t').rename(columns={'#sample':'sample_id'})
meta=pd.concat([d,d.sample_id.apply(parse_sample)],axis=1)
# cluster abundance matrix (global cross-panel abundance)
S=d.pivot_table(index='sample_id',columns=['species','cluster'],values='global_abundance',aggfunc='max').fillna(0)
S.columns=[f'{a}|{b}' for a,b in S.columns]
meta=meta.groupby('sample_id').first()
# bile matrices
BAcat=pd.read_csv(F/'bileacid_category_log2_matrix.tsv',sep='\t',index_col=0)
BAgrp=pd.read_csv(F/'bileacid_group_log2_matrix.tsv',sep='\t',index_col=0)
FUN=pd.read_csv(F/'humann_selected_function_matrix.tsv',sep='\t',index_col=0)
# normalize BA sample O1_LI_10cm -> O1_LI_10cm_L_MG_1
def ba_to_strain(x): return x+'_L_MG_1'
BAcat.index=BAcat.index.map(ba_to_strain); BAgrp.index=BAgrp.index.map(ba_to_strain)
# LI only
li=meta[(meta.gi_region=='LI')&(meta.compartment=='luminal')].index
common_S_BA=li.intersection(S.index).intersection(BAcat.index)
common_S_F=li.intersection(S.index).intersection(FUN.index)
print('LI samples in strain',len(li),'strain+BA',len(common_S_BA),'strain+FUN',len(common_S_F))

def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); ranked=p[o]
    q=np.minimum.accumulate((ranked*n/np.arange(1,n+1))[::-1])[::-1]
    out=np.empty(n); out[o]=np.minimum(q,1); return out

def resid_vars(df,cols,covs):
    # residualize columns on intercept + categorical monkey + continuous position
    Z=pd.get_dummies(df[covs[0]],drop_first=True).astype(float)
    Z['intercept']=1.0; Z['position']=df[covs[1]]
    Z=Z.loc[:,Z.std()!=0]
    out={}
    for c in cols:
        y=df[c].astype(float)
        beta=np.linalg.lstsq(Z.values,y.values,rcond=None)[0]
        out[c]=y.values-Z.values@beta
    return out

def associate(strain_cols, feature_cols, feature_df, label, min_detect=8):
    subs=meta.loc[common_S_F if label.startswith('function') else common_S_BA]
    matrix_S=S.loc[subs.index,strain_cols]
    X=feature_df.loc[subs.index,feature_cols]
    covs=['monkey_id','position_cm']
    covdf=meta.loc[subs.index,covs]
    s_in=pd.concat([matrix_S,covdf],axis=1)
    f_in=pd.concat([X,covdf],axis=1)
    sr=resid_vars(s_in,strain_cols,covs)
    fr=resid_vars(f_in,feature_cols,covs)
    rows=[]
    for sc in strain_cols:
        if (matrix_S[sc]>0).sum()<min_detect: continue
        for fc in feature_cols:
            r,p=spearmanr(sr[sc],fr[fc])
            rows.append({'strain_cluster':sc,'feature':fc,'n':len(subs),'n_detected':int((matrix_S[sc]>0).sum()),'rho':r,'p':p})
    t=pd.DataFrame(rows)
    t['q']=bh(t.p)
    t.to_csv(OUT/f'coupling_{label}.tsv',sep='\t',index=False)
    return t

# selected clusters: detected in >=8 LI samples
Sli=S.loc[common_S_BA]
det=(Sli>0).sum()
strain_cols=det[det>=8].index.tolist()
print('clusters tested',len(strain_cols),'of',Sli.shape[1])
# category BA tests
ba_cols=list(BAcat.columns)
tba=associate(strain_cols,ba_cols,BAcat,'bileacid_categories_li')
print('\nTOP BILE-ACID CATEGORY ASSOCIATIONS')
print(tba.sort_values('p').head(40).to_string(index=False))
# group BA tests
tgr=associate(strain_cols,list(BAgrp.columns),BAgrp,'bileacid_groups_li')
print('\nTOP BILE-ACID GROUP ASSOCIATIONS')
print(tgr.sort_values('p').head(40).to_string(index=False))
# functional tests
fun_cols=list(FUN.columns)
tf=associate(strain_cols,fun_cols,FUN,'functions_li')
print('\nTOP FUNCTION ASSOCIATIONS')
print(tf.sort_values('p').head(60).to_string(index=False))
# summarize significant
for label,t in [('bileacid_categories_li',tba),('bileacid_groups_li',tgr),('functions_li',tf)]:
    print(label,'tests',len(t),'q<0.10',sum(t.q<.10),'q<0.05',sum(t.q<.05),'p<0.001',sum(t.p<.001))
