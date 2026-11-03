import pandas as pd, numpy as np
from pathlib import Path
P=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/analysis/2brad_vendor')
d=pd.read_csv(P/'species_abundance_long.csv.gz',low_memory=False)
for c in ['Species','Genus','Family','sample_category','compartment','gi_region']:
    d[c]=d[c].astype('string').str.replace('-','_',regex=False)
gi=d[(d.sample_category=='GI') & (d.Species!='human')].copy()
# species-level stats
sp=gi.groupby(['Species','Genus','Family'],as_index=False).agg(
 mean_abundance=('abundance','mean'), median_abundance=('abundance','median'),
 max_abundance=('abundance','max'), prevalence=('abundance',lambda x:(x>0).mean()))
for comp in ['luminal','mucosal']:
    q=gi[gi.compartment==comp].groupby('Species').abundance.agg(**{f'mean_{comp}':'mean',f'prevalence_{comp}':lambda x:(x>0).mean()})
    sp=sp.merge(q,left_on='Species',right_index=True,how='left')
sp['log2FC_mucosal_vs_luminal_mean']=np.log2((sp.mean_mucosal.fillna(0)+1e-6)/(sp.mean_luminal.fillna(0)+1e-6))
sp['min_prevalence']=sp[['prevalence_luminal','prevalence_mucosal']].min(axis=1)
sp['max_prevalence']=sp[['prevalence_luminal','prevalence_mucosal']].max(axis=1)
# region specificity
reg=gi.groupby(['Species','gi_region']).abundance.mean().unstack(fill_value=0)
sp=sp.merge(reg.max(axis=1).rename('max_region_mean'),left_on='Species',right_index=True,how='left')
sp=sp.merge((reg.max(axis=1)/reg.sum(axis=1).replace(0,np.nan)).rename('region_specificity'),left_on='Species',right_index=True,how='left')
sp=sp.merge(reg.idxmax(axis=1).rename('peak_region'),left_on='Species',right_index=True,how='left')
# control detections
ctl=d[d.sample_category.isin(['body_control','negative_control','vendor_extra'])]
cv=ctl.groupby(['sample_category','Species']).abundance.agg(control_mean='mean',control_prevalence=lambda x:(x>0).mean(),control_max='max').reset_index()
for cat,sub in cv.groupby('sample_category'):
    cols={c:f'{c}_{cat}' for c in sub.columns if c not in ['sample_category','Species']}
    sp=sp.merge(sub[['Species']+list(cols)].rename(columns=cols),on='Species',how='left')
# Candidate tiers
sp['panel_tier']=np.select(
    [ (sp.max_prevalence>=.5)&(sp.max_abundance>=.1),
      (sp.max_prevalence>=.3)&(sp.max_abundance>=.01),
      (sp.max_prevalence>=.2)&(sp.max_abundance>=.001)],
    ['tier1_high','tier2_medium','tier3_low'],default='not_candidate')
sp['priority']=np.where(sp.Species.str.contains('Helicobacter|Lactobacillus|Bifidobacterium|Alloscardovia|Prevotella|Succinivibrio|Turicibacter|Escherichia|Bacteroides|Ligilactobacillus|Limosilactobacillus|Streptococcus|Enterococcus|Clostridium_',case=False),'yes','general')
sp.sort_values(['panel_tier','max_prevalence','mean_abundance'],ascending=[True,False,False]).to_csv(P/'panel_candidate_species.tsv',sep='\t',index=False)
# per-genus summary for panel design
gen=gi.groupby(["Genus","Family"],as_index=False).agg(n_species=("Species","nunique"),mean_abundance=("abundance","mean"),max_abundance=("abundance","max"))
pv=gi.groupby(["Genus","Species"]).abundance.agg(lambda x:(x>0).mean()).reset_index(name="prev")
gen=gen.merge(pv.groupby("Genus").prev.max().rename("max_species_prevalence"),on="Genus",how="left")
gen.to_csv(P/'panel_candidate_genera.tsv',sep='\t',index=False)
print('TIERS')
print(sp.panel_tier.value_counts().to_string())
print('\nTIER1/TIER2 PRIORITY SPECIES')
sel=sp[(sp.panel_tier.isin(['tier1_high','tier2_medium']))&(sp.priority=='yes')]
print(sel.sort_values(['max_prevalence','mean_abundance'],ascending=False).head(80).to_string(index=False))
print('\nTOP GENERA')
print(gen.head(50).to_string(index=False))
