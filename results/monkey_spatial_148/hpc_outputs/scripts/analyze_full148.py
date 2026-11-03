from pathlib import Path
import re, warnings
import numpy as np,pandas as pd
from scipy.stats import wilcoxon,spearmanr
BASE=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan')
WORK=BASE/'work/full_148_20261027'
F=BASE/'analysis/functional_coupling_20261027'
OUT=BASE/'analysis/full_148_20261027'; OUT.mkdir(parents=True,exist_ok=True)
warnings.filterwarnings('ignore')

meta=pd.read_csv(OUT/'full_148_sample_metadata.tsv',sep='\t').set_index('sample_id')
d=pd.read_csv(WORK/'profile_2brad/merged.tsv',sep='\t').rename(columns={'#sample':'sample_id'})
d=d.merge(meta.reset_index(),on='sample_id',how='left',suffixes=('','_meta'))
d['cluster_id']=d.species+'|'+d.cluster
d.to_csv(OUT/'full148_2brad_strain_calls_long.tsv',sep='\t',index=False)

# sample-level summary
def shannon(x):
    x=x[x>0].values.astype(float)
    if len(x)==0:return 0
    p=x/x.sum(); return float(-(p*np.log2(p)).sum())
rows=[]
for sid,g in d.groupby('sample_id'):
    rows.append({'sample_id':sid,'n_strain_calls':g.cluster_id.nunique(),'n_species':g.species.nunique(),
                 'summed_global_abundance':g.global_abundance.sum(),'strain_shannon':shannon(g.global_abundance),
                 'median_support':g.support.median(),'median_coverage':g.coverage.median(),
                 'max_depth':g.depth.max()})
samp=pd.DataFrame(rows)
samp=meta.join(samp.set_index('sample_id'))
for c in ['n_strain_calls','n_species','summed_global_abundance','strain_shannon','median_support','median_coverage','max_depth']:
    samp[c]=samp[c].fillna(0).astype(float)
samp['has_2brad_strain_call']=samp.n_strain_calls>0
samp.to_csv(OUT/'full148_sample_strain_summary.tsv',sep='\t')

# integrated availability
shot=pd.read_csv(BASE/'work/monkey_panel_20261023/profile_shotgun/merged.tsv',sep='\t').rename(columns={'#sample':'sample_id'})
if 'cluster_id' not in shot: shot['cluster_id']=shot.species+'|'+shot.cluster
shotgrp=shot.groupby('sample_id').agg(n_shotgun_strain_calls=('cluster_id','nunique'),n_shotgun_species=('species','nunique'))
shot_manifest=pd.read_csv(BASE/'work/monkey_panel_20261023/manifest_shotgun.csv')
# join zero counts for sequenced samples without calls
shotgrp=shot_manifest.set_index('sample').join(shotgrp,how='left').fillna({'n_shotgun_strain_calls':0,'n_shotgun_species':0})
FUN=pd.read_csv(F/'humann_selected_function_matrix.tsv',sep='\t',index_col=0)
BA=pd.read_csv(F/'bileacid_category_log2_matrix.tsv',sep='\t',index_col=0)
BA.index=BA.index.map(lambda x:x+'_L_MG_1')
samp['has_shotgun']=samp.index.isin(shotgrp.index)
samp['has_humann_pathway']=samp.index.isin(FUN.index)
samp['has_bileacid_srml']=samp.index.isin(BA.index)
samp['analysis_tier']=np.select(
    [samp.sample_group=='body_control',
     (samp.sample_group=='GI_spatial') & samp.has_bileacid_srml,
     (samp.sample_group=='GI_spatial') & samp.has_shotgun & samp.has_humann_pathway,
     samp.sample_group=='GI_spatial'],
    ['body_control','GI_multimodal_with_BA','GI_paired_2bRAD_shotgun_function','GI_2bRAD_only'],default='other')
samp.to_csv(OUT/'integrated_sample_master_148.tsv',sep='\t')

# summaries by sample group/compartment
gi=meta.index[meta.sample_group=='GI_spatial']
dc=d[d.sample_id.isin(gi)].copy()
# species prevalence by sample set
def prev_table(dd,label,total=None):
    nd=total if total is not None else dd.sample_id.nunique()
    return dd.groupby('species').agg(**{f'n_{label}':('sample_id','nunique'),f'prevalence_{label}':('sample_id',lambda x:x.nunique()/nd),
                                        f'mean_global_abundance_{label}':('global_abundance','mean'),f'median_global_abundance_{label}':('global_abundance','median'),
                                        f'max_global_abundance_{label}':('global_abundance','max')})
pall=prev_table(d[d.sample_group=='GI_spatial'],'GI',128)
plum=prev_table(d[(d.sample_group=='GI_spatial')&(d.compartment=='luminal')],'luminal',66)
pmuc=prev_table(d[(d.sample_group=='GI_spatial')&(d.compartment=='mucosal')],'mucosal',62)
pctl=prev_table(d[d.sample_group=='body_control'],'body_control',20)
sp=pall.join(plum,how='outer').join(pmuc,how='outer').join(pctl,how='outer').fillna(0)
sp['prevalence_ratio_GI_over_body']=(sp.prevalence_GI+1e-6)/(sp.prevalence_body_control+1e-6)
sp['potential_background_flag']=(sp.prevalence_body_control>=0.20)&(sp.prevalence_ratio_GI_over_body<2)
sp.sort_values('prevalence_GI',ascending=False).to_csv(OUT/'full148_species_GI_control_summary.tsv',sep='\t')

# cluster GI/control summary
def cprev(dd,label,total=None):
    nd=total if total is not None else dd.sample_id.nunique()
    return dd.groupby(['species','cluster']).agg(**{f'n_{label}':('sample_id','nunique'),f'prevalence_{label}':('sample_id',lambda x:x.nunique()/nd),
                                                    f'mean_global_abundance_{label}':('global_abundance','mean'),f'median_global_abundance_{label}':('global_abundance','median')})
c_gi=cprev(dc,'GI',128); c_lum=cprev(d[(d.sample_group=='GI_spatial')&(d.compartment=='luminal')],'luminal',66)
c_muc=cprev(d[(d.sample_group=='GI_spatial')&(d.compartment=='mucosal')],'mucosal',62); c_ctl=cprev(d[d.sample_group=='body_control'],'body_control',20)
cl=c_gi.join(c_lum,how='outer').join(c_muc,how='outer').join(c_ctl,how='outer').fillna(0).reset_index()
cl['cluster_id']=cl.species+'|'+cl.cluster
cl['prevalence_ratio_GI_over_body']=(cl.prevalence_GI+1e-6)/(cl.prevalence_body_control+1e-6)
cl['potential_background_flag']=(cl.prevalence_body_control>=0.20)&(cl.prevalence_ratio_GI_over_body<2)
# region specificity across GI
reg=dc.groupby(['species','cluster','gi_region']).agg(n=('sample_id','nunique'),mean_abund=('global_abundance','mean')).reset_index()
denom=dc.groupby('gi_region').sample_id.nunique().rename('denom').reset_index()
reg=reg.merge(denom,on='gi_region'); reg['prevalence']=reg.n/reg.denom
piv=reg.pivot(index=['species','cluster'],columns='gi_region',values='prevalence').fillna(0)
meanpiv=reg.pivot(index=['species','cluster'],columns='gi_region',values='mean_abund').fillna(0)
cl=cl.set_index(['species','cluster'])
cl['max_region_prevalence']=piv.max(axis=1); cl['peak_region']=piv.idxmax(axis=1)
cl['region_specificity']=piv.max(axis=1)/piv.sum(axis=1).replace(0,np.nan)
for c in piv.columns:
    cl[f'prevalence_{c}']=piv[c]
    cl[f'mean_global_abundance_{c}']=meanpiv[c]
cl=cl.reset_index()
cl.sort_values(['max_region_prevalence','region_specificity'],ascending=False).to_csv(OUT/'full148_2brad_cluster_GI_control_region_summary.tsv',sep='\t',index=False)


def _pairs_index(pairs,col):
    return pairs[col].map(lambda x:S.index.get_loc(x)).values

# mucosal/luminal paired tests on matched monkey-region-position
gi_meta=meta.loc[gi]
pair_records=[]
for (monkey,region,pos),g in gi_meta.groupby(['monkey_id','gi_region','position_cm']):
    gl=g[g.compartment=='luminal']; gm=g[g.compartment=='mucosal']
    if len(gl)==1 and len(gm)==1:
        pair_records.append({'monkey_id':monkey,'gi_region':region,'position_cm':pos,
                             'sample_luminal':gl.index[0],'sample_mucosal':gm.index[0]})
pairs=pd.DataFrame(pair_records)
S=dc.pivot_table(index='sample_id',columns='cluster_id',values='global_abundance',aggfunc='max').reindex(gi_meta.index).fillna(0)
S=meta.loc[S.index,['monkey_id','gi_region','position_cm','compartment']].join(S)
tests=[]
for clu in [c for c in S.columns if c.startswith('species_hash') or True]:
    if clu in ['monkey_id','gi_region','position_cm','compartment']: continue
    lv=S.loc[pairs.sample_luminal.values,clu].values; mv=S.loc[pairs.sample_mucosal.values,clu].values
    if (lv>0).sum()>=4 and (mv>0).sum()>=4:
        try:
            st,p=wilcoxon(mv,lv,zero_method='wilcox')
            parts=clu.split('|'); species='|'.join(parts[:-1]); cluster=parts[-1]
            tests.append({'species':species,'cluster':cluster,'n_pairs':len(pairs),
                          'n_mucosal_positive':int((mv>0).sum()),'n_luminal_positive':int((lv>0).sum()),
                          'median_mucosal':float(np.median(mv)),'median_luminal':float(np.median(lv)),
                          'median_log2FC_mucosal_luminal':np.log2((np.median(mv)+1e-8)/(np.median(lv)+1e-8)),
                          'W':st,'p_two_sided':p})
        except: pass
tt=pd.DataFrame(tests)
if len(tt):
    tt['p_bonferroni']=np.minimum(tt.p_two_sided*len(tt),1)
    tt=tt.sort_values('p_two_sided')
tt.to_csv(OUT/'full148_mucosal_luminal_paired_tests.tsv',sep='\t',index=False)


# report
with open(OUT/'FULL148_ANALYSIS.md','w') as out:
    out.write('# Full 148-sample monkey 2bRAD strain analysis\n\n')
    out.write('## Cohort definition\n\n- 148 samples = 128 GI spatial samples + 20 body/control samples.\n')
    out.write('- 128 GI samples = 66 luminal + 62 mucosal.\n')
    out.write('- Four extra vendor throat/trachea files were excluded from the 148-sample set.\n\n')
    out.write('## Sample-level profiling\n\n```text\n')
    out.write(samp.groupby('sample_group')[['n_strain_calls','n_species','strain_shannon']].agg(['count','mean','median','min','max']).round(3).to_string()+'\n```\n\n')
    out.write('## GI vs body-control species prevalence\n\n```text\n')
    cols=['prevalence_GI','prevalence_luminal','prevalence_mucosal','prevalence_body_control','prevalence_ratio_GI_over_body','potential_background_flag']
    out.write(sp.sort_values('prevalence_GI',ascending=False)[cols].round(3).to_string()+'\n```\n\n')
    out.write('## LI region and compartment signals\n\n')
    nli=cl[(cl.peak_region=='LI')&(cl.region_specificity>=0.8)&(cl.max_region_prevalence>=0.3)]
    out.write(f'- LI-enriched clusters (specificity >= 0.8, peak prevalence >= 0.3): **{len(nli)}**.\n')
    if len(tt):
        out.write(f'- Matched mucosal/luminal cluster tests: **{len(tt)}**; Bonferroni q<0.05: **{(tt.p_bonferroni<0.05).sum()}**.\n')
        out.write(f'- Nominal p<0.05: **{(tt.p_two_sided<0.05).sum()}**.\n')
    out.write('\nSee table outputs for full results.\n')
print('Wrote',OUT)
print('\nSAMPLE SUMMARY')
print(samp.groupby('sample_group')[['n_strain_calls','n_species','strain_shannon']].agg(['count','mean','median','min','max']).round(3).to_string())
print('\nAVAILABILITY')
print(samp.groupby(['sample_group','analysis_tier']).size().to_string())
print('\nTOP GI SPECIES')
print(sp.sort_values('prevalence_GI',ascending=False)[cols].head(30).round(3).to_string())
print('\nREGION/COMPARTMENT')
print('LI-enriched clusters',len(nli),'paired tests',len(tt),'bonf sig',0 if len(tt)==0 else (tt.p_bonferroni<0.05).sum())
