from pathlib import Path
import pandas as pd,re
BASE=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan')
OUT=BASE/'analysis/full_148_20261027'
WORK=BASE/'work/full_148_20261027'
for p in [OUT,WORK/'profile_2brad']:
    p.mkdir(parents=True,exist_ok=True)
m=pd.read_csv(BASE/'inventory/sample_manifest_2brad.tsv',sep='\t')
# 148-sample definition: all GI plus B_* body/control samples; exclude 4 vendor extra.
full=m[(m.sample_category=='GI')|((m.sample_category=='body_control')&m.sample_id.str.contains('_B_'))].copy()
full['sample_group']=full.sample_category.map({'GI':'GI_spatial','body_control':'body_control'})
full['body_site_code']=full.sample_id.str.extract(r'_B_([^_]+)_')
full['data_type']='2bRAD'
full=full.sort_values(['sample_group','monkey_id','gi_region','position_cm','compartment','sample_id'])
full[['sample_id','data_type','sample_group','sample_category','source','monkey_id','gi_region','position_cm','compartment','body_site_code','path']].to_csv(OUT/'full_148_sample_metadata.tsv',sep='\t',index=False)
pd.DataFrame({'sample':full.sample_id,'reads1':full.path}).to_csv(WORK/'manifest_full148_2brad.csv',index=False)
print('full cohort',len(full))
print(full.groupby('sample_group').size().to_string())
print(full.groupby(['sample_group','sample_category','compartment'],dropna=False).size().to_string())
