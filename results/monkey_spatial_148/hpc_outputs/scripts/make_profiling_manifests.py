from pathlib import Path
import pandas as pd,re
OUT=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/work/monkey_panel_20261023')
# 2bRAD manifest from current inventory
m=pd.read_csv('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/inventory/sample_manifest_2brad.tsv',sep='\t')
m=m[m.sample_category=='GI'][['sample_id','path']]
m.to_csv(OUT/'manifest_2brad.csv',index=False,columns=['sample_id','path'],header=['sample','reads1'])
# shotgun kneaddata GI manifest
sg=[]
for p in Path('/lustre1/g/aos_shihuang/data/monkey/shotgun/kneaddata_results').glob('*_L_MG_1.R1_kneaddata.fastq.gz'):
    sid=p.name.replace('.R1_kneaddata.fastq.gz','')
    if re.match(r'^(O1|O2|Y1|Y2)_(E|S|SI|LI|K)_',sid): sg.append((sid,str(p)))
pd.DataFrame(sg,columns=['sample','reads1']).sort_values('sample').to_csv(OUT/'manifest_shotgun.csv',index=False)
print('2bRAD',len(m),'shotgun',len(sg))
