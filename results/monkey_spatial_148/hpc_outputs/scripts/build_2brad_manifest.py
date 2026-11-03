from pathlib import Path
import pandas as pd,re
BASE=Path('/lustre1/g/aos_shihuang/data/monkey/2bRAD')
OUT=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/inventory')
files={}
for p in (BASE/'content/DZQD2025110503-b1/Clean_data').glob('*/*.fq.gz'):
    files[p.name.removesuffix(".fq.gz")]=('luminal_content',p)
for p in (BASE/'mucosa/Cleandata').glob('*/*.fq.gz'):
    files[p.name.removesuffix(".fq.gz")]=('mucosa',p)
rows=[]
for sid,(source,p) in files.items():
    m=re.match(r'^(O1|O2|Y1|Y2)_(E|S|SI|LI|K)(?:_(\d+)(?:cm)?)?_(L|M)_MG_1$',sid)
    if m: monkey,region,pos,comp=m.groups(); cat='GI'; compartment={'L':'luminal','M':'mucosal'}[comp]; pos=int(pos) if pos else None
    elif re.match(r'^(O1|O2|Y1|Y2)_B_',sid): monkey=region=pos=compartment=None; cat='body_control'
    elif '__MG_2' in sid: monkey=region=pos=compartment=None; cat='negative_control'
    else: monkey=region=pos=compartment=None; cat='vendor_extra'
    st=p.stat()
    rows.append(dict(sample_id=sid,data_type='2bRAD',source=source,sample_category=cat,monkey_id=monkey,gi_region=region,position_cm=pos,compartment=compartment,path=str(p),size_bytes=st.st_size))
pd.DataFrame(rows).sort_values(['sample_category','monkey_id','gi_region','compartment','sample_id']).to_csv(OUT/'sample_manifest_2brad.tsv',sep='\t',index=False)
# vendor luminal stat list
sl=BASE/'content/DZQD2025110503-b1/Clean_data/2bRAD_result/list/Abundance_Stat.list'
if sl.exists():
    d=pd.read_csv(sl,sep='\t',header=None,names=['sample_id','path'])
    d.to_csv(OUT/'vendor_file_manifest_luminal.tsv',sep='\t',index=False)
print('manifest samples',len(rows),'GI',sum(r['sample_category']=='GI' for r in rows))
print(pd.DataFrame(rows).groupby(['sample_category','source']).size().to_string())
