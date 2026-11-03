from pathlib import Path
import pandas as pd,re,sys
p=Path('/lustre1/g/aos_shihuang/data/monkey/shotgun/assembly/checkm/checkm_output.txt')
rows=[]; rx=re.compile(r"^\s*(\S+)\s+(.*?)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$")
for line in open(p):
    m=rx.match(line)
    if m:
        rows.append({"Bin Id":m.group(1),"Completeness":float(m.group(12)),"Contamination":float(m.group(13)),"Strain":float(m.group(14))})

d=pd.DataFrame(rows)
if d.empty: print('no checkm rows'); sys.exit(1)
d.to_csv('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/inventory/mag_checkm_qc.tsv',sep='\t',index=False)
print('rows',len(d),'unique bins',d['Bin Id'].nunique())
print(d[['Completeness','Contamination','Strain']].describe().to_string())
for comp,cont in [(90,5),(90,10),(80,10),(70,10)]:
 print('HQ90_5' if comp==90 and cont==5 else f'comp{comp}_cont{cont}',((d.Completeness>=comp)&(d.Contamination<=cont)).sum())
# join gtdb
g=pd.read_csv('/lustre1/g/aos_shihuang/data/monkey/shotgun/assembly/gtdbtk/classify/gtdbtk.bac120.summary.tsv',sep='\t')
a=pd.read_csv('/lustre1/g/aos_shihuang/data/monkey/shotgun/assembly/gtdbtk/classify/gtdbtk.ar53.summary.tsv',sep='\t') if Path('/lustre1/g/aos_shihuang/data/monkey/shotgun/assembly/gtdbtk/classify/gtdbtk.ar53.summary.tsv').exists() else pd.DataFrame()
g=pd.concat([g,a],ignore_index=True)
g['species']=g.classification.str.rsplit(';',n=1).str[-1].str.replace('^s__','',regex=True)
m=d.merge(g[['user_genome','species','closest_genome_ani']],left_on='Bin Id',right_on='user_genome',how='left')
m['HQ']= (m.Completeness>=90)&(m.Contamination<=5)
print('GTDB annotated',m.species.notna().sum(),'HQ annotated',m[m.HQ].species.notna().sum())
hq=m[m.HQ]
print('\nTOP HQ MAG species by genome count')
print(hq.groupby('species').agg(n_MAG=('Bin Id','nunique'),median_completeness=('Completeness','median'),median_contamination=('Contamination','median'),median_closest_ANI=('closest_genome_ani','median')).sort_values('n_MAG',ascending=False).head(80).to_string())
m.to_csv('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/inventory/mag_checkm_gtdb.tsv',sep='\t',index=False)
