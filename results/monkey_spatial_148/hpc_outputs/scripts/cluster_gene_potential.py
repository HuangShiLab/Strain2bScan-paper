from pathlib import Path
import re,pandas as pd
BASE=Path('/lustre1/g/aos_shihuang/data/monkey/shotgun/assembly')
DB=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/work/monkey_panel_20261023/dbs')
OUT=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/analysis/functional_coupling_20261027')
# Map genome key `sample_L_bin_N` to DB genome `sample.contigs.fa.metabat-bins_bin.N`
def mag_key(genome):
    return genome.replace('.contigs.fa.metabat-bins_bin.','_bin_')

# Build genome->key from all members files
rows=[]
for mp in DB.glob('*.members.tsv'):
    sp=mp.name.replace('.members.tsv','')
    m=pd.read_csv(mp,sep='\t',comment='#',names=['genome','cluster'])
    for _,r in m.iterrows():
        key=mag_key(r.genome)
        rows.append({'species':sp,'cluster':r.cluster,'genome':r.genome,'mag_key':key})
g=pd.DataFrame(rows)
# PF results
# Simpler robust parse
pf_rows=[]
for path,label in [(BASE/'prodigal/steroid_scan/PF01073_result.tsv','3Beta_HSD_PF01073'),(BASE/'prodigal/steroid_scan/PF02544_result.tsv','Steroid_dh_PF02544')]:
    for line in open(path):
        head,rest=line.split(':',1)
        mag=head.replace('.pfam.tbl','')
        f=rest.split()
        if len(f)>=6:
            pf_rows.append({'mag_key':mag,'function':label,'target':f[2],'evalue':float(f[4]),'score':float(f[5])})
pf=pd.DataFrame(pf_rows)
pfagg=pf.groupby(['mag_key','function']).size().unstack(fill_value=0).reset_index()
# 5BR hits
br=[]
for p in (BASE/'isoalloLCA/prodigal/5BR_hits').glob('*.5BR.di.tsv'):
    mag=p.name.replace('.5BR.di.tsv','')
    n=sum(1 for line in open(p) if line.strip())
    if n:
        # first line fields: mag, gene, subject, pident, length, ...
        fields=line.strip().split('\t') if False else next(open(p)).split('\t')
        br.append({'mag_key':mag,'5BR_hit_count':n,'5BR_subject':fields[2] if len(fields)>2 else '', '5BR_identity':float(fields[3]) if len(fields)>3 else np.nan, '5BR_aln':float(fields[4]) if len(fields)>4 else np.nan})
br=pd.DataFrame(br)
# merge per genome
g=g.merge(pfagg,on='mag_key',how='left').merge(br,on='mag_key',how='left')
for c in ['3Beta_HSD_PF01073','Steroid_dh_PF02544','5BR_hit_count']:
    if c not in g: g[c]=0
g[['3Beta_HSD_PF01073','Steroid_dh_PF02544','5BR_hit_count']]=g[['3Beta_HSD_PF01073','Steroid_dh_PF02544','5BR_hit_count']].fillna(0)
g['has_3Beta_HSD']=g['3Beta_HSD_PF01073']>0
g['has_Steroid_dh']=g['Steroid_dh_PF02544']>0
g['has_5BR_hit']=g['5BR_hit_count']>0
# cluster aggregate
agg=g.groupby(['species','cluster']).agg(
    n_members=('genome','nunique'),
    members_3Beta_HSD=('has_3Beta_HSD','sum'),
    members_Steroid_dh=('has_Steroid_dh','sum'),
    members_5BR_hit=('has_5BR_hit','sum'),
    genes_3Beta_HSD=('3Beta_HSD_PF01073','sum'),
    genes_Steroid_dh=('Steroid_dh_PF02544','sum'),
    genes_5BR=('5BR_hit_count','sum')).reset_index()
agg['cluster_id']=agg.species+'|'+agg.cluster
agg.to_csv(OUT/'strain_cluster_ba_gene_potential.tsv',sep='\t',index=False)
print('genomes',len(g),'clusters',len(agg),'members with 3Beta_HSD',g.has_3Beta_HSD.sum(),'Steroid_dh',g.has_Steroid_dh.sum(),'5BR',g.has_5BR_hit.sum())
