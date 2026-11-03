from pathlib import Path
import subprocess, gzip, re
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd

BASE=Path('/lustre1/g/aos_shihuang/data/monkey/2bRAD')
OUT=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/analysis/2brad_provisional_clusters')
OUT.mkdir(parents=True,exist_ok=True)
CNM=Path('/lustre1/g/aos_shihuang/Strain2b/databases/CNM_0.01')
SPECIES=[
 'Helicobacter_B_macacae','Prevotella_sp000434975','Succinivibrio_sp000431835',
 'Prevotella_copri_A','Prevotella_hominis','Ligilactobacillus_animalis',
 'Limosilactobacillus_reuteri','Lactobacillus_johnsonii','Streptococcus_anginosus',
 'Escherichia_coli','Enterococcus_B_hirae','Alloscardovia_macacae',
 'Prevotella_pectinovora','Prevotella_intermedia'
]
# sample paths
files={}
for p in (BASE/'content/DZQD2025110503-b1/Clean_data').glob('*/*.fq.gz'):
    files[p.stem]=p
for p in (BASE/'mucosa/Cleandata').glob('*/*.fq.gz'):
    files[p.stem]=p
samples=[]
for monkey in ["O1","O2","Y1","Y2"]:
    for compartment in ["L","M"]:
        wanted=[("S","1"),("SI","1" if monkey!="O2" else "49"),("LI","1"),("LI","55")]
        for region,pos in wanted:
            sep="_" if region=="S" else "cm_"
            prefix=f"{monkey}_{region}_{pos}{sep}{compartment}_MG_1"
            matches=[]
            for source in ["content/DZQD2025110503-b1/Clean_data","mucosa/Cleandata"]:
                matches.extend((BASE/source).glob(f"{prefix}/*.fq.gz"))
            if matches:
                p=sorted(matches)[0]; sid=p.name.removesuffix(".fq.gz"); samples.append((sid,str(p)))
if len(samples)!=32:
    print("WARNING: pilot sample count is",len(samples),[s for s,_ in samples])

def parse_sample(sid):
    parts=sid.split('_')
    monkey=parts[0]; region=parts[1]
    raw=parts[2].removesuffix('cm')
    comp=parts[3] if parts[3] in ('L','M') else parts[4]
    return dict(sample_id=sid,monkey_id=monkey,gi_region=region,position_cm=int(raw),
                compartment={'L':'luminal','M':'mucosal'}[comp])

def process_dbs():
    COMP=str.maketrans("ACGT","TGCA")
    def canon(x): return min(x,x.translate(COMP))
    tag_to_species={}; tag_species_clusters={}; cluster_members={}
    for sp in SPECIES:
        p=CNM/f'{sp}.CNM.xls'
        df=pd.read_csv(p,sep='\t',index_col=0)
        df.index=df.index.str.slice(0,32)
        clusters=list(df.columns); cluster_members[sp]=clusters
        for tag,row in df.iterrows():
            tag=canon(str(tag).upper())
            present=tuple(c for c in clusters if row[c]==1)
            if tag in tag_to_species and tag_to_species[tag]!=sp:
                tag_to_species[tag]='__multi__'
            else:
                tag_to_species[tag]=sp
            tag_species_clusters.setdefault(tag,{}).setdefault(sp,set()).update(present)
    tag_clusters={tag:{sp:tuple(sorted(cs)) for sp,cs in ds.items()} for tag,ds in tag_species_clusters.items()}
    return tag_to_species,tag_clusters,cluster_members

def process_sample(arg):
    sid,path=arg
    meta=parse_sample(sid)
    proc=subprocess.Popen(['pigz','-dc',path],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
    total_reads=0; matched_reads=0
    species_reads=Counter(); species_tags=Counter(); tag_reads=Counter()
    f=proc.stdout
    while True:
        h=f.readline()
        if not h: break
        seq=f.readline().strip(); f.readline(); qual=f.readline()
        if not seq or not qual: break
        total_reads += 1
        seq=seq.decode('ascii').upper()
        COMP=str.maketrans('ACGT','TGCA')
        seq=min(seq,seq.translate(COMP))
        if seq in tag_to_species:
            sp=tag_to_species[seq]
            if sp!='__multi__':
                matched_reads += 1; species_reads[sp]+=1; tag_reads[(sp,seq)]+=1
    proc.stdout.close(); proc.wait()
    for (sp,tag),n in tag_reads.items():
        species_tags[sp]+=1
    rows=[]
    for sp in SPECIES:
        clusters=cluster_members[sp]
        assigned_reads=0; assigned_tags=0
        for cluster in clusters:
            cr=0; ct=0
            for (s2,tag),n in tag_reads.items():
                if s2==sp:
                    cs=tag_clusters[tag][sp]
                    if cs==(cluster,):
                        cr+=n; ct+=1
            assigned_reads+=cr; assigned_tags+=ct
            rows.append({**meta,'species':sp,'cluster':cluster,'total_reads':total_reads,
                         'species_db_reads':species_reads[sp],'species_db_tags':species_tags[sp],
                         'cluster_tag_reads':cr,'cluster_distinct_tags':ct,
                         'total_assigned_reads':None,'total_assigned_tags':None})
    for r in rows:
        r['total_assigned_reads']=assigned_reads; r['total_assigned_tags']=assigned_tags
    return rows

if __name__=='__main__':
    tag_to_species,tag_clusters,cluster_members=process_dbs()
    print('loaded',len(tag_to_species),'species-specific tags across',len(SPECIES),'species',flush=True)
    allrows=[]
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs={ex.submit(process_sample,a):a[0] for a in samples}
        for fut in as_completed(futs):
            sid=futs[fut]
            try: rows=fut.result(); allrows.extend(rows); print('done',sid,flush=True)
            except Exception as e: print('ERROR',sid,e,flush=True)
    df=pd.DataFrame(allrows)
    df.to_csv(OUT/'provisional_cluster_calls_pilot_raw.tsv',sep='\t',index=False)
    # Detection gate: at least two distinct cluster-unique tags observed
    det=df[(df.cluster_distinct_tags>=2)&(df.cluster_tag_reads>=2)].copy()
    # rough within-species relative marker-read abundance
    det['cluster_relative_abundance']=det.cluster_tag_reads/det.groupby(['sample_id','species']).cluster_tag_reads.transform('sum').replace(0,pd.NA)
    det.to_csv(OUT/'provisional_cluster_calls_pilot_detected.tsv',sep='\t',index=False)
    # summarize
    summ=det.groupby(['species','cluster','compartment','gi_region']).agg(
        n_samples=('sample_id','nunique'),mean_relative_abundance=('cluster_relative_abundance','mean'),
        median_relative_abundance=('cluster_relative_abundance','median'),
        mean_cluster_tag_reads=('cluster_tag_reads','mean')).reset_index()
    summ.to_csv(OUT/'provisional_cluster_spatial_summary.tsv',sep='\t',index=False)
    print('wrote',OUT,flush=True)
