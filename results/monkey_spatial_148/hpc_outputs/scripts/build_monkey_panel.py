from pathlib import Path
import pandas as pd, os, shutil
BASE=Path('/lustre1/g/aos_shihuang/data/monkey/shotgun/assembly')
WORK=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/work/monkey_panel_20261023')
GEN=WORK/'genome_dirs'; GEN.mkdir(parents=True,exist_ok=True)
SPECIES=[
 'Succinivibrio sp000431835','WG-1 sp900539665','UBA4372 sp900543815',
 'Helicobacter_B sp003364135','Oribacterium sp900768305',
 'Treponema_D succinifaciens','Treponema_D berlinense',
 'Ornithospirochaeta sp934247875','UBA1724 sp900548225',
 'Cryptobacteroides sp900546445','UBA636 sp900546285','Onthomorpha sp900545215',
 'Ruminococcus_E sp934476515','Ligilactobacillus animalis','CAG-873 sp001701165',
 'Cryptobacteroides sp000433355','Lachnospira eligens_A','CAG-349 sp003539515'
]
# Load GTDB
p=BASE/'gtdbtk/classify'
d=pd.concat([pd.read_csv(p/'gtdbtk.bac120.summary.tsv',sep='\t'),pd.read_csv(p/'gtdbtk.ar53.summary.tsv',sep='\t')],ignore_index=True)
def taxval(x,pre):
    x=str(x).split(';')
    return next((z[len(pre):] for z in reversed(x) if z.startswith(pre)),None)
d['species']=d.classification.apply(lambda x:taxval(x,'s__'))
d['genus']=d.classification.apply(lambda x:taxval(x,'g__'))
def mag_path(u):
    a,b=u.rsplit('_bin.',1)
    return BASE/'pass_filter'/a/f'bin.{b}.fa'
d['path']=d.user_genome.apply(mag_path); d['exists']=d.path.apply(os.path.exists)
rows=[]
for i,sp in enumerate(SPECIES):
    q=d[(d.species==sp)&d.exists].copy()
    safe=sp.replace(' ','_').replace('/','_').replace('-','_')
    outdir=GEN/safe
    if outdir.exists(): shutil.rmtree(outdir)
    outdir.mkdir(parents=True)
    links=[]
    for _,r in q.iterrows():
        dst=outdir/(r.user_genome+'.fa')
        if dst.exists() or dst.is_symlink(): dst.unlink()
        dst.symlink_to(r.path)
        links.append(str(dst))
    rows.append({'species':sp,'safe_species':safe,'n_mags':len(q),'n_samples':q.user_genome.str.split('.contigs',regex=False).str[0].nunique(),'genome_dir':str(outdir)})
pd.DataFrame(rows).to_csv(WORK/'selected_species_manifest.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False))
