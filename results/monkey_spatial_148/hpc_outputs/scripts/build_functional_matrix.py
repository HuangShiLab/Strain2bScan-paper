from pathlib import Path
import xml.etree.ElementTree as ET
import base64,zlib,re,os
import numpy as np,pandas as pd
from scipy.signal import find_peaks,savgol_filter
BA=Path('/lustre1/g/aos_shihuang/data/monkey/bile acid/Target')
OUT=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/analysis/functional_coupling_20261027')
OUT.mkdir(parents=True,exist_ok=True)

def arrays(ch):
    out={}
    for arr in ch.iter():
        if not arr.tag.endswith('binaryDataArray'): continue
        ps=[x.attrib.get('name') for x in arr.findall('{*}cvParam')]
        b=base64.b64decode(arr.find('{*}binary').text or '')
        try: b=zlib.decompress(b)
        except: pass
        typ=np.float64 if '64-bit float' in ps else np.float32
        a=np.frombuffer(b,dtype=typ).astype(float)
        if 'time array' in ps: out['time']=a
        elif 'intensity array' in ps: out['int']=a
    return out

def peak_metrics(t,y,lo=1.5,hi=16.0,halfwidth=0.30):
    m=(t>=lo)&(t<=hi); t=t[m]; y=y[m].astype(float)
    if len(y)<10: return None
    sm=savgol_filter(y,11,2) if len(y)>=11 else y
    im=int(np.argmax(sm))
    rt=float(t[im]); apex=float(sm[im])
    baseline=float(np.percentile(sm,10))
    if apex-baseline < max(100,0.05*np.max(sm)): return None
    left=max(0,np.searchsorted(t,rt-halfwidth)); right=min(len(t),np.searchsorted(t,rt+halfwidth,'right'))
    seg=sm[left:right]-baseline
    area=float(np.sum(seg*np.diff(t[left:right],append=t[right-1]+np.median(np.diff(t)))))
    snr=(apex-baseline)/(np.std(np.diff(sm))+1e-9)
    return rt,apex,baseline,area,max(0,snr)

def parse_file(p):
    root=ET.parse(p).getroot()
    rows=[]
    internal=[]
    for ch in root.iter():
        if not ch.tag.endswith('chromatogram'): continue
        cid=ch.attrib.get('id','')
        if 'SRM' not in cid: continue
        a=arrays(ch)
        if 'time' not in a or 'int' not in a: continue
        pm=peak_metrics(a['time'],a['int'])
        if pm is None: continue
        rt,apex,baseline,area,snr=pm
        q1=re.search(r'Q1=([0-9.]+)',cid); q3=re.search(r'Q3=([0-9.]+)',cid)
        exp=re.search(r'experiment=(\d+)',cid); trans=re.search(r'transition=(\d+)',cid)
        name=re.search(r'name=LM_baile acid_([^ ]+)',cid)
        row={'sample_id':re.sub(r'^031500-b1-','',p.stem),'chrom_id':cid,'Q1':float(q1.group(1)) if q1 else np.nan,
             'Q3':float(q3.group(1)) if q3 else np.nan,'experiment':int(exp.group(1)) if exp else np.nan,
             'transition':int(trans.group(1)) if trans else np.nan,'transition_name':name.group(1) if name else '',
             'rt':rt,'apex':apex,'baseline':baseline,'area':area,'snr':snr}
        rows.append(row)
        if 'IS0' in row['transition_name'] and row['experiment']==2:
            internal.append(row)
    return rows,internal

rows=[]; istds=[]
for p in sorted(BA.glob('*.mzML')):
    rr,ii=parse_file(p); rows+=rr; istds+=ii
print('chromatograms parsed',len(rows),'IS',len(istds))
allrows=pd.DataFrame(rows); allrows.to_csv(OUT/'srm_peak_metrics.tsv',sep='\t',index=False)
# BA samples only; exclude QC
allrows=allrows[~allrows.sample_id.str.contains('^QC')]
# Aggregate negative-mode target transition groups (experiment 1)
neg=allrows[allrows.experiment==1].copy()
neg['ba_group']=neg.transition_name.str.extract(r'^(\d+)_')
# internal standard total per sample from exp2
issum=pd.DataFrame(istds).groupby('sample_id').area.sum().rename('is_area') if istds else None
if issum is not None and len(issum): neg=neg.merge(issum,left_on='sample_id',right_index=True,how='left')
else: neg['is_area']=np.nan
neg['is_area']=neg.is_area.fillna(neg.groupby('sample_id').area.transform('sum'))
neg['area_isnorm']=neg.area/neg.is_area*1e6
# sum across transitions per named target group; if duplicate groups are isomers, retain by Q1 category+group
def category(q):
    if q<380: return 'LCA_like_Q1_375'
    if q<396: return 'dihydroxy_Q1_391'
    if q>=400 and q<412: return 'cholate_like_Q1_407'
    return 'conjugated_or_highmass'
neg['ba_category']=neg.Q1.apply(category)
# two matrices: target-group and chemical category
gm=neg.groupby(['sample_id','transition_name']).area_isnorm.sum().unstack(fill_value=0)
cm=neg.groupby(['sample_id','ba_category']).area_isnorm.sum().unstack(fill_value=0)
# log2 transform after 1
gm2=np.log2(gm+1); cm2=np.log2(cm+1)
gm2.to_csv(OUT/'bileacid_group_log2_matrix.tsv',sep='\t')
cm2.to_csv(OUT/'bileacid_category_log2_matrix.tsv',sep='\t')
neg.to_csv(OUT/'bileacid_transition_features.tsv',sep='\t',index=False)
print('BA groups',gm.shape,'categories',cm.shape,'samples',len(gm))
print(cm2.describe().T[['count','mean','50%','min','max']].to_string())

# HUMAnN aggregate
H=Path('/lustre1/g/aos_shihuang/data/monkey/shotgun/metaphlan/humann/July_run')
series={}
for p in H.glob('*pathabundance-cpm.tsv'):
    d=pd.read_csv(p,sep='\t')
    sid=re.sub(r'\.R1_kneaddata.*$','',p.name)
    series[sid]=pd.Series(pd.to_numeric(d.iloc[:,1],errors='coerce').fillna(0).values,index=d.iloc[:,0].astype(str))
X=pd.DataFrame(series).T.fillna(0)
# unstratified only
un=X[[i for i in X.columns if '|' not in i]]
un.to_csv(OUT/'humann_pathabundance_cpm_unstratified.tsv',sep='\t')
# selected pathways
sel_names=[
 'CENTFERM-PWY: pyruvate fermentation to butanoate',
 'PWY-5676: acetyl-CoA fermentation to butanoate II',
 'PWY-5677: succinate fermentation to butanoate',
 'P163-PWY: L-lysine fermentation to acetate and butanoate',
 'P108-PWY: pyruvate fermentation to propanoate I',
 'FERMENTATION-PWY: mixed acid fermentation',
 'ANAEROFRUCAT-PWY: homolactic fermentation',
 'P122-PWY: heterolactic fermentation',
 'P441-PWY: superpathway of N-acetylneuraminate degradation',
 'GLCMANNANAUT-PWY: superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation',
 'FUC-RHAMCAT-PWY: superpathway of fucose and rhamnose degradation',
 'FUCCAT-PWY: fucose degradation',
 'GLUCARDEG-PWY: D-glucarate degradation I',
 'GLUCARGALACTSUPER-PWY: superpathway of D-glucarate and D-galactarate degradation',
]
sel_names=[x for x in sel_names if x in un.columns]
un[sel_names].to_csv(OUT/'humann_selected_function_matrix.tsv',sep='\t')
print('HUMAnN samples/pathways',un.shape,'selected',len(sel_names))
