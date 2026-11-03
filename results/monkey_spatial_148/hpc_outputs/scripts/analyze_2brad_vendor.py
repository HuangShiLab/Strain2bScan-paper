from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

BASE=Path("/lustre1/g/aos_shihuang/data/monkey/2bRAD")
OUT=Path("/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/analysis/2brad_vendor")
OUT.mkdir(parents=True,exist_ok=True)
paths={
 "luminal":BASE/"content/DZQD2025110503-b1/Clean_data/2bRAD_result/quantitative/Abundance_Stat.filtered.xls",
 "mucosa":BASE/"mucosa/Abundance.filtered.anno.xls",
}
meta_cols={"#Kingdom","Phylum","Class","Order","Family","Genus","Species","ChineseName","Isolated"}
frames=[]
for dtype,p in paths.items():
    df=pd.read_csv(p,sep="\t",encoding="gb18030")
    tax_cols=[c for c in df.columns if c in ["#Kingdom","Phylum","Class","Order","Family","Genus","Species"]]
    samples=[c for c in df.columns if c not in meta_cols]
    long=df.melt(id_vars=tax_cols,value_vars=samples,var_name="sample_id",value_name="abundance")
    long["data_type"]="2bRAD"; long["vendor_compartment"]=dtype
    frames.append(long)
long=pd.concat(frames,ignore_index=True)
long["abundance"]=pd.to_numeric(long.abundance,errors="coerce").fillna(0)

def parse_sample(s):
    m=re.match(r"^(O1|O2|Y1|Y2)_(E|S|SI|LI|K)(?:_(\d+)(?:cm)?)?_(L|M)_MG",s)
    if m:
        monkey,region,pos,compartment=m.groups()
        return pd.Series({"monkey_id":monkey,"gi_region":region,"position_cm":float(pos) if pos else np.nan,
                          "compartment":{"L":"luminal","M":"mucosal"}[compartment],"sample_category":"GI"})
    if re.match(r"^(O1|O2|Y1|Y2)_B_",s): category="body_control"
    elif "__MG_2" in s: category="negative_control"
    elif any(x in s for x in ["Lip","shoulder","Leg"]): category="body_control"
    elif s.startswith(("O1_Tg","Y1_Th","Y2_Tg","S3779425")): category="vendor_extra"
    else: category="other"
    return pd.Series({"monkey_id":None,"gi_region":None,"position_cm":np.nan,"compartment":None,"sample_category":category})

long=pd.concat([long,long.sample_id.apply(parse_sample)],axis=1)
long.to_csv(OUT/"species_abundance_long.csv.gz",index=False,compression="gzip")
gi=long[long.sample_category=="GI"].copy()
tax=gi[["Species","#Kingdom","Phylum","Class","Order","Family","Genus"]].drop_duplicates().set_index("Species")
tax=tax[~tax.index.duplicated(keep="first")]
wide=gi.pivot_table(index="sample_id",columns="Species",values="abundance",aggfunc="sum").fillna(0)
wide=wide[tax.reindex(wide.columns).index]
wide.to_csv(OUT/"species_abundance_wide_gi.tsv",sep="\t")
sample_meta=long[long.sample_category=="GI"].groupby("sample_id").first()[["monkey_id","gi_region","position_cm","compartment","vendor_compartment"]].sort_index()
sample_meta.to_csv(OUT/"sample_metadata_gi.tsv",sep="\t")
alpha=pd.DataFrame(index=wide.index)
alpha["observed_species"]=(wide>0).sum(axis=1)
p=wide.div(wide.sum(axis=1).replace(0,np.nan),axis=0).fillna(0)
alpha["shannon"]=-(p.replace(0,np.nan)*np.log2(p.replace(0,np.nan))).sum(axis=1).fillna(0)
alpha["simpson"]=1-(p**2).sum(axis=1)
alpha=pd.concat([sample_meta,alpha],axis=1)
alpha.to_csv(OUT/"alpha_diversity_gi.tsv",sep="\t")
counts=sample_meta.groupby(["monkey_id","gi_region","compartment"]).size().reset_index(name="n").sort_values(["monkey_id","gi_region","compartment"])
counts.to_csv(OUT/"group_sample_counts.tsv",sep="\t",index=False)
sp=gi.groupby(["Species","compartment"]).abundance.agg(mean_abundance="mean",median_abundance="median",prevalence=lambda x:(x>0).mean(),max_abundance="max").reset_index()
sp_p=sp.pivot(index="Species",columns="compartment",values=["mean_abundance","median_abundance","prevalence","max_abundance"])
sp_p.columns=[f"{a}_{b}" for a,b in sp_p.columns]
sp_p=sp_p.join(tax)
sp_p["log2FC_mucosal_vs_luminal_mean"]=np.log2((sp_p["mean_abundance_mucosal"]+1e-6)/(sp_p["mean_abundance_luminal"]+1e-6))
sp_p["prevalence_diff_mucosal_minus_luminal"]=sp_p["prevalence_mucosal"]-sp_p["prevalence_luminal"]
sp_p.sort_values("mean_abundance_mucosal",ascending=False).to_csv(OUT/"species_compartment_summary.tsv",sep="\t")
reg=gi.groupby(["Species","gi_region"]).abundance.agg(mean_abundance="mean",median_abundance="median",prevalence=lambda x:(x>0).mean()).reset_index()
reg_p=reg.pivot(index="Species",columns="gi_region",values=["mean_abundance","prevalence"])
reg_p.columns=[f"{a}_{b}" for a,b in reg_p.columns]
region_mean_cols=[x for x in reg_p if x.startswith("mean_abundance_")]
reg_p["max_region_mean_abundance"]=reg_p[region_mean_cols].max(axis=1)
reg_p["region_specificity"]=reg_p["max_region_mean_abundance"]/reg_p[region_mean_cols].sum(axis=1).replace(0,np.nan)
reg_p=reg_p.join(tax).sort_values("max_region_mean_abundance",ascending=False)
reg_p.to_csv(OUT/"species_region_summary.tsv",sep="\t")
D=squareform(pdist(wide.values,metric="braycurtis"))
A=-0.5*D**2
A=A-np.mean(A,axis=0,keepdims=True)-np.mean(A,axis=1,keepdims=True)+np.mean(A)
vals,vecs=np.linalg.eigh(A); order=np.argsort(vals)[::-1]; vals=vals[order]; vecs=vecs[:,order]
nz=vals>1e-10; coords=vecs[:,nz][:,:4]*np.sqrt(vals[nz][:4])
pcoa=pd.DataFrame(coords,index=wide.index,columns=[f"PCoA{i+1}" for i in range(coords.shape[1])])
pcoa=pd.concat([sample_meta,pcoa],axis=1); pcoa.to_csv(OUT/"bray_curtis_pcoa_gi.tsv",sep="\t")
pd.DataFrame({"eigenvalue":vals[:10],"explained_ratio":vals[:10]/max(vals.sum(),1e-12)},index=[f"axis{i+1}" for i in range(10)]).to_csv(OUT/"bray_curtis_pcoa_eigenvalues.tsv",sep="\t")
rng=np.random.default_rng(20261023)
pos_v=vecs[:,nz]*np.sqrt(vals[nz])
def permanova(groups,nperm=999):
    groups=np.asarray(groups)
    def pseudo_f(g):
        grand=pos_v.mean(axis=0); ssb=0; ssw=0
        for lv in np.unique(g):
            x=pos_v[g==lv]; mu=x.mean(axis=0)
            ssb += len(x)*np.sum((mu-grand)**2); ssw += np.sum((x-mu)**2)
        k=len(np.unique(g)); n=len(g)
        return (ssb/(k-1))/(ssw/(n-k)) if n>k else 0
    obs=pseudo_f(groups); perm=np.array([pseudo_f(rng.permutation(groups)) for _ in range(nperm)])
    return obs,(np.sum(perm>=obs)+1)/(nperm+1)
tests=[]
for var in ["compartment","monkey_id","gi_region"]:
    ok=sample_meta[var].notna().values
    if len(np.unique(sample_meta.loc[ok,var]))>1:
        obs,pv=permanova(sample_meta.loc[ok,var].values)
        tests.append({"variable":var,"n":int(ok.sum()),"pseudo_F":obs,"p_perm_999":pv})
pd.DataFrame(tests).to_csv(OUT/"permanova_bray_curtis_gi.tsv",sep="\t",index=False)
with open(OUT/"README_analysis.md","w") as out:
    out.write("# Monkey 2bRAD vendor species-level pre-analysis\n\n")
    out.write("Input: vendor-filtered BcgI species abundance tables for luminal/content and mucosa deliveries.\n")
    out.write("Scope: exploratory species-level QC only; no strain-level claims are made here.\n\n")
    out.write("## Sample counts\n\n```text\n"); out.write(counts.to_string(index=False)+"\n```\n\n")
    ctrl=long.groupby(["vendor_compartment","sample_category"]).sample_id.nunique().reset_index(name="n")
    out.write("## Controls / non-GI samples\n\n```text\n"); out.write(ctrl.to_string(index=False)+"\n```\n\n")
    out.write("## Alpha diversity summary\n\n```text\n")
    out.write(alpha.groupby("compartment")[["observed_species","shannon","simpson"]].agg(["count","mean","median","min","max"]).to_string()+"\n```\n\n")
    out.write("## PERMANOVA (Bray-Curtis, 999 permutations)\n\n```text\n"); out.write(pd.DataFrame(tests).to_string(index=False)+"\n```\n\n")
    cols=["Genus","Species","mean_abundance_luminal","mean_abundance_mucosal","prevalence_luminal","prevalence_mucosal","log2FC_mucosal_vs_luminal_mean"]
    out.write("## Most abundant species, by compartment\n\n```text\n")
    out.write(sp_p.sort_values("mean_abundance_luminal",ascending=False).head(30).reset_index()[cols].to_string()+"\n```\n\n")
    out.write("## Most mucosal-enriched prevalent species\n\n```text\n")
    out.write(sp_p[(sp_p.prevalence_luminal>=0.1)|(sp_p.prevalence_mucosal>=0.1)].sort_values(["log2FC_mucosal_vs_luminal_mean","prevalence_mucosal"],ascending=[False,False]).head(30).reset_index()[cols].to_string()+"\n```\n\n")
    prev_cols=[x for x in reg_p if x.startswith("prevalence_")]
    out.write("## Most region-enriched prevalent species\n\n```text\n")
    out.write(reg_p[(reg_p[prev_cols].max(axis=1)>=0.2)].head(40).to_string()+"\n```\n")
print("Wrote",OUT)
print("\nSAMPLE COUNTS"); print(counts.to_string(index=False))
ctrl=long.groupby(["vendor_compartment","sample_category"]).sample_id.nunique().reset_index(name="n")
print("\nCONTROLS / NON-GI"); print(ctrl.to_string(index=False))
print("\nPERMANOVA"); print(pd.DataFrame(tests).to_string(index=False))
cols=["Genus","Species","mean_abundance_luminal","mean_abundance_mucosal","prevalence_luminal","prevalence_mucosal","log2FC_mucosal_vs_luminal_mean"]
print("\nTOP LUMINAL/SHARED SPECIES"); print(sp_p.sort_values("mean_abundance_luminal",ascending=False).head(20)[cols].to_string(index=False))
print("\nMOST MUCOSAL-ENRICHED"); print(sp_p[(sp_p.prevalence_luminal>=0.1)|(sp_p.prevalence_mucosal>=0.1)].sort_values("log2FC_mucosal_vs_luminal_mean",ascending=False).head(20)[cols].to_string(index=False))
