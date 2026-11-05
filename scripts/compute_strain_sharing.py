#!/usr/bin/env python3
"""Calculate pairwise strain-sharing rates from Strain2bScan cluster calls.

Primary, StrainPhlAn-compatible metric:
  shared dominant 0.95 strain cluster / number of species called in both samples.

Coexistence-aware metrics:
  - any shared cluster / species called in both samples
  - mean per-species cluster-set Jaccard / species called in both samples
"""
from pathlib import Path
from itertools import combinations
import csv, math, re, statistics

ROOT=Path(__file__).resolve().parents[1]
MONKEY=ROOT/"results/monkey_spatial_148/hpc_outputs/analysis/full_148_20261027"
SALIVA=ROOT/"results/wms32_rebuild/native20_matched/native_strain_calls.tsv"
OUT=ROOT/"results/strain_sharing_20261105"
OUT.mkdir(parents=True,exist_ok=True)

def read_tsv(path):
    with open(path,newline="") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    # Normalize the optional leading "#" used by Strain2bScan output headers.
    return [{k.lstrip("#"):v for k,v in r.items()} for r in rows]

def write_tsv(path,fields,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with open(path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerows(rows)

def clusters(x):
    return set(re.split(r"\|",str(x).strip()))

def pair_metrics(calls_a,calls_b):
    """calls_* = {species: list[(cluster_label, abundance)]}."""
    common=sorted(set(calls_a)&set(calls_b))
    n_any=0; n_dom=0; jac_values=[]; shared_clusters=0; union_clusters=0
    shared_details=[]
    for sp in common:
        a=calls_a[sp]; b=calls_b[sp]
        ca={x[0] for x in a}; cb={x[0] for x in b}
        inter=ca&cb; union=ca|cb
        dom_a=max(a,key=lambda x:float(x[1]))[0]
        dom_b=max(b,key=lambda x:float(x[1]))[0]
        is_any=bool(inter); is_dom=dom_a==dom_b
        n_any += int(is_any); n_dom += int(is_dom)
        jac_values.append(len(inter)/len(union) if union else 0)
        shared_clusters += len(inter); union_clusters += len(union)
        shared_details.append({"species":sp,"clusters_a":"|".join(sorted(ca)),
                               "clusters_b":"|".join(sorted(cb)),
                               "dominant_a":dom_a,"dominant_b":dom_b,
                               "shared":int(is_any),"dominant_shared":int(is_dom),
                               "cluster_jaccard":len(inter)/len(union) if union else 0})
    denom=len(common)
    return {"n_analyzable_species":denom,
            "n_shared_species_any_cluster":n_any,
            "n_shared_species_dominant_cluster":n_dom,
            "strain_sharing_rate_any_cluster":(n_any/denom if denom else None),
            "strain_sharing_rate_dominant_cluster":(n_dom/denom if denom else None),
            "mean_per_species_cluster_jaccard":(statistics.mean(jac_values) if jac_values else None),
            "shared_cluster_intersections":shared_clusters,
            "union_cluster_labels":union_clusters,
            "details":shared_details}

def summarize(rows,metrics=("strain_sharing_rate_any_cluster","strain_sharing_rate_dominant_cluster","mean_per_species_cluster_jaccard")):
    out=[]
    grouped={}
    for r in rows:
        grouped.setdefault(r["comparison_group"],[]).append(r)
    for group,rs in grouped.items():
        rr=[r for r in rs if r["n_analyzable_species"]>0]
        rec={"comparison_group":group,"n_pairs":len(rs),"n_pairs_analyzable":len(rr),
             "n_pairs_with_shared_strain":sum(int(r["n_shared_species_any_cluster"]>0) for r in rr)}
        for m in metrics:
            vals=[float(r[m]) for r in rr if r[m]!=""]
            rec[f"{m}_mean"]=statistics.mean(vals) if vals else ""
            rec[f"{m}_median"]=statistics.median(vals) if vals else ""
        rec["pooled_shared_species"]=sum(int(r["n_shared_species_any_cluster"]) for r in rr)
        rec["pooled_analyzable_species_pairs"]=sum(int(r["n_analyzable_species"]) for r in rr)
        rec["pooled_strain_sharing_rate"]=(rec["pooled_shared_species"]/rec["pooled_analyzable_species_pairs"]
                                            if rec["pooled_analyzable_species_pairs"] else "")
        out.append(rec)
    return sorted(out,key=lambda x:x["comparison_group"])

# ---------- saliva ----------
sal_rows=read_tsv(SALIVA)
sal={}
for r in sal_rows:
    sample=r["sample"].lstrip("#")
    species=r["species"]; cluster=r["cluster"]; abundance=float(r["abundance"])
    sal.setdefault(sample,{}).setdefault(species,[]).append((cluster,abundance))
saliva_rows=[]
for a,b in combinations(sorted(sal),2):
    def sample_info(s):
        # Format is S<timepoint>-<subject>, e.g. S11-1.
        stem=s.split(".")[0]
        time=stem.split("-",1)[0].lstrip("S")
        subject=stem.split("-",1)[1]
        return time,subject
    ta,sa=sample_info(a); tb,sb=sample_info(b)
    if sa==sb:
        group="same_subject_between_timepoints"
    elif ta==tb:
        group="same_timepoint_between_subjects"
    else:
        group="different_subject_between_timepoints"
    m=pair_metrics(sal[a],sal[b]); m.pop("details")
    saliva_rows.append({"sample_a":a,"sample_b":b,"subject_a":sa,"subject_b":sb,
                        "timepoint_a":ta,"timepoint_b":tb,"comparison_group":group,
                        "same_subject":int(sa==sb),"same_timepoint":int(ta==tb),**m})
sal_fields=["sample_a","sample_b","subject_a","subject_b","timepoint_a","timepoint_b","comparison_group",
            "same_subject","same_timepoint","n_analyzable_species","n_shared_species_any_cluster",
            "n_shared_species_dominant_cluster","strain_sharing_rate_any_cluster",
            "strain_sharing_rate_dominant_cluster","mean_per_species_cluster_jaccard",
            "shared_cluster_intersections","union_cluster_labels"]
write_tsv(OUT/"saliva_pairwise_strain_sharing.tsv",sal_fields,saliva_rows)
write_tsv(OUT/"saliva_strain_sharing_group_summary.tsv",
          ["comparison_group","n_pairs","n_pairs_analyzable","n_pairs_with_shared_strain",
           "strain_sharing_rate_any_cluster_mean","strain_sharing_rate_any_cluster_median",
           "strain_sharing_rate_dominant_cluster_mean","strain_sharing_rate_dominant_cluster_median",
           "mean_per_species_cluster_jaccard_mean","mean_per_species_cluster_jaccard_median",
           "pooled_shared_species","pooled_analyzable_species_pairs","pooled_strain_sharing_rate"],
          summarize(saliva_rows))

# ---------- monkey ----------
mon_rows=read_tsv(MONKEY/"full148_2brad_strain_calls_long.tsv")
meta={r["sample_id"]:r for r in read_tsv(MONKEY/"full_148_sample_metadata.tsv")}
mon={}
for r in mon_rows:
    sid=r["sample_id"]
    if meta.get(sid,{}).get("sample_group")!="GI_spatial": continue
    species=r["species"]; cluster=r["cluster"]; abundance=float(r["abundance"])
    mon.setdefault(sid,{}).setdefault(species,[]).append((cluster,abundance))
monkey_rows=[]
for a,b in combinations(sorted(mon),2):
    ma,mb=meta[a],meta[b]
    same_monkey=ma["monkey_id"]==mb["monkey_id"]
    same_region=ma["gi_region"]==mb["gi_region"]
    same_comp=ma["compartment"]==mb["compartment"]
    same_pos=ma["position_cm"]==mb["position_cm"]
    age_a="old" if ma["monkey_id"].startswith("O") else "young"
    age_b="old" if mb["monkey_id"].startswith("O") else "young"
    same_age=age_a==age_b
    if same_monkey and same_region and same_comp and not same_pos:
        group="within_monkey_between_position_same_region"
    elif same_monkey and same_region and not same_comp:
        group="within_monkey_mucosal_luminal_same_region"
    elif same_monkey:
        group="within_monkey_between_region"
    elif same_region and same_comp and same_age:
        group="between_monkey_same_age_same_region_compartment"
    elif same_region and same_comp and not same_age:
        group="between_age_corresponding_region_compartment"
    elif same_comp:
        group="between_monkey_between_region_same_compartment"
    else:
        group="between_compartment"
    m=pair_metrics(mon[a],mon[b]); m.pop("details")
    monkey_rows.append({"sample_a":a,"sample_b":b,"monkey_a":ma["monkey_id"],"monkey_b":mb["monkey_id"],
                        "age_a":age_a,"age_b":age_b,"region_a":ma["gi_region"],"region_b":mb["gi_region"],
                        "compartment_a":ma["compartment"],"compartment_b":mb["compartment"],
                        "position_a":ma["position_cm"],"position_b":mb["position_cm"],
                        "comparison_group":group,"same_monkey":int(same_monkey),"same_region":int(same_region),
                        "same_compartment":int(same_comp),"same_age":int(same_age),
                        "same_position":int(same_pos),**m})
mon_fields=["sample_a","sample_b","monkey_a","monkey_b","age_a","age_b","region_a","region_b",
            "compartment_a","compartment_b","position_a","position_b","comparison_group",
            "same_monkey","same_region","same_compartment","same_age","same_position",
            "n_analyzable_species","n_shared_species_any_cluster","n_shared_species_dominant_cluster",
            "strain_sharing_rate_any_cluster","strain_sharing_rate_dominant_cluster",
            "mean_per_species_cluster_jaccard","shared_cluster_intersections","union_cluster_labels"]
write_tsv(OUT/"monkey_pairwise_strain_sharing.tsv",mon_fields,monkey_rows)
write_tsv(OUT/"monkey_strain_sharing_group_summary.tsv",
          ["comparison_group","n_pairs","n_pairs_analyzable","n_pairs_with_shared_strain",
           "strain_sharing_rate_any_cluster_mean","strain_sharing_rate_any_cluster_median",
           "strain_sharing_rate_dominant_cluster_mean","strain_sharing_rate_dominant_cluster_median",
           "mean_per_species_cluster_jaccard_mean","mean_per_species_cluster_jaccard_median",
           "pooled_shared_species","pooled_analyzable_species_pairs","pooled_strain_sharing_rate"],
          summarize(monkey_rows))

# age-corresponding-site compact summary
age=[]
for r in monkey_rows:
    if r["comparison_group"]=="between_age_corresponding_region_compartment":
        age.append(r)
summary=[]
for group in sorted({f"{r['region_a']}|{r['compartment_a']}" for r in age}):
    region,comp=group.split("|")
    rr=[r for r in age if r["region_a"]==region and r["compartment_a"]==comp]
    rec={"region":region,"compartment":comp,"n_pairs":len(rr),
         "n_pairs_analyzable":sum(int(r["n_analyzable_species"])>0 for r in rr)}
    for m in ["strain_sharing_rate_any_cluster","strain_sharing_rate_dominant_cluster",
              "mean_per_species_cluster_jaccard"]:
        vals=[float(r[m]) for r in rr if r["n_analyzable_species"]>0]
        rec[f"{m}_mean"]=statistics.mean(vals) if vals else ""
        rec[f"{m}_median"]=statistics.median(vals) if vals else ""
    rec["pooled_shared_species"]=sum(int(r["n_shared_species_any_cluster"]) for r in rr)
    rec["pooled_analyzable_species_pairs"]=sum(int(r["n_analyzable_species"]) for r in rr)
    rec["pooled_strain_sharing_rate"]=(rec["pooled_shared_species"]/rec["pooled_analyzable_species_pairs"]
                                       if rec["pooled_analyzable_species_pairs"] else "")
    summary.append(rec)
write_tsv(OUT/"monkey_between_age_corresponding_site_summary.tsv",
          ["region","compartment","n_pairs","n_pairs_analyzable",
           "strain_sharing_rate_any_cluster_mean","strain_sharing_rate_any_cluster_median",
           "strain_sharing_rate_dominant_cluster_mean","strain_sharing_rate_dominant_cluster_median",
           "mean_per_species_cluster_jaccard_mean","mean_per_species_cluster_jaccard_median",
           "pooled_shared_species","pooled_analyzable_species_pairs","pooled_strain_sharing_rate"],summary)

readme='''# Strain-sharing analysis

Inputs are already-gated Strain2bScan 0.95 strain-cluster calls.

For each sample pair and each species detected in both samples:
- `any_cluster`: at least one identical cluster label was observed in both samples.
- `dominant_cluster`: the highest within-species abundance cluster was identical in both samples.
- The pair-level primary sharing rate is `n_shared_species_dominant_cluster / n_analyzable_species`.
- A coexistence-aware secondary rate is `n_shared_species_any_cluster / n_analyzable_species`.
- Per-species cluster-set Jaccard is also reported for multi-strain species.

Saliva comparisons are within-subject timepoint contrasts and between-subject controls.
Monkey O1/O2 are labelled old and Y1/Y2 young; GI region and compartment are used for
corresponding-site contrasts. Age assignment should be checked against the final sample metadata.
'''
(OUT/"README.md").write_text(readme)
print("wrote",OUT)
