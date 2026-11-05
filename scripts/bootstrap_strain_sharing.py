#!/usr/bin/env python3
"""Bootstrap contrasts for Strain2bScan strain-sharing rates."""
from pathlib import Path
import csv, random, statistics
ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/"results/strain_sharing_20261105"
OUT=ROOT/"results/strain_sharing_20261105/bootstrap_contrasts.tsv"
RNG=random.Random(20261105)

def read(path):
    with open(path,newline="") as f: return list(csv.DictReader(f,delimiter="\t"))

def boot(values_a,values_b,n=10000):
    vals=[]
    for _ in range(n):
        a=[RNG.choice(values_a) for _ in values_a]
        b=[RNG.choice(values_b) for _ in values_b]
        vals.append(statistics.mean(a)-statistics.mean(b))
    vals.sort()
    return vals[int(0.025*n)],vals[int(0.975*n)-1]

def contrast(rows,group_a,group_b,metric,label,n=10000):
    aa=[float(r[metric]) for r in rows if r["comparison_group"]==group_a and int(r["n_analyzable_species"])>0]
    bb=[float(r[metric]) for r in rows if r["comparison_group"]==group_b and int(r["n_analyzable_species"])>0]
    if not aa or not bb: return None
    lo,hi=boot(aa,bb,n)
    return {"comparison":label,"metric":metric,"group_a":group_a,"group_b":group_b,
            "n_a":len(aa),"n_b":len(bb),"mean_a":statistics.mean(aa),"mean_b":statistics.mean(bb),
            "mean_difference_A_minus_B":statistics.mean(aa)-statistics.mean(bb),
            "bootstrap_ci_low":lo,"bootstrap_ci_high":hi}

sal=read(IN/"saliva_pairwise_strain_sharing.tsv")
mon=read(IN/"monkey_pairwise_strain_sharing.tsv")
out=[]
# Saliva: within-subject temporal sharing versus all between-subject controls.
for metric in ["strain_sharing_rate_any_cluster","strain_sharing_rate_dominant_cluster"]:
    aa=[float(r[metric]) for r in sal if r["same_subject"]=="1"]
    bb=[float(r[metric]) for r in sal if r["same_subject"]=="0"]
    lo,hi=boot(aa,bb)
    out.append({"comparison":"saliva_within_subject_vs_between_subject","metric":metric,
                "group_a":"same_subject_between_timepoints","group_b":"all_between_subject",
                "n_a":len(aa),"n_b":len(bb),"mean_a":statistics.mean(aa),"mean_b":statistics.mean(bb),
                "mean_difference_A_minus_B":statistics.mean(aa)-statistics.mean(bb),
                "bootstrap_ci_low":lo,"bootstrap_ci_high":hi})
# Monkey: local positional stability versus between-monkey same-age same-site.
for metric in ["strain_sharing_rate_any_cluster","strain_sharing_rate_dominant_cluster"]:
    r=contrast(mon,"within_monkey_between_position_same_region",
               "between_monkey_same_age_same_region_compartment",metric,
               "monkey_within_monkey_positions_vs_between_monkey_same_site",10000)
    if r: out.append(r)
# Monkey old-young corresponding site: LI luminal versus LI mucosal.
li_lum=[r for r in mon if r["comparison_group"]=="between_age_corresponding_region_compartment"
        and r["region_a"]=="LI" and r["compartment_a"]=="luminal"]
li_muc=[r for r in mon if r["comparison_group"]=="between_age_corresponding_region_compartment"
        and r["region_a"]=="LI" and r["compartment_a"]=="mucosal"]
for metric in ["strain_sharing_rate_any_cluster","strain_sharing_rate_dominant_cluster"]:
    aa=[float(r[metric]) for r in li_lum if int(r["n_analyzable_species"])>0]
    bb=[float(r[metric]) for r in li_muc if int(r["n_analyzable_species"])>0]
    lo,hi=boot(aa,bb)
    out.append({"comparison":"monkey_old_young_LI_luminal_vs_LI_mucosal","metric":metric,
                "group_a":"LI_luminal","group_b":"LI_mucosal","n_a":len(aa),"n_b":len(bb),
                "mean_a":statistics.mean(aa),"mean_b":statistics.mean(bb),
                "mean_difference_A_minus_B":statistics.mean(aa)-statistics.mean(bb),
                "bootstrap_ci_low":lo,"bootstrap_ci_high":hi})
fields=["comparison","metric","group_a","group_b","n_a","n_b","mean_a","mean_b",
        "mean_difference_A_minus_B","bootstrap_ci_low","bootstrap_ci_high"]
with open(OUT,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
    w.writeheader(); w.writerows(out)
print(OUT)
