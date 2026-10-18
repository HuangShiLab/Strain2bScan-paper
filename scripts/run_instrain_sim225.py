#!/usr/bin/env python3
"""Run inStrain on the 225 matched single-species simulation samples and score in 0.95 cluster space."""
from __future__ import annotations
import argparse, csv, json, os, shutil, subprocess, time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

PAPER=Path(__file__).resolve().parents[1]
MATCHED=PAPER/"work/strainscan_headtohead_rerun/run_manifest.tsv"
POOL=PAPER/"figure_raw_data/sim_pool_manifest.tsv"
READS=PAPER/"archive/large_intermediates/figure_raw_data/sim_single_species"
WORK=PAPER/"work/instrain_sim225"
OUT=PAPER/"results/instrain_sim225"
FASTA=PAPER/"archive/large_intermediates/figure_raw_data/sim_genome_pool"
INSTRAIN=PAPER/"tools/instrain191-venv-20261018/bin/instrain"
BOWTIE2="bowtie2-align-s"; BOWTIE2_BUILD="bowtie2-build-s"; SAMTOOLS="samtools"

def read_tsv(p): return list(csv.DictReader(open(p),delimiter="\t"))
def matched():
 rows=[]
 for r in read_tsv(MATCHED):
  if r["kind"]!="single" or "_rep1_" not in r["sample"]: continue
  if r["species"]=="Mycobacterium_tuberculosis":
   if "__same_k" in r["sample"]: rows.append(r)
  elif "__diff_k" in r["sample"]: rows.append(r)
 if len(rows)!=225: raise SystemExit(f"expected 225 got {len(rows)}")
 return rows

def build_reference(species, log):
 refdir=WORK/"references"; refdir.mkdir(parents=True,exist_ok=True)
 fa=refdir/f"{species}.fna"; stb=refdir/f"{species}.stb"; index=refdir/species
 # Rebuild atomically if the FASTA is absent or contains duplicate headers (a failure mode from
 # the first concurrent implementation).
 def fasta_ok():
  if not fa.exists(): return False
  seen=set()
  for line in fa.open():
   if line.startswith(">"):
    h=line[1:].strip().split()[0]
    if h in seen: return False
    seen.add(h)
  return bool(seen)
 index_ok=all((refdir/f"{species}.{suffix}").exists() for suffix in ("1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"))
 if fa.exists() and fasta_ok() and index_ok:
  return fa,stb,index
 with log.open("a") as _lh: _lh.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} REBUILD_REFERENCE {species}\n")
 tmpfa=refdir/f"{species}.fna.tmp"; tmpstb=refdir/f"{species}.stb.tmp"
 with tmpfa.open("w") as out, tmpstb.open("w") as st:
  for ordinal,genome in enumerate(sorted((FASTA/species).glob("*.fna")),1):
   acc=genome.stem
   with genome.open() as gh:
    contig=0
    for line in gh:
     if line.startswith(">"):
      contig+=1
      sid=line[1:].strip().split()[0]
      name=f"{acc}_{ordinal}_{contig}__{sid}"
      out.write(f">{name}\n"); st.write(f"{name}\t{acc}\n")
     else: out.write(line)
 tmpfa.replace(fa); tmpstb.replace(stb)
 for suffix in ("1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2",
                "1.bt2l","2.bt2l","3.bt2l","4.bt2l","rev.1.bt2l","rev.2.bt2l"):
  old=refdir/f"{species}.{suffix}"
  if old.exists(): old.unlink()
 with log.open("a") as lh:
  subprocess.run([BOWTIE2_BUILD,"--threads","8",str(fa),str(index)],check=True,stdout=lh,stderr=subprocess.STDOUT)
 return fa,stb,index

def detect_clusters(path, genome_cluster):
 det=set()
 if not path.exists(): return det
 for row in csv.DictReader(path.open()):
  genome=row.get("genome","")
  try: cov=float(row.get("coverage",0)); breadth=float(row.get("breadth",0))
  except ValueError: continue
  if cov>=0.1 and breadth>=0.05 and genome in genome_cluster: det.add(genome_cluster[genome])
 return det

def metrics(truth,pred):
 tp=len(truth&pred); fp=len(pred-truth); fn=len(truth-pred)
 p=tp/(tp+fp) if tp+fp else 0; r=tp/(tp+fn) if tp+fn else 0; f=2*p*r/(p+r) if p+r else 0
 return tp,fp,fn,p,r,f

def run_one(row,genome_cluster,log):
 sample=row["sample"]; species=row["species"]
 bam=WORK/"bam"/f"{sample}.bam"; out=WORK/"instrain"/f"{sample}.instrain"
 bam.parent.mkdir(parents=True,exist_ok=True); out.parent.mkdir(parents=True,exist_ok=True)
 fa,stb,index=build_reference(species,log)
 cov=out/"genomeCoverage.csv"
 if not cov.exists(): cov=Path(str(out)+".genomeCoverage.csv")
 if cov.exists() and bam.exists():
  with log.open("a") as _lh: _lh.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} SKIP {sample}\n")
  truth=set()
  for line in Path(row["truth"]).read_text().splitlines():
   if not line.startswith("#") and line.strip():
    fields=line.split("\t"); truth.add(fields[1])
  pred=detect_clusters(cov,genome_cluster)
  tp,fp,fn,p,r,f=metrics(truth,pred)
  return {"sample":sample,"species":species,"n_truth":len(truth),"n_pred":len(pred),"TP":tp,"FP":fp,"FN":fn,"precision":p,"recall":r,"f1":f}
 if not bam.exists():
  cmd=f'{BOWTIE2} -x "{index}" -1 "{row["r1"]}" -2 "{row["r2"]}" -p 2 | {SAMTOOLS} sort -@ 2 -o "{bam}"'
  with log.open("a") as lh:
   subprocess.run(cmd,shell=True,check=True,stdout=lh,stderr=subprocess.STDOUT)
 if not out.with_suffix("").with_name(out.name+".genomeCoverage.csv").exists() and not Path(str(out)+".genomeCoverage.csv").exists():
  with log.open("a") as lh:
   subprocess.run([INSTRAIN,"quick_profile","-p","2","-s",str(stb),"-o",str(out),str(bam),str(fa)],check=True,stdout=lh,stderr=subprocess.STDOUT)
 cov=out/"genomeCoverage.csv"
 if not cov.exists(): cov=Path(str(out)+".genomeCoverage.csv")
 truth=set()
 for line in Path(row["truth"]).read_text().splitlines():
  if not line.startswith("#") and line.strip():
   fields=line.split("\t"); truth.add(fields[1])
 pred=detect_clusters(cov,genome_cluster)
 tp,fp,fn,p,r,f=metrics(truth,pred)
 return {"sample":sample,"species":species,"n_truth":len(truth),"n_pred":len(pred),"TP":tp,"FP":fp,"FN":fn,"precision":p,"recall":r,"f1":f}

def main():
 ap=argparse.ArgumentParser(); ap.add_argument("--workers",type=int,default=3); args=ap.parse_args()
 WORK.mkdir(parents=True,exist_ok=True); (WORK/"logs").mkdir(exist_ok=True); OUT.mkdir(parents=True,exist_ok=True)
 rows=matched(); pool=read_tsv(POOL); gc={r["accession"]:r["cluster"] for r in pool}
 start=time.time(); results=[]
 print('[references] building/validating all species indexes serially',flush=True)
 for species in sorted({r["species"] for r in rows}):
  build_reference(species,WORK/"logs"/f"reference_{species}.log")
 def task(row):
  sample=row["sample"]; log=WORK/"logs"/f"{sample}.log"
  try:
   return run_one(row,gc,log)
  except Exception as exc:
   with log.open("a") as lh: lh.write(f"ERROR {sample}: {exc}\n")
   return {"sample":sample,"species":row["species"],"n_truth":0,"n_pred":0,"TP":0,"FP":0,"FN":0,"precision":0.0,"recall":0.0,"f1":0.0}
 with ThreadPoolExecutor(max_workers=args.workers) as ex:
  futures={ex.submit(task,r):r for r in rows}
  for i,fut in enumerate(as_completed(futures),1):
   results.append(fut.result()); print(f"[{i}/225] {results[-1]['sample']} P={results[-1]['precision']:.3f} R={results[-1]['recall']:.3f} F1={results[-1]['f1']:.3f}",flush=True)
 cols=["sample","species","n_truth","n_pred","TP","FP","FN","precision","recall","f1"]
 with (OUT/"persample.tsv").open("w",newline="") as fh:
  w=csv.DictWriter(fh,fieldnames=cols,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(results)
 prov={"tool":"inStrain","version":subprocess.check_output([INSTRAIN,"--version"],text=True).strip(),"n_runs":len(results),"worker_threads":args.workers,"elapsed_s":time.time()-start}
 (OUT/"provenance.json").write_text(json.dumps(prov,indent=2)+"\n")
 print(prov)
if __name__=="__main__": main()
