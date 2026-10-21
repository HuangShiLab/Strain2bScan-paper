#!/usr/bin/env python3
"""Rerun the enzyme-count sweep on 225 matched single-species simulations.

Five cumulative enzyme ladders are applied to the same 15-species genome pool and the same
225 matched samples. Every prediction is scored against the Strain2bScan 0.95 database's own
truth-cluster partition, preserving one per-sample provenance row.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, os, re, subprocess, time, statistics
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

PAPER=Path(__file__).resolve().parents[1]
MATCHED=PAPER/"work/strainscan_headtohead_rerun/run_manifest.tsv"
POOL=PAPER/"figure_raw_data/sim_pool_manifest.tsv"
GENOMES=PAPER/"figure_raw_data/sim_genome_pool"
WORK=PAPER/"work/enzyme_sweep_225"
OUT=PAPER/"results/enzyme_sweep_225"
BINARY=PAPER.parent/"Strain2bScan/target/release/strain2bscan"
LADDER={
 1:"BcgI",
 2:"BcgI,CspCI",
 4:"BcgI,CspCI,BaeI,AlfI",
 8:"BcgI,CspCI,BaeI,AlfI,AloI,BsaXI,CjeI,PpiI",
 14:"CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI",
}
FIELDS=["species","sample","strategy","k","rep","depth","enzyme_count","enzyme_set",
        "n_truth_clusters","n_pred_clusters","TP","FP","FN","precision","recall","f1",
        "prediction","prediction_sha256","truth","database"]

def read_tsv(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))

def sha256(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
 return h.hexdigest()

def run(cmd,log):
 with log.open("a") as f:
  subprocess.run(cmd,check=True,stdout=f,stderr=subprocess.STDOUT)

def metrics(truth,pred):
 tp=len(truth&pred); fp=len(pred-truth); fn=len(truth-pred)
 p=tp/(tp+fp) if tp+fp else 0.0
 r=tp/(tp+fn) if tp+fn else 0.0
 f=2*p*r/(p+r) if p+r else 0.0
 return tp,fp,fn,p,r,f

def matched():
 rows=[]
 for r in read_tsv(MATCHED):
  if r["kind"]!="single" or "_rep1_" not in r["sample"]:continue
  if r["species"]=="Mycobacterium_tuberculosis":
   if "__same_k" in r["sample"]:rows.append(r)
  elif "__diff_k" in r["sample"]:rows.append(r)
 if len(rows)!=225:raise RuntimeError(f"expected 225 matched samples, got {len(rows)}")
 return rows

def build_db(binary,species,enzymes,count,work,threads):
 dbdir=work/"databases"/species/f"enzymes_{count}";dbdir.mkdir(parents=True,exist_ok=True)
 db=dbdir/"db.tsv";members=db.with_suffix(".members.tsv");log=work/"logs"/f"build_{species}_{count}.log"
 if db.exists() and members.exists() and db.stat().st_size:
  return db,members
 run([str(binary),"cluster","--genomes",str(GENOMES/species),"--enzyme",enzymes,
      "--similarity","0.95","--out",str(db)],log)
 if not db.exists() or not members.exists():raise RuntimeError(f"database build failed for {species}, {enzymes}")
 return db,members

def parse_calls(path):
 return {line.split("\t",1)[0] for line in path.read_text().splitlines() if line and not line.startswith("#")}

def profile_one(binary,db,row,out,log,threads):
 out.parent.mkdir(parents=True,exist_ok=True)
 if not out.exists():
  fifo=out.with_suffix(".fq")
  if fifo.exists():os.unlink(fifo)
  os.mkfifo(fifo)
  env={**os.environ,"STRAIN2BSCAN_THREADS":str(threads)}
  proc=subprocess.Popen([str(binary),"profile","--db",str(db),"--reads",str(fifo),"--out",str(out)],
    stdout=log.open("a"),stderr=subprocess.STDOUT,env=env)
  with fifo.open("wb") as stream:
   writer=subprocess.Popen(["gzip","-dc",row["r1"],row["r2"]],stdout=stream,stderr=subprocess.DEVNULL)
   writer.wait()
  code=proc.wait()
  fifo.unlink(missing_ok=True)
  if code:raise subprocess.CalledProcessError(code,proc.args)
 return parse_calls(out)

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--threads",type=int,default=4)
 ap.add_argument("--workers",type=int,default=4)
 args=ap.parse_args()
 WORK.mkdir(parents=True,exist_ok=True);(WORK/"logs").mkdir(exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
 rows=matched()
 if not BINARY.exists():raise FileNotFoundError(BINARY)
 rows=matched()
 pool=read_tsv(POOL);truth_cluster={(r["species"],r["accession"]):r["cluster"] for r in pool}
 byspecies=defaultdict(list)
 for row in rows:byspecies[row["species"]].append(row)
 start=time.time();dbprov=[];futures={};results=[]
 with ThreadPoolExecutor(max_workers=args.workers) as ex:
  for species in sorted(byspecies):
   for count,enzymes in LADDER.items():
    db,members=build_db(BINARY,species,enzymes,count,WORK,args.threads)
    text=(WORK/"logs"/f"build_{species}_{count}.log").read_text()
    def marker(pattern,default=""):
     m=re.search(pattern,text);return int(m.group(1)) if m else default
    dbprov.append({"species":species,"enzyme_count":count,"enzyme_set":enzymes,
      "strain_specific_markers":marker(r"strain_specific=(\d+)"),
      "cluster_specific_markers":marker(r"cluster_specific=(\d+)"),
      "n_clusters":marker(r"into (\d+) cluster"),
      "database":str(db),"database_sha256":sha256(db),"members":str(members)})
    g2c={}
    for line in members.read_text().splitlines():
     if not line.startswith("#"):
      genome,cluster=line.rstrip("\n").split("\t");g2c[genome]=cluster
    for row in byspecies[species]:
     sample=row["sample"]
     out=WORK/"predictions"/f"enzymes_{count}"/f"{sample}.pred"
     log=WORK/"logs"/f"profile_{count}_{sample}.log"
     truth={truth_cluster[(species,line.split("\t")[0])]
            for line in Path(row["truth"]).read_text().splitlines()
            if line and not line.startswith("#")}
     fut=ex.submit(profile_one,BINARY,db,row,out,log,args.threads)
     futures[fut]=(species,row,count,enzymes,truth,g2c,db)
  for done,fut in enumerate(as_completed(futures),1):
   species,row,count,enzymes,truth,g2c,db=futures[fut]
   pred=fut.result()
   tp,fp,fn,p,r,f=metrics(truth,pred)
   parts=row["sample"].split("__",1)[1].split("_")
   results.append({"species":species,"sample":row["sample"],"strategy":parts[0],
     "k":parts[1].removeprefix("k"),"rep":parts[2].removeprefix("rep"),"depth":parts[3].removeprefix("d"),
     "enzyme_count":count,"enzyme_set":enzymes,"n_truth_clusters":len(truth),"n_pred_clusters":len(pred),
     "TP":tp,"FP":fp,"FN":fn,"precision":p,"recall":r,"f1":f,
     "prediction":str(out),"prediction_sha256":sha256(out),"truth":row["truth"],"database":str(db)})
   print(f"[{done}/{len(futures)}] {sample} enzymes={count} P={p:.3f} R={r:.3f} F1={f:.3f}",flush=True)
 fields=FIELDS
 with (OUT/"persample.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
  w.writeheader();w.writerows(results)
 with (OUT/"database_provenance.tsv").open("w",newline="") as f:
  fields=list(dbprov[0]);w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
  w.writeheader();w.writerows(dbprov)
 # Aggregate medians by species and enzyme count.
 agg=defaultdict(list)
 for r in results:agg[(r["species"],r["enzyme_count"])].append(r)
 with (OUT/"by_species_enzyme.tsv").open("w",newline="") as f:
  fields=["species","enzyme_count","n_samples","median_precision","median_recall","median_f1",
          "mean_precision","mean_recall","mean_f1"]
  w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader()
  for (species,count),rs in sorted(agg.items()):
   vals=[tuple(float(x[m]) for m in ("precision","recall","f1")) for x in rs]
   med=[statistics.median(v[i] for v in vals) for i in range(3)]
   mean=[sum(v[i] for v in vals)/len(vals) for i in range(3)]
   w.writerow({"species":species,"enzyme_count":count,"n_samples":len(rs),
    **{f"median_{m}":med[i] for i,m in enumerate(("precision","recall","f1"))},
    **{f"mean_{m}":mean[i] for i,m in enumerate(("precision","recall","f1"))}})
 prov={"binary":str(BINARY),"binary_sha256":sha256(BINARY),"enzyme_ladder":LADDER,
       "n_samples_per_ladder":len(rows),"n_results":len(results),"worker_threads":args.workers,
       "profile_threads":args.threads,"elapsed_s":time.time()-start,
       "scoring":"Strain2bScan own 0.95-cluster space; truth clusters from sim_pool_manifest.tsv"}
 (OUT/"provenance.json").write_text(json.dumps(prov,indent=2)+"\n")
 print(prov,flush=True)
if __name__=="__main__":main()
