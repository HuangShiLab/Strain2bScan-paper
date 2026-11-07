#!/usr/bin/env python3
"""Build FracMinHash (k-mer) databases and profile the 225 matched simulations."""
from __future__ import annotations
import argparse,csv,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BIN=Path("/Users/macstudio/Downloads/Strain2bScan/target/release/strain2bscan")
GENOMES=ROOT/"archive/large_intermediates/figure_raw_data/sim_genome_pool"
READS=ROOT/"archive/large_intermediates/figure_raw_data/sim_single_species"
MATCHED=ROOT/"results/strainscan_rerun/run_manifest.tsv"
POOL=ROOT/"figure_raw_data/sim_pool_manifest.tsv"
OUT=ROOT/"results/kmer_sim225"; WORK=ROOT/"work/kmer_sim225"
def read(p): return list(csv.DictReader(open(p,newline=""),delimiter="\t"))
def matched():
 rows=[]
 for r in read(MATCHED):
  if r["kind"]!="single" or "_rep1_" not in r["sample"]: continue
  if r["species"]=="Mycobacterium_tuberculosis":
   if "__same_k" in r["sample"]: rows.append(r)
  elif "__diff_k" in r["sample"]: rows.append(r)
 if len(rows)!=225: raise SystemExit(len(rows))
 return rows
def run(cmd,log):
 with log.open("w") as h:
  subprocess.run(cmd,check=True,stdout=h,stderr=subprocess.STDOUT)
def main():
 ap=argparse.ArgumentParser(); ap.add_argument("--sketch-scale",type=int,default=100); ap.add_argument("--jobs",type=int,default=4); args=ap.parse_args()
 WORK.mkdir(parents=True,exist_ok=True); (WORK/"logs").mkdir(exist_ok=True); (WORK/"dbs").mkdir(exist_ok=True); (WORK/"predictions").mkdir(exist_ok=True)
 rows=matched(); species=sorted({r["species"] for r in rows})
 start=time.time()
 for sp in species:
  db=WORK/"dbs"/f"{sp}.kmer{args.sketch_scale:03d}.tsv"; log=WORK/"logs"/f"build_{sp}.log"
  if not db.exists(): run([BIN,"cluster","--genomes",str(GENOMES/sp),"--marker-source","kmer","--kmer-size","31","--sketch-scale",str(args.sketch_scale),"--similarity","0.95","--out",str(db)],log)
 print("databases complete",flush=True)
 gc={r["accession"]:r["cluster"] for r in read(POOL)}
 results=[]
 for i,row in enumerate(rows,1):
  sample=row["sample"]; sp=row["species"]; db=WORK/"dbs"/f"{sp}.kmer{args.sketch_scale:03d}.tsv"
  pred=WORK/"predictions"/f"{sample}.pred"; log=WORK/"logs"/f"profile_{sample}.log"
  if not pred.exists():
   tmp=WORK/f"{sample}.fq.gz"
   if not tmp.exists():
    with open(row["r1"],"rb") as a, open(row["r2"],"rb") as b, open(tmp,"wb") as o:
     for chunk in iter(lambda:a.read(1024*1024),b""): o.write(chunk)
     for chunk in iter(lambda:b.read(1024*1024),b""): o.write(chunk)
   t=time.time(); run([BIN,"profile","--db",str(db),"--reads",str(tmp),"--marker-source","kmer","--out",str(pred)],log); wall=time.time()-t
  else: wall=None
  truth={line.split("\t")[1] for line in Path(row["truth"]).read_text().splitlines()[1:] if line.strip()}
  pred_clusters=set()
  for line in pred.read_text().splitlines()[1:]: pred_clusters.add(line.split("\t")[0])
  tp=len(truth&pred_clusters); fp=len(pred_clusters-truth); fn=len(truth-pred_clusters)
  precision=tp/(tp+fp) if tp+fp else 0; recall=tp/(tp+fn) if tp+fn else 0; f1=2*precision*recall/(precision+recall) if precision+recall else 0
  depth=float(sample.rsplit("_d",1)[1]); results.append({"sample":sample,"species":sp,"depth":depth,"n_truth":len(truth),"n_pred":len(pred_clusters),"TP":tp,"FP":fp,"FN":fn,"precision":precision,"recall":recall,"f1":f1,"wall_s":wall})
  if i%25==0 or i==1: print(f"{i}/225 {sample} F1={f1:.3f}",flush=True)
 OUT.mkdir(parents=True,exist_ok=True)
 fields=list(results[0])
 with (OUT/"persample.tsv").open("w",newline="") as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(results)
 (OUT/"provenance.json").write_text(f'{{"marker_source":"kmer","kmer_size":31,"sketch_scale":{args.sketch_scale},"n_runs":{len(results)}}}\n')
 print(OUT/"persample.tsv",flush=True)
if __name__=="__main__": main()
