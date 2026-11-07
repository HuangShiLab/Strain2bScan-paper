#!/usr/bin/env python3
"""Time inStrain quick_profile on the existing 225 matched simulation BAMs."""
from __future__ import annotations
import csv, re, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
MATCHED=ROOT/"results/strainscan_rerun/run_manifest.tsv"
BAM=ROOT/"work/instrain_sim225/bam"; REF=ROOT/"work/instrain_sim225/references"
OUT=ROOT/"results/instrain_sim225"; WORK=Path("/tmp/instrain_sim225_timing")
INSTRAIN=ROOT/"tools/instrain191-venv-20261018/bin/instrain"

def matched():
 rows=[]
 for r in csv.DictReader(MATCHED.open(),delimiter="\t"):
  if r["kind"]!="single" or "_rep1_" not in r["sample"]: continue
  if r["species"]=="Mycobacterium_tuberculosis":
   if "__same_k" in r["sample"]: rows.append(r)
  elif "__diff_k" in r["sample"]: rows.append(r)
 if len(rows)!=225: raise SystemExit(f"expected 225 rows, got {len(rows)}")
 return rows

def parse_time(stderr):
 m=re.search(r"real\s+([0-9.]+)",stderr)
 rss_values=[int(x) for x in re.findall(r"^\s*(\d+)\s+",stderr,re.M)]
 return float(m.group(1)) if m else None,(max(rss_values)/1024 if rss_values else None)

def run_one(row):
 sample=row["sample"]; species=row["sample"].split("__",1)[0]
 bam=BAM/f"{sample}.bam"; stb=REF/f"{species}.stb"; fasta=REF/f"{species}.fna"
 out=WORK/sample
 cmd=[INSTRAIN,"quick_profile","-p","2","-s",str(stb),"-o",str(out),str(bam),str(fasta)]
 timed=["/usr/bin/time","-l",*cmd]
 proc=subprocess.run(timed,capture_output=True,text=True)
 wall,rss=parse_time(proc.stderr)
 if proc.returncode!=0: raise RuntimeError(f"{sample}: {proc.stderr[-1000:]}")
 depth=float(sample.rsplit("_d",1)[1])
 return {"sample":sample,"species":species,"depth":depth,"wall_s":wall,"maxrss_mb":rss,"returncode":proc.returncode}

def main():
 rows=matched(); WORK.mkdir(parents=True,exist_ok=True); records=[]
 with ThreadPoolExecutor(max_workers=4) as ex:
  futures={ex.submit(run_one,r):r for r in rows}
  for i,fut in enumerate(as_completed(futures),1):
   rec=fut.result(); records.append(rec); print(f"[{i}/225] {rec['sample']} wall={rec['wall_s']:.2f}s rss={rec['maxrss_mb']:.0f}MB",flush=True)
 fields=["sample","species","depth","wall_s","maxrss_mb","returncode"]
 with (OUT/"profiling_runtime.tsv").open("w",newline="") as fh:
  w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(records)
 print(OUT/"profiling_runtime.tsv")
if __name__=="__main__": main()
