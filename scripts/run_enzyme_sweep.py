import subprocess, re, os
BIN=os.environ.get("STRAIN2BSCAN_BIN","strain2bscan")
GEN,TRUTH,READS,WORK="acnes/genomes","acnes/truth","acnes/reads","enzsweep"
os.makedirs(WORK,exist_ok=True)
env={**os.environ,"STRAIN2BSCAN_THREADS":"16"}
LADDER={1:"BcgI",2:"BcgI,CspCI",4:"BcgI,CspCI,BaeI,AlfI",
        8:"BcgI,CspCI,BaeI,AlfI,AloI,BsaXI,CjeI,PpiI",
        14:"CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"}
# per-sample truth
stru={}
for i in range(1,6):
    rows=[]
    for l in open(f"{TRUTH}/sample{i}.truth.tsv"):
        if not l.startswith("#"): g,a=l.strip().split("\t"); rows.append((g,float(a)))
    stru[i]=rows
def timed(cmd):
    r=subprocess.run(["/usr/bin/time","-l",BIN]+cmd,capture_output=True,text=True,env=env)
    t=re.search(r'([0-9.]+) real',r.stderr); rss=re.search(r'(\d+)\s+maximum resident',r.stderr)
    return r.stdout, (float(t.group(1)) if t else 0), (int(rss.group(1))/1048576 if rss else 0)
TARGET="GCF_009737125.1"  # depth-sweep target strain
out=open(f"{WORK}/enzyme_sweep.tsv","w")
out.write("n_enzymes\tstrain_specific\tcluster_specific\tn_clusters\tbuild_s\tbuild_RSS_MB\tprofile_s\tprofile_RSS_MB\tprecision\trecall\tbray_curtis\tdetect_target_0.5x\tdetect_target_1x\n")
for k in sorted(LADDER):
    enz=LADDER[k]; db=f"{WORK}/db_{k}.tsv"
    so,bt,brss=timed(["cluster","--genomes",GEN,"--enzyme",enz,"--out",db,"--similarity","0.95"])
    ss=int(re.search(r'strain_specific=(\d+)',so).group(1)); cs=int(re.search(r'cluster_specific=(\d+)',so).group(1))
    ncl=int(re.search(r'into (\d+) cluster',so).group(1))
    g2c={}
    for l in open(db.replace(".tsv",".members.tsv")):
        if not l.startswith("#"): g,c=l.strip().split("\t"); g2c[g]=c
    TP=FP=FN=0; bcs=[]; ptimes=[]; prss=0
    for i in range(1,6):
        pred=f"{WORK}/s{i}_{k}.pred"
        po,pt,pr=timed(["profile","--db",db,"--reads",f"{READS}/sample{i}.fq","--out",pred])
        ptimes.append(pt); prss=max(prss,pr)
        ab={}
        for g,a in stru[i]: c=g2c.get(g,g); ab[c]=ab.get(c,0)+a
        tf=f"{WORK}/s{i}_{k}.truth"; open(tf,"w").write("".join(f"{c}\t{a}\n" for c,a in ab.items()))
        ev=subprocess.run([BIN,"evaluate","--pred",pred,"--truth",tf,"--present","0.01"],capture_output=True,text=True).stdout
        TP+=int(re.search(r'TP=(\d+)',ev).group(1)); FP+=int(re.search(r'FP=(\d+)',ev).group(1)); FN+=int(re.search(r'FN=(\d+)',ev).group(1))
        bcs.append(float(re.search(r'Bray-Curtis=([0-9.]+)',ev).group(1)))
    P=TP/(TP+FP) if TP+FP else 0; R=TP/(TP+FN) if TP+FN else 0; BC=sum(bcs)/len(bcs)
    # low-depth detection of TARGET at 0.5x and 1x
    tcl=g2c.get(TARGET,TARGET); det={}
    for dx in ["0.5","1.0"]:
        fq=f"depthsweep/cae_{dx}X.fq"
        r=subprocess.run([BIN,"profile","--db",db,"--reads",fq],capture_output=True,text=True,env=env).stdout
        det[dx]=1 if re.search(rf'^\s+{re.escape(tcl)}\s',r,re.M) else 0
    out.write(f"{k}\t{ss}\t{cs}\t{ncl}\t{bt:.2f}\t{brss:.0f}\t{sum(ptimes)/5:.2f}\t{prss:.0f}\t{P:.3f}\t{R:.3f}\t{BC:.3f}\t{det['0.5']}\t{det['1.0']}\n")
    print(f"{k} enzymes: strain-spec={ss} clusters={ncl} build={bt:.2f}s/{brss:.0f}MB profile={sum(ptimes)/5:.2f}s P={P:.2f} R={R:.2f} BC={BC:.3f} det@0.5x={det['0.5']} @1x={det['1.0']}",flush=True)
out.close(); print("done -> enzsweep/enzyme_sweep.tsv")
