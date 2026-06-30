import os, glob, subprocess, re, shutil, sys
sys.path.insert(0,'.')
from degrade import degrade, read_seq
BIN=os.environ.get("STRAIN2BSCAN_BIN","strain2bscan")
SET="CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"
GEN,TRUTH,READS,WORK="acnes/genomes","acnes/truth","acnes/reads","refqual"
os.makedirs(WORK,exist_ok=True)
env={**os.environ,"STRAIN2BSCAN_THREADS":"16"}
# truth strains + per-sample truth
truth_acc=set(); sample_truth={}
for i in range(1,6):
    rows=[]
    for l in open(f"{TRUTH}/sample{i}.truth.tsv"):
        if l.startswith("#"): continue
        g,ab=l.strip().split("\t"); truth_acc.add(g); rows.append((g,float(ab)))
    sample_truth[i]=rows
all_g=[os.path.basename(p)[:-4] for p in glob.glob(f"{GEN}/*.fna")]
background=[g for g in all_g if g not in truth_acc]
contaminant=read_seq(sorted(glob.glob("multispecies/genomes/Escherichia_coli/*.fna"))[0])
levels=[(1.00,0.00,1),(0.95,0.01,20),(0.90,0.02,50),(0.80,0.05,100),(0.70,0.08,200),(0.50,0.10,400)]
def write_fasta(frags,path):
    with open(path,"w") as f:
        for i,fr in enumerate(frags):
            f.write(f">ctg{i}\n")
            for j in range(0,len(fr),70): f.write(fr[j:j+70]+"\n")
out=open(f"{WORK}/refqual_degradation.tsv","w")
out.write("completeness\tcontamination\tn_contigs\tn_clusters\tprecision\trecall\tbray_curtis\n")
print(f"{len(truth_acc)} truth strains degraded; {len(background)} background kept complete")
for (comp,contam,nc) in levels:
    tag=int(comp*100); pdir=f"{WORK}/panel_{tag}"; os.makedirs(pdir,exist_ok=True)
    for f in glob.glob(f"{pdir}/*.fna"): os.remove(f)
    for acc in truth_acc:
        if comp>=1.0: shutil.copy(f"{GEN}/{acc}.fna",f"{pdir}/{acc}.fna")
        else:
            seq=read_seq(f"{GEN}/{acc}.fna")
            write_fasta(degrade(seq,comp,contam,nc,contaminant,sum(ord(c) for c in acc)+tag),f"{pdir}/{acc}.fna")
    for g in background: os.symlink(os.path.abspath(f"{GEN}/{g}.fna"),f"{pdir}/{g}.fna")
    db=f"{WORK}/db_{tag}.tsv"
    r=subprocess.run([BIN,"cluster","--genomes",pdir,"--enzyme",SET,"--out",db,"--similarity","0.95"],capture_output=True,text=True,env=env)
    m=re.search(r'into (\d+) cluster',r.stdout); ncl=int(m.group(1)) if m else 0
    g2c={}
    for l in open(db.replace(".tsv",".members.tsv")):
        if not l.startswith("#"): g,c=l.strip().split("\t"); g2c[g]=c
    TP=FP=FN=0; bcs=[]
    for i in range(1,6):
        pred=f"{WORK}/s{i}_{tag}.pred"
        subprocess.run([BIN,"profile","--db",db,"--reads",f"{READS}/sample{i}.fq","--out",pred],capture_output=True,text=True,env=env)
        ab={}
        for g,a in sample_truth[i]: c=g2c.get(g,g); ab[c]=ab.get(c,0)+a
        tf=f"{WORK}/s{i}_{tag}.truth"
        open(tf,"w").write("".join(f"{c}\t{a}\n" for c,a in ab.items()))
        ev=subprocess.run([BIN,"evaluate","--pred",pred,"--truth",tf,"--present","0.01"],capture_output=True,text=True,env=env).stdout
        TP+=int(re.search(r'TP=(\d+)',ev).group(1)); FP+=int(re.search(r'FP=(\d+)',ev).group(1)); FN+=int(re.search(r'FN=(\d+)',ev).group(1))
        bcs.append(float(re.search(r'Bray-Curtis=([0-9.]+)',ev).group(1)))
    P=TP/(TP+FP) if TP+FP else 0; R=TP/(TP+FN) if TP+FN else 0; BC=sum(bcs)/len(bcs)
    out.write(f"{comp:.2f}\t{contam:.2f}\t{nc}\t{ncl}\t{P:.3f}\t{R:.3f}\t{BC:.3f}\n")
    print(f"completeness={comp:.2f} contam={contam:.2f} contigs={nc}: clusters={ncl} P={P:.2f} R={R:.2f} BC={BC:.3f}",flush=True)
out.close(); print("done -> refqual/refqual_degradation.tsv")
