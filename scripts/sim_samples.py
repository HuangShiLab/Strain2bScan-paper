import random, glob, os, sys
random.seed(17)
GEN="multispecies/genomes"; OUT="multispecies"; 
os.makedirs(f"{OUT}/samples",exist_ok=True); os.makedirs(f"{OUT}/truth",exist_ok=True)
NSAMP=int(sys.argv[1]) if len(sys.argv)>1 else 20
SPP=int(sys.argv[2]) if len(sys.argv)>2 else 8       # species per sample
comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}
species=[d for d in sorted(os.listdir(GEN)) if os.path.isdir(f"{GEN}/{d}") and glob.glob(f"{GEN}/{d}/*.fna")]
def load(fa):
    return "".join(l.strip() for l in open(fa) if not l.startswith(">")).upper()
def sim(seq,depth,fh,tag):
    L=len(seq); n=int(depth*L/150)
    for i in range(n):
        p=random.randint(0,max(0,L-150)); r=seq[p:p+150]
        if random.random()<0.5: r="".join(comp.get(c,"N") for c in reversed(r))
        fh.write(f"@{tag}_{i}\n{r}\n+\n{'I'*150}\n")
    return n
for s in range(NSAMP):
    chosen=random.sample(species, min(SPP,len(species)))
    truth=[]
    fq=open(f"{OUT}/samples/sample_{s:02d}.fq","w"); total=0
    for sp in chosen:
        strains=glob.glob(f"{GEN}/{sp}/*.fna")
        k=2 if (len(strains)>=2 and random.random()<0.4) else 1   # ~40% species get 2 strains
        for fa in random.sample(strains,k):
            acc=os.path.basename(fa)[:-4]
            depth=round(random.lognormvariate(0.7,0.5),2)  # median ~2x, all >~0.8x
            depth=max(depth,1.0)
            n=sim(load(fa),depth,fq,acc); total+=n
            truth.append((sp,acc,depth))
    fq.close()
    with open(f"{OUT}/truth/sample_{s:02d}.truth.tsv","w") as t:
        t.write("#species\tstrain\tdepth\n")
        for sp,acc,d in truth: t.write(f"{sp}\t{acc}\t{d}\n")
    print(f"sample_{s:02d}: {len(chosen)} species, {len(truth)} strains, {total} reads",flush=True)
print(f"generated {NSAMP} samples")
