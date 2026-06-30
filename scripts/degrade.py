import random, sys
def read_seq(path):
    s=[]
    for l in open(path):
        if not l.startswith(">"): s.append(l.strip())
    return "".join(s).upper()
def split_into(seq, n):
    n=max(1,n); k=max(1,len(seq)//n)
    return [seq[i:i+k] for i in range(0,len(seq),k)][:n] or [seq]
def degrade(seq, completeness, contam_frac, n_contigs, contaminant, seed, win=5000):
    rng=random.Random(seed)
    wins=[seq[i:i+win] for i in range(0,len(seq),win)]
    idx=list(range(len(wins))); rng.shuffle(idx)
    target=completeness*len(seq); cum=0; keep=set()
    for i in idx:
        if cum>=target: break
        keep.add(i); cum+=len(wins[i])
    kept="".join(wins[i] for i in sorted(keep))
    frags=split_into(kept, n_contigs)
    clen=int(contam_frac*len(kept))
    if clen>0 and contaminant:
        st=rng.randint(0,max(0,len(contaminant)-clen))
        frags+=split_into(contaminant[st:st+clen], max(1,n_contigs//10))
    return frags
if __name__=="__main__":
    inp,out,comp,contam,nc,cpath,seed=sys.argv[1:8]
    seq=read_seq(inp); cont=read_seq(cpath) if cpath!="-" else ""
    frags=degrade(seq,float(comp),float(contam),int(nc),cont,int(seed))
    with open(out,"w") as f:
        for i,fr in enumerate(frags):
            f.write(f">ctg{i}\n")
            for j in range(0,len(fr),70): f.write(fr[j:j+70]+"\n")
    tot=sum(len(x) for x in frags)
    print(f"  {out}: {len(frags)} contigs, {tot} bp (orig {len(seq)}, kept~{comp}, contam~{contam})")
