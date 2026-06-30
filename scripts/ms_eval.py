import subprocess, glob, os, collections
BIN=os.environ.get("STRAIN2BSCAN_BIN","strain2bscan")
SET="CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"
DBS="multispecies/dbs"
def members(species):
    m={}; f=f"{DBS}/{species}.members.tsv"
    if os.path.exists(f):
        for line in open(f):
            if not line.startswith('#'):
                g,c=line.strip().split('\t'); m[g]=c
    return m
samples=sorted(glob.glob("multispecies/samples/*.fq"))
agg=collections.Counter()
print(f"{'sample':<12} {'sp_P':>5} {'sp_R':>5} {'strainR':>8} {'strainR>=1x':>11}")
for s in samples:
    name=os.path.basename(s)[:-3]
    truth=[]
    for line in open(f"multispecies/truth/{name}.truth.tsv"):
        if not line.startswith('#'):
            sp,acc,depth=line.strip().split('\t'); truth.append((sp,acc,float(depth)))
    truth_sp=set(t[0] for t in truth)
    env={**os.environ,"STRAIN2BSCAN_THREADS":"16"}
    out=subprocess.run([BIN,"multi-profile","--dbs",DBS,"--reads",s,"--enzyme",SET],
                       capture_output=True,text=True,env=env).stdout
    det=collections.defaultdict(set)
    for line in out.splitlines():
        if line.startswith('  ') and '\t' in line:
            p=line.strip().split('\t')
            if len(p)>=2: det[p[0]].add(p[1])
    det_sp=set(det.keys())
    tp_sp=len(truth_sp&det_sp); fp_sp=len(det_sp-truth_sp); fn_sp=len(truth_sp-det_sp)
    # strain-level recall
    tp_st=fn_st=0; tp_hi=fn_hi=0
    for sp,acc,depth in truth:
        cl=members(sp).get(acc); hit=(cl in det.get(sp,set()))
        if hit: tp_st+=1
        else: fn_st+=1
        if depth>=1.0:
            if hit: tp_hi+=1
            else: fn_hi+=1
    for k,v in [('tp_sp',tp_sp),('fp_sp',fp_sp),('fn_sp',fn_sp),('tp_st',tp_st),('fn_st',fn_st),('tp_hi',tp_hi),('fn_hi',fn_hi)]: agg[k]+=v
    sP=tp_sp/(tp_sp+fp_sp) if tp_sp+fp_sp else 0; sR=tp_sp/(tp_sp+fn_sp) if tp_sp+fn_sp else 0
    stR=tp_st/(tp_st+fn_st) if tp_st+fn_st else 0; hiR=tp_hi/(tp_hi+fn_hi) if tp_hi+fn_hi else 0
    print(f"{name:<12} {sP:>5.2f} {sR:>5.2f} {stR:>8.2f} {hiR:>11.2f}")
sP=agg['tp_sp']/(agg['tp_sp']+agg['fp_sp']); sR=agg['tp_sp']/(agg['tp_sp']+agg['fn_sp'])
print(f"\nTOTAL  species precision={sP:.3f} recall={sR:.3f} | strain recall={agg['tp_st']/(agg['tp_st']+agg['fn_st']):.3f} | strain recall(>=1x)={agg['tp_hi']/(agg['tp_hi']+agg['fn_hi']):.3f}")
