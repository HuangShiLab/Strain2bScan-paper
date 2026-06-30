import subprocess, glob, os, collections, sys
BIN=os.environ.get("STRAIN2BSCAN_BIN","strain2bscan")
SET="CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"
DBS="multispecies/dbs"; GATE=sys.argv[1] if len(sys.argv)>1 else "30"
def members(sp):
    m={}; f=f"{DBS}/{sp}.members.tsv"
    if os.path.exists(f):
        for line in open(f):
            if not line.startswith('#'): g,c=line.strip().split('\t'); m[g]=c
    return m
agg=collections.Counter()
for s in sorted(glob.glob("multispecies/samples/*.fq")):
    name=os.path.basename(s)[:-3]; truth=[]
    for line in open(f"multispecies/truth/{name}.truth.tsv"):
        if not line.startswith('#'): sp,acc,d=line.strip().split('\t'); truth.append((sp,acc,float(d)))
    truth_sp=set(t[0] for t in truth)
    env={**os.environ,"STRAIN2BSCAN_THREADS":"16"}
    out=subprocess.run([BIN,"multi-profile","--dbs",DBS,"--reads",s,"--enzyme",SET,"--min-species-markers",GATE],
                       capture_output=True,text=True,env=env).stdout
    det=collections.defaultdict(set)
    for line in out.splitlines():
        if line.startswith('  ') and '\t' in line:
            p=line.strip().split('\t')
            if len(p)>=2: det[p[0]].add(p[1])
    det_sp=set(det.keys())
    agg['tp_sp']+=len(truth_sp&det_sp); agg['fp_sp']+=len(det_sp-truth_sp); agg['fn_sp']+=len(truth_sp-det_sp)
    for sp,acc,depth in truth:
        cl=members(sp).get(acc); hit=(cl in det.get(sp,set()))
        agg['tp_st' if hit else 'fn_st']+=1
        if depth>=1.0: agg['tp_hi' if hit else 'fn_hi']+=1
sP=agg['tp_sp']/(agg['tp_sp']+agg['fp_sp']); sR=agg['tp_sp']/(agg['tp_sp']+agg['fn_sp'])
print(f"gate={GATE:>4}: species P={sP:.3f} R={sR:.3f} | strain recall={agg['tp_st']/(agg['tp_st']+agg['fn_st']):.3f} (>=1x {agg['tp_hi']/(agg['tp_hi']+agg['fn_hi']):.3f}) | FP_sp={agg['fp_sp']}")
