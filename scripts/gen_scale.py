import random, os
random.seed(3)
base="scaling"; os.makedirs(base, exist_ok=True)
GEN=os.path.join(base,"genomes"); os.makedirs(GEN, exist_ok=True)
N=400; GROUP=5  # N genomes in N/GROUP clusters of ~5 near-identical strains
def rseq(n): return ''.join(random.choice("ACGT") for _ in range(n))
g=0
for grp in range(N//GROUP):
    core=rseq(200000)  # shared within group (cluster)
    for s in range(GROUP):
        seq=core+rseq(30000)  # + strain-private
        with open(os.path.join(GEN,f"g{g:04d}.fna"),"w") as f:
            f.write(f">g{g:04d}\n")
            for i in range(0,len(seq),70): f.write(seq[i:i+70]+"\n")
        g+=1
print(f"generated {g} genomes in {N//GROUP} groups -> {GEN}")
# subset dirs via symlink
import glob
allf=sorted(glob.glob(os.path.join(GEN,"*.fna")))
for K in [50,100,200,400]:
    d=os.path.join(base,f"g{K}"); os.makedirs(d, exist_ok=True)
    for f in glob.glob(os.path.join(d,"*.fna")): os.remove(f)
    for f in allf[:K]: os.symlink(os.path.abspath(f), os.path.join(d, os.path.basename(f)))
    print(f"  subset g{K}: {K} genomes")
