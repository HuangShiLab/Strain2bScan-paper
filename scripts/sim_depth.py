import random, os
random.seed(11)
G="acnes/genomes/GCF_009737125.1.fna"
seq=[]
for line in open(G):
    if not line.startswith(">"): seq.append(line.strip())
seq="".join(seq).upper()
L=len(seq); comp={"A":"T","T":"A","C":"G","G":"C","N":"N"}
od="depthsweep"; os.makedirs(od, exist_ok=True)
for X in [0.1,0.5,1.0,5.0]:
    n=int(round(X*L/150))
    fn=f"{od}/cae_{X}X.fq"
    with open(fn,"w") as f:
        for i in range(n):
            p=random.randint(0,L-150); r=seq[p:p+150]
            if random.random()<0.5: r="".join(comp.get(c,"N") for c in reversed(r))
            f.write(f"@r{i}\n{r}\n+\n{'I'*150}\n")
    print(f"{X}X -> {n} reads -> {fn}")
print("genome length:",L)
