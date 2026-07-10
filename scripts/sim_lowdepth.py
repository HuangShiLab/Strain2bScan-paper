#!/usr/bin/env python3
"""Low-depth variant of sim_samples.py: species depths span sub-1x (the low-abundance regime),
so the calibration can see how recall and the three-tier split behave when species are faint.
Writes multispecies/samples_low/*.fq + multispecies/truth_low/*.truth.tsv. Args: NSAMP SPP."""
import random, glob, os, sys
random.seed(23)
GEN = "multispecies/genomes"
OUT = "multispecies"
os.makedirs(f"{OUT}/samples_low", exist_ok=True)
os.makedirs(f"{OUT}/truth_low", exist_ok=True)
NSAMP = int(sys.argv[1]) if len(sys.argv) > 1 else 20
SPP = int(sys.argv[2]) if len(sys.argv) > 2 else 12
comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
species = [d for d in sorted(os.listdir(GEN)) if os.path.isdir(f"{GEN}/{d}") and glob.glob(f"{GEN}/{d}/*.fna")]


def load(fa):
    return "".join(l.strip() for l in open(fa) if not l.startswith(">")).upper()


def sim(seq, depth, fh, tag):
    L = len(seq)
    n = int(depth * L / 150)
    for i in range(n):
        p = random.randint(0, max(0, L - 150))
        r = seq[p:p + 150]
        if random.random() < 0.5:
            r = "".join(comp.get(c, "N") for c in reversed(r))
        fh.write(f"@{tag}_{i}\n{r}\n+\n{'I' * 150}\n")
    return n


for s in range(NSAMP):
    chosen = random.sample(species, min(SPP, len(species)))
    truth = []
    fq = open(f"{OUT}/samples_low/sample_{s:02d}.fq", "w")
    total = 0
    for sp in chosen:
        strains = glob.glob(f"{GEN}/{sp}/*.fna")
        k = 2 if (len(strains) >= 2 and random.random() < 0.4) else 1
        for fa in random.sample(strains, k):
            acc = os.path.basename(fa)[:-4]
            # low-abundance: median ~0.67x, spread ~0.2-2.5x, floored at 0.1x
            depth = round(random.lognormvariate(-0.4, 0.7), 2)
            depth = max(depth, 0.1)
            n = sim(load(fa), depth, fq, acc)
            total += n
            truth.append((sp, acc, depth))
    fq.close()
    with open(f"{OUT}/truth_low/sample_{s:02d}.truth.tsv", "w") as t:
        t.write("#species\tstrain\tdepth\n")
        for sp, acc, d in truth:
            t.write(f"{sp}\t{acc}\t{d}\n")
    print(f"sample_{s:02d}: {len(chosen)} species, {len(truth)} strains, {total} reads", flush=True)
print(f"generated {NSAMP} low-depth samples")
