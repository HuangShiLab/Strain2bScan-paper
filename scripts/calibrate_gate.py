#!/usr/bin/env python3
"""Calibrate the Layer-1 breadth gate on the 55-species community.

Sweeps --min-species-marker-frac (at fixed floor + detect) and reports, aggregated over all
simulated samples: species-level precision/recall where a species counts as a POSITIVE call only
if it is strain-RESOLVED (tier Resolved) -- that is exactly what the gate controls -- plus
strain-level recall (and recall restricted to >=1x truth strains) to expose the accuracy cost of
raising the bar. Recommends the smallest frac that keeps species precision >= TARGET.

Run from the workdir that contains multispecies/{dbs,samples,truth}. Env: FLOOR, DETECT, TARGET.
"""
import subprocess, glob, os, collections

BIN = os.environ.get("STRAIN2BSCAN_BIN", os.path.expanduser("~/Downloads/Strain2bScan/target/release/strain2bscan"))
SET = "CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"
DBS = os.environ.get("DBS_DIR", "multispecies/dbs")
SAMPLES_DIR = os.environ.get("SAMPLES_DIR", "multispecies/samples")
TRUTH_DIR = os.environ.get("TRUTH_DIR", "multispecies/truth")
SAMPLES = sorted(glob.glob(f"{SAMPLES_DIR}/*.fq"))
FLOOR = int(os.environ.get("FLOOR", "200"))
DETECT = int(os.environ.get("DETECT", "10"))
TARGET = float(os.environ.get("TARGET", "0.98"))
FRACS = [float(x) for x in os.environ.get("FRACS", "0.0,0.02,0.04,0.06,0.08,0.10,0.15,0.20,0.30").split(",")]


def members(species):
    m = {}
    f = f"{DBS}/{species}.members.tsv"
    if os.path.exists(f):
        for line in open(f):
            if not line.startswith('#'):
                p = line.strip().split('\t')
                if len(p) >= 2:
                    m[p[0]] = p[1]
    return m


truth_by_sample = {}
for s in SAMPLES:
    name = os.path.basename(s)[:-3]
    tr = []
    for line in open(f"{TRUTH_DIR}/{name}.truth.tsv"):
        if not line.startswith('#'):
            sp, acc, depth = line.strip().split('\t')
            tr.append((sp, acc, float(depth)))
    truth_by_sample[s] = tr


def run(sample, frac):
    env = {**os.environ, "STRAIN2BSCAN_THREADS": "16"}
    out = subprocess.run(
        [BIN, "multi-profile", "--dbs", DBS, "--reads", sample, "--enzyme", SET,
         "--min-species-markers", str(FLOOR), "--min-species-detect", str(DETECT),
         "--min-species-marker-frac", str(frac)],
        capture_output=True, text=True, env=env).stdout
    resolved = collections.defaultdict(set)   # sp -> clusters called
    resolved_sp, detected_sp = set(), set()
    for line in out.splitlines():
        if not (line.startswith('  ') and '\t' in line):
            continue
        p = line.strip().split('\t')
        if len(p) < 2:
            continue
        sp, f2 = p[0], p[1]
        if f2.startswith('[detected'):
            detected_sp.add(sp)
        elif f2.startswith('[strain-resolved'):
            resolved_sp.add(sp)
        else:
            resolved_sp.add(sp)
            resolved[sp].add(f2)
    return resolved_sp, detected_sp, resolved


print(f"# calibration  FLOOR={FLOOR}  DETECT={DETECT}  samples={len(SAMPLES)}  target_precision={TARGET}")
print("frac\tsp_precision\tsp_recall\tstrain_recall\tstrain_recall_ge1x\tFP_species\tdet_not_resolved\tresolved_calls")
rows = []
for frac in FRACS:
    agg = collections.Counter()
    for s in SAMPLES:
        truth = truth_by_sample[s]
        truth_sp = set(t[0] for t in truth)
        resolved_sp, detected_sp, resolved = run(s, frac)
        agg['tp'] += len(resolved_sp & truth_sp)
        agg['fp'] += len(resolved_sp - truth_sp)
        agg['fn'] += len(truth_sp - resolved_sp)
        agg['dnr'] += len(detected_sp)
        agg['ncalls'] += sum(len(v) for v in resolved.values())
        for sp, acc, depth in truth:
            cl = members(sp).get(acc)
            hit = cl in resolved.get(sp, set())
            agg['st_tp'] += 1 if hit else 0
            agg['st_n'] += 1
            if depth >= 1.0:
                agg['hi_tp'] += 1 if hit else 0
                agg['hi_n'] += 1
    P = agg['tp'] / (agg['tp'] + agg['fp']) if (agg['tp'] + agg['fp']) else 0.0
    R = agg['tp'] / (agg['tp'] + agg['fn']) if (agg['tp'] + agg['fn']) else 0.0
    stR = agg['st_tp'] / agg['st_n'] if agg['st_n'] else 0.0
    hiR = agg['hi_tp'] / agg['hi_n'] if agg['hi_n'] else 0.0
    print(f"{frac:.2f}\t{P:.3f}\t{R:.3f}\t{stR:.3f}\t{hiR:.3f}\t{agg['fp']}\t{agg['dnr']}\t{agg['ncalls']}", flush=True)
    rows.append((frac, P, R, stR, hiR))

cand = [r for r in rows if r[1] >= TARGET]
if cand:
    best = min(cand, key=lambda r: r[0])
    print(f"\n# RECOMMEND --min-species-marker-frac {best[0]:.2f}  "
          f"(species precision {best[1]:.3f}, recall {best[2]:.3f}, strain_recall {best[3]:.3f}, >=1x {best[4]:.3f})")
else:
    print(f"\n# no frac reached precision {TARGET}; best = {max(r[1] for r in rows):.3f}")
