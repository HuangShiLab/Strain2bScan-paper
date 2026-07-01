#!/usr/bin/env python3
"""Same-panel head-to-head: StrainScan built on the IDENTICAL C. acnes 64-genome panel, profiled
on the same 5 mock samples, for a fair accuracy + resource comparison against Strain2bScan.

LINUX ONLY: StrainScan_build needs `dashing` (no osx-arm64 build). Run from the workdir after
`make refqual-data` (provides acnes/{genomes,reads,truth}).

Prereqs:
  conda env create -f env/strainscan-linux.yml && conda activate strainscan
  export STRAINSCAN=/path/to/StrainScan            # cloned repo (has StrainScan.py, StrainScan_build.py)
  export STRAINSCAN_PY=$(which python)             # the strainscan env python

Outputs: refqual/headtohead_strainscan.tsv  (per-sample StrainScan time/RSS + detection vs truth)
and appends a comparison note. Strain2bScan side is in results/headtohead_performance.tsv.
"""
import os, re, glob, subprocess, statistics

SS = os.environ["STRAINSCAN"]; PY = os.environ.get("STRAINSCAN_PY", "python")
GEN, READS, TRUTH, WORK = "acnes/genomes", "acnes/reads", "acnes/truth", "refqual"
os.makedirs(WORK, exist_ok=True)
DB = f"{WORK}/ss_db_panel"


def timed(cmd):
    r = subprocess.run(["/usr/bin/time", "-v"] + cmd, capture_output=True, text=True)
    t = re.search(r"Elapsed .*?: ([0-9:.]+)", r.stderr)
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", r.stderr)
    secs = 0.0
    if t:
        p = t.group(1).split(":"); secs = float(p[-1]) + (float(p[-2]) * 60 if len(p) > 1 else 0)
    return r, secs, (int(m.group(1)) / 1024 if m else 0)


# 1) build StrainScan DB on the identical 64-genome panel (needs dashing)
if not os.path.isdir(DB):
    print("building StrainScan DB on the 64-genome panel ...", flush=True)
    subprocess.run([PY, f"{SS}/StrainScan_build.py", "-i", GEN, "-o", DB], check=True)

# 2) profile each sample; record time, RSS, and detected strains
out = open(f"{WORK}/headtohead_strainscan.tsv", "w")
out.write("sample\tstrainscan_time_s\tstrainscan_maxRSS_MB\tdetected_strains\n")
ts, ms = [], []
for i in range(1, 6):
    od = f"{WORK}/ss_out_s{i}"
    r, t, mrss = timed([PY, f"{SS}/StrainScan.py", "-i", f"{READS}/sample{i}.fq", "-d", DB, "-o", od])
    ts.append(t); ms.append(mrss)
    rep = f"{od}/final_report.txt"
    strains = []
    if os.path.exists(rep):
        for l in open(rep).read().splitlines()[1:]:
            f = l.split("\t")
            if len(f) > 1 and f[1]:
                strains.append(f[1])
    out.write(f"sample{i}\t{t:.2f}\t{mrss:.0f}\t{';'.join(strains)}\n")
    print(f"sample{i}: {t:.2f}s {mrss:.0f}MB -> {strains}", flush=True)
out.write(f"mean\t{statistics.mean(ts):.2f}\t{statistics.mean(ms):.0f}\t-\n")
out.close()
print("done -> refqual/headtohead_strainscan.tsv")
print("NOTE: compare against Strain2bScan in results/headtohead_performance.tsv (same 5 samples,")
print("      same 64-genome panel). Also map StrainScan strain names to truth clusters for accuracy.")
