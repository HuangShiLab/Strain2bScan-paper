#!/usr/bin/env python3
"""Profile the 32 diurnal saliva 2bRAD samples (4 subjects x 8 timepoints) against the
62-species BcgI DB. Emit a long-format strain-abundance table for the individual-discrimination
analysis (Fig 10)."""
import subprocess, os, glob, time, gzip, sys

PAPER = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.environ.get("S2BS_SALIVA_WORK", f"{PAPER}/work/saliva_rerun")
BIN = os.environ.get("S2BS_BIN", f"{os.path.dirname(PAPER)}/Strain2bScan/target/release/strain2bscan")
DBS  = os.environ.get("DBS", f"{WORK}/dbs/bcgi")
GATE = os.environ.get("GATE", "50")
DETECT = os.environ.get("DETECT", "5")
READS = sorted(glob.glob(f"{WORK}/reads/saliva_2brad/*_1.fastq.gz"))
OUTDIR = f"{WORK}/results/saliva_preds"; os.makedirs(OUTDIR, exist_ok=True)

def profile(reads):
    if subprocess.run(["gzip", "-t", reads]).returncode != 0:
        return None, "corrupt"
    tmp = f"{WORK}/results/_tmp_saliva.fq"
    with gzip.open(reads, "rb") as fi, open(tmp, "wb") as fo:
        while True:
            b = fi.read(1 << 22)
            if not b: break
            fo.write(b)
    cmd = [BIN, "multi-profile", "--dbs", DBS, "--reads", tmp, "--enzyme", "BcgI",
           "--min-species-markers", GATE, "--min-species-detect", DETECT]
    env = {**os.environ, "STRAIN2BSCAN_THREADS": "16"}
    out = subprocess.run(cmd, capture_output=True, text=True, env=env)
    os.remove(tmp)
    return out.stdout, out.stderr

# Columns are looked up by NAME from multi-profile's own header line, never by position.
# The previous positional parse stopped at index 4 and so silently dropped depth,
# global_abundance and sample_fraction -- which left the downstream community matrices with
# no cross-species abundance to use, and they fell back to `support` (a marker count).
WANT = ["cluster", "abundance", "coverage", "support", "depth",
        "global_abundance", "sample_fraction"]

long_rows = []   # sample, subject, timepoint, species, then WANT in order
detected_rows = []
print(f"{'sample':10s} {'markers':>8s} {'resolved':>8s} {'strains':>7s} {'detect':>6s} {'t_s':>5s}")
for r in READS:
    base = os.path.basename(r)                       # S11-1__SRR29712951_1.fastq.gz
    alias = base.split("__")[0]                       # S11-1
    # alias = S<timepoint>-<subject>: prefix S9/S11/S13/S17 = time of day (9AM/11AM/1PM/5PM),
    # suffix 1..14 = the 8 subjects. (Verified against PERMANOVA: subject R^2~0.83, timepoint ns.)
    tp, subj = alias.split("-")
    t0 = time.time()
    so, se = profile(r)
    dt = time.time() - t0
    if so is None:
        print(f"{alias:10s}  {se}")
        continue
    open(f"{OUTDIR}/{alias}.txt", "w").write(so)
    nmark = "?"; nres = 0; nstrain = 0; ndet = 0
    hdr = None
    for line in so.splitlines():
        if line.startswith("sample:"):
            nmark = line.split()[1]
        if line.startswith("#species\t"):
            hdr = {n: i for i, n in enumerate(line.lstrip("#").split("\t"))}
            missing = [c for c in WANT if c not in hdr]
            if missing:
                sys.exit(f"{alias}: multi-profile output has no {missing} column(s); "
                         f"header was {sorted(hdr)}. Refusing to write a table that would "
                         f"silently lose them.")
            continue
        if not line.startswith("  "):
            continue
        parts = line.strip().split("\t")
        sp = parts[0]
        if "[detected, not strain-resolvable]" in line:
            ndet += 1
            detected_rows.append((alias, subj, tp, sp))
        elif "[strain-resolved, no cluster" in line:
            nres += 1
        elif hdr is not None and len(parts) == len(hdr):
            long_rows.append((alias, subj, tp, sp, *(parts[hdr[c]] for c in WANT)))
            nstrain += 1
    if hdr is None:
        sys.exit(f"{alias}: no '#species' header in multi-profile output -- the output format "
                 f"changed. Fix the parser rather than letting rows be dropped.")
    # count distinct resolved species
    nres = len(set(row[3] for row in long_rows if row[0] == alias))
    print(f"{alias:10s} {nmark:>8s} {nres:>8d} {nstrain:>7d} {ndet:>6d} {dt:>5.1f}")

with open(f"{WORK}/results/saliva_strain_long.tsv", "w") as w:
    # Append-only: the first eight columns keep their historical positions and names, so a
    # table written by an older run stays readable; the three that were being dropped go on
    # the end. `sample_fraction` is the one to use for a cross-species community matrix.
    w.write("sample\tsubject\ttimepoint\tspecies\tcluster\twithin_abund\tcoverage\tsupport"
            "\tdepth\tglobal_abundance\tsample_fraction\n")
    for r in long_rows:
        w.write("\t".join(str(x) for x in r) + "\n")
with open(f"{WORK}/results/saliva_detected_long.tsv", "w") as w:
    w.write("sample\tsubject\ttimepoint\tspecies\n")
    for r in detected_rows:
        w.write("\t".join(str(x) for x in r) + "\n")
print(f"\nwrote saliva_strain_long.tsv ({len(long_rows)} strain calls across {len(set(r[0] for r in long_rows))} samples)")
