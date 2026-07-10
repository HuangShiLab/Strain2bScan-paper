#!/usr/bin/env python3
"""Profile the MSA-1002 mock host-contamination series (2bRAD + shotgun WMS) with Strain2bScan
against the 62-species BcgI DB, and score species precision/recall vs the 20-species truth.
Shotgun WMS is digested in-silico with BcgI -> like-for-like vs native 2bRAD."""
import subprocess, os, sys, time, gzip

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
BIN  = "/Users/macstudio/Downloads/Strain2bScan/target/release/strain2bscan"
DBS  = f"{WORK}/dbs/bcgi"
TRUTH = set(l.strip() for l in open(f"{WORK}/truth_msa1002.txt") if l.strip())
GATE = os.environ.get("GATE", "50")          # --min-species-markers (Fig9 used 50)
DETECT = os.environ.get("DETECT", "5")        # --min-species-detect

# sample_alias -> (data_type, host%, rep, R1 path)
SAMPLES = [
    ("2bRAD", 90, 1, "mock_2brad/SRR29727660_1.fastq.gz"),
    ("2bRAD", 90, 2, "mock_2brad/SRR29727659_1.fastq.gz"),
    ("2bRAD", 99, 1, "mock_2brad/SRR29727658_1.fastq.gz"),
    ("2bRAD", 99, 2, "mock_2brad/SRR29727657_1.fastq.gz"),
    ("WMS",    0, 1, "mock_wms/SRR29727654_1.fastq.gz"),   # control (no host spike)
    ("WMS",   90, 1, "mock_wms/SRR29727656_1.fastq.gz"),
    ("WMS",   99, 1, "mock_wms/SRR29727655_1.fastq.gz"),
]

def profile(reads):
    # Strain2bScan reads plain FASTQ only -> decompress .gz to a scratch file first.
    tmp = None
    if reads.endswith(".gz"):
        if subprocess.run(["gzip", "-t", reads]).returncode != 0:
            raise RuntimeError(f"corrupt gz: {reads}")
        tmp = f"{WORK}/results/_tmp_reads.fq"
        with gzip.open(reads, "rb") as fi, open(tmp, "wb") as fo:
            while True:
                b = fi.read(1 << 22)
                if not b: break
                fo.write(b)
        reads = tmp
    cmd = [BIN, "multi-profile", "--dbs", DBS, "--reads", reads, "--enzyme", "BcgI",
           "--min-species-markers", GATE, "--min-species-detect", DETECT]
    env = {**os.environ, "STRAIN2BSCAN_THREADS": "16"}
    t0 = time.time()
    out = subprocess.run(cmd, capture_output=True, text=True, env=env)
    dt = time.time() - t0
    if tmp and os.path.exists(tmp):
        os.remove(tmp)
    resolved, detected = set(), set()
    for line in out.stdout.splitlines():
        if not line.startswith("  "):
            continue
        sp = line.strip().split("\t")[0]
        if "[detected, not strain-resolvable]" in line:
            detected.add(sp)
        else:
            resolved.add(sp)  # strain-resolved (cluster line or "no cluster above threshold")
    return resolved, detected, dt, out.stdout, out.stderr

rows = []
os.makedirs(f"{WORK}/results/mock_preds", exist_ok=True)
print(f"{'sample':22s} {'reads_or_err':>12s} {'TP':>3s} {'FP':>3s} {'prec':>5s} {'recall':>6s} {'det':>4s} {'t_s':>6s}")
for dt_type, host, rep, rel in SAMPLES:
    reads = f"{WORK}/reads/{rel}"
    name = f"{dt_type}_host{host}_rep{rep}"
    if not os.path.exists(reads) or os.path.getsize(reads) == 0:
        print(f"{name:22s} {'MISSING':>12s}")
        continue
    resolved, detected, dt, so, se = profile(reads)
    open(f"{WORK}/results/mock_preds/{name}.txt", "w").write(so)
    tp = len(resolved & TRUTH); fp = len(resolved - TRUTH)
    prec = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / len(TRUTH)
    # pull sample tag-marker count from the header line
    nmark = "?"
    for l in so.splitlines():
        if l.startswith("sample:"):
            nmark = l.split()[1]; break
    print(f"{name:22s} {nmark:>12s} {tp:>3d} {fp:>3d} {prec:>5.2f} {recall:>6.2f} {len(detected):>4d} {dt:>6.1f}")
    rows.append((dt_type, host, rep, nmark, tp, fp, f"{prec:.4f}", f"{recall:.4f}", len(detected), f"{dt:.1f}",
                 ",".join(sorted(resolved - TRUTH)) or "-"))

with open(f"{WORK}/results/mock_hostcontam.tsv", "w") as w:
    w.write("data_type\thost_pct\trep\tsample_markers\tTP\tFP\tprecision\trecall\tn_detected_not_resolved\ttime_s\tfalse_positives\n")
    for r in rows:
        w.write("\t".join(str(x) for x in r) + "\n")
print(f"\nwrote {WORK}/results/mock_hostcontam.tsv  ({len(rows)} samples)")
