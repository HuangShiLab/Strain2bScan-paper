#!/usr/bin/env python3
"""Fig 10 validation: paired shotgun (WMS, in-silico BcgI) vs native 2bRAD on the SAME saliva
samples. Saliva WMS is host-dominated, so shotgun recovers far fewer bacterial strain markers than
native 2bRAD (reinforcing the enrichment advantage). The concordance question is directional: are
the strains WMS *does* call confirmed by native 2bRAD on the SAME sample (and more so than on other
samples)? Metric: fraction of WMS strain/species calls present in same-sample vs cross-sample 2bRAD."""
import subprocess, os, glob, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
BIN  = "/Users/macstudio/Downloads/Strain2bScan/target/release/strain2bscan"
DBS  = f"{WORK}/dbs/oral"
MAX_READS = 400_000_000
os.makedirs(f"{WORK}/results/wms_preds", exist_ok=True)

def parse_calls(stdout):
    calls = {}
    for line in stdout.splitlines():
        if not line.startswith("  "): continue
        p = line.strip().split("\t")
        if "[detected" in line or "[strain-resolved, no cluster" in line: continue
        if len(p) >= 5:
            calls[f"{p[0]}|{p[1]}"] = calls.get(f"{p[0]}|{p[1]}", 0.0) + float(p[4])
    return calls

def profile_wms(reads):
    alias = os.path.basename(reads).split("__")[0]
    cache = f"{WORK}/results/wms_preds/{alias}.txt"
    if os.path.exists(cache) and os.path.getsize(cache) > 0:
        return parse_calls(open(cache).read())
    tmp = f"{WORK}/results/_tmp_wms.fq"
    with open(tmp, "wb") as fo:
        p1 = subprocess.Popen(["gzip", "-dc", reads], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(["head", "-n", str(MAX_READS * 4)], stdin=p1.stdout, stdout=fo)
        p1.stdout.close(); p2.communicate(); p1.wait()
    if os.path.getsize(tmp) == 0:
        os.remove(tmp); return {}
    out = subprocess.run([BIN, "multi-profile", "--dbs", DBS, "--reads", tmp, "--enzyme", "BcgI",
                          "--min-species-markers", "50", "--min-species-detect", "5"],
                         capture_output=True, text=True, env={**os.environ, "STRAIN2BSCAN_THREADS": "16"})
    os.remove(tmp)
    open(cache, "w").write(out.stdout)
    return parse_calls(out.stdout)

# native 2bRAD strain calls per sample
rows = [l.rstrip("\n").split("\t") for l in open(f"{WORK}/results/saliva_strain_long.tsv")][1:]
brad = {}
for r in rows:
    brad.setdefault(r[0], {})[f"{r[3]}|{r[4]}"] = 1.0
brad_sp = {s: set(k.split("|")[0] for k in d) for s, d in brad.items()}

wms = {}
for f in sorted(glob.glob(f"{WORK}/reads/saliva_wms/*_1.fastq.gz")):
    alias = os.path.basename(f).split("__")[0]
    c = profile_wms(f)
    wms[alias] = c
    print(f"WMS {alias}: {len(c)} strain calls, {len(set(k.split('|')[0] for k in c))} species", flush=True)
wms = {a: c for a, c in wms.items() if c}
wms_sp = {s: set(k.split("|")[0] for k in d) for s, d in wms.items()}
aliases = sorted(set(wms) & set(brad))
print(f"paired usable samples: {aliases}")

def frac_in(wset, bset):
    return len(wset & bset) / len(wset) if wset else float("nan")

lines = ["sample\tn_wms_strains\tstrain_conf_same\tstrain_conf_cross\tspecies_conf_same\tspecies_conf_cross"]
ss_same, ss_cross, sp_same, sp_cross = [], [], [], []
for al in aliases:
    ws, wsp = set(wms[al]), wms_sp[al]
    s_same = frac_in(ws, set(brad[al])); s_cross = np.nanmean([frac_in(ws, set(brad[b])) for b in aliases if b != al])
    p_same = frac_in(wsp, brad_sp[al]);  p_cross = np.nanmean([frac_in(wsp, brad_sp[b]) for b in aliases if b != al])
    lines.append(f"{al}\t{len(ws)}\t{s_same:.3f}\t{s_cross:.3f}\t{p_same:.3f}\t{p_cross:.3f}")
    print(lines[-1])
    ss_same.append(s_same); ss_cross.append(s_cross); sp_same.append(p_same); sp_cross.append(p_cross)
lines.append(f"#mean\t-\t{np.mean(ss_same):.3f}\t{np.mean(ss_cross):.3f}\t{np.mean(sp_same):.3f}\t{np.mean(sp_cross):.3f}")
open(f"{WORK}/results/saliva_concordance.tsv", "w").write("\n".join(lines) + "\n")
print(f"\nstrain-call confirmation by 2bRAD: same-sample {np.mean(ss_same):.2f} vs cross-sample {np.mean(ss_cross):.2f}")
print(f"species-call confirmation by 2bRAD: same-sample {np.mean(sp_same):.2f} vs cross-sample {np.mean(sp_cross):.2f}")

# figure: per-sample WMS strain-call confirmation, same vs cross
fig, ax = plt.subplots(figsize=(7.6, 4.6))
x = np.arange(len(aliases)); w = 0.36
ax.bar(x - w/2, ss_same, w, color="#2ca02c", label="confirmed by SAME-sample 2bRAD")
ax.bar(x + w/2, ss_cross, w, color="#bbbbbb", label="confirmed by other samples' 2bRAD (mean)")
ax.set_xticks(x); ax.set_xticklabels([f"{a}\n(subj {a.split('-')[1]})" for a in aliases], fontsize=9)
ax.set_ylabel("fraction of WMS strain calls confirmed"); ax.set_ylim(0, 1.05)
ax.set_title("Shotgun↔2bRAD saliva concordance: WMS strain calls are\nconfirmed by native 2bRAD on the same sample", fontsize=11)
ax.legend(fontsize=8.5, loc="upper right"); ax.grid(alpha=0.25, axis="y")
fig.tight_layout()
fig.savefig(f"{WORK}/results/saliva_concordance.png", dpi=150)
fig.savefig(f"{WORK}/results/saliva_concordance.pdf")
print("wrote results/saliva_concordance.tsv + .png/.pdf")
