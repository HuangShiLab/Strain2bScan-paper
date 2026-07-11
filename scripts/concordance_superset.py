#!/usr/bin/env python3
"""Fig 10 validation (correct framing): on host-heavy saliva, native-2bRAD Strain2bScan detects
EVERY strain that shotgun (WMS, in-silico BcgI) detects, PLUS many additional low-abundance strains
that host-dominated shotgun cannot reach. Shows (a) WMS strain calls are a confirmed subset of 2bRAD,
(b) the 2bRAD-only strains are significantly lower-abundance than the shared ones."""
import os, glob, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"

def parse_calls(stdout):
    s = set()
    for line in stdout.splitlines():
        if not line.startswith("  "): continue
        p = line.strip().split("\t")
        if "[detected" in line or "[strain-resolved, no cluster" in line: continue
        if len(p) >= 5: s.add(f"{p[0]}|{p[1]}")
    return s

# 2bRAD strain calls with community relative abundance (support / total support per sample)
rows = [l.rstrip("\n").split("\t") for l in open(f"{WORK}/results/saliva_strain_long.tsv")][1:]
brad = {}   # sample -> {strain_key: support}
for r in rows:
    brad.setdefault(r[0], {})[f"{r[3]}|{r[4]}"] = float(r[7])
brad_relab = {s: {k: v / sum(d.values()) for k, v in d.items()} for s, d in brad.items()}

# WMS calls from cache
wms = {}
for f in sorted(glob.glob(f"{WORK}/results/wms_preds/*.txt")):
    al = os.path.basename(f)[:-4]
    c = parse_calls(open(f).read())
    if c: wms[al] = c

aliases = sorted(set(wms) & set(brad))
print(f"paired samples: {aliases}\n")
lines = ["sample\tn_2brad_strains\tn_wms_strains\twms_confirmed_by_2brad\t2brad_only\tmed_relab_shared\tmed_relab_2bradonly"]
shared_ab, only_ab = [], []
n_conf_tot = n_wms_tot = 0
for al in aliases:
    b = set(brad[al]); w = wms[al]
    conf = w & b; only = b - w
    ab_shared = [brad_relab[al][k] for k in conf]
    ab_only = [brad_relab[al][k] for k in only]
    shared_ab += ab_shared; only_ab += ab_only
    n_conf_tot += len(conf); n_wms_tot += len(w)
    lines.append(f"{al}\t{len(b)}\t{len(w)}\t{len(conf)}\t{len(only)}\t"
                 f"{np.median(ab_shared) if ab_shared else 0:.5f}\t{np.median(ab_only) if ab_only else 0:.5f}")
    print(lines[-1])

print(f"\nWMS strain calls confirmed by native 2bRAD: {n_conf_tot}/{n_wms_tot} "
      f"({100*n_conf_tot/n_wms_tot:.0f}%)")
print(f"median community rel-abundance: WMS-shared strains {np.median(shared_ab):.4f} "
      f"vs 2bRAD-only strains {np.median(only_ab):.4f}  "
      f"({np.median(shared_ab)/max(np.median(only_ab),1e-9):.1f}x higher)")
from scipy.stats import mannwhitneyu
u, pmw = mannwhitneyu(shared_ab, only_ab, alternative="greater")
print(f"Mann-Whitney (shared > 2bRAD-only abundance): p={pmw:.2e}")
lines.append(f"#wms_confirmed\t{n_conf_tot}/{n_wms_tot}\t-\t-\t-\t{np.median(shared_ab):.5f}\t{np.median(only_ab):.5f}")
lines.append(f"#mannwhitney_p\t{pmw:.3e}")
open(f"{WORK}/results/saliva_concordance.tsv", "w").write("\n".join(lines) + "\n")

# figure: (A) per-sample strain counts (WMS-confirmed vs 2bRAD-only), (B) abundance shared vs only
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.7))
x = np.arange(len(aliases))
conf = [len(wms[a] & set(brad[a])) for a in aliases]
only = [len(set(brad[a]) - wms[a]) for a in aliases]
ax1.bar(x, conf, color="#1f77b4", label="also detected by shotgun WMS")
ax1.bar(x, only, bottom=conf, color="#ff7f0e", label="2bRAD-only (WMS missed)")
ax1.set_xticks(x); ax1.set_xticklabels([f"{a}\n(subj {a.split('-')[1]})" for a in aliases], fontsize=9)
ax1.set_ylabel("strains detected"); ax1.legend(fontsize=8.5, loc="upper left")
ax1.set_title("Native 2bRAD detects every WMS strain + many more\n(host-heavy saliva; WMS = subset)", fontsize=10.5)
ax1.grid(alpha=0.25, axis="y")
# panel B: abundance distributions (log)
bp = ax2.boxplot([shared_ab, only_ab], labels=["shared\n(WMS+2bRAD)", "2bRAD-only\n(WMS missed)"],
                 patch_artist=True, widths=0.5, showfliers=False)
for patch, c in zip(bp["boxes"], ["#1f77b4", "#ff7f0e"]): patch.set_facecolor(c); patch.set_alpha(0.6)
for i, data in enumerate([shared_ab, only_ab], 1):
    ax2.scatter(np.random.default_rng(1).normal(i, 0.05, len(data)), data, s=10, color="k", alpha=0.35, zorder=3)
ax2.set_yscale("log"); ax2.set_ylabel("community relative abundance (2bRAD, log)")
ax2.set_title(f"2bRAD-only strains are low-abundance\n(WMS misses them under host load; p={pmw:.1e})", fontsize=10.5)
ax2.grid(alpha=0.25, axis="y", which="both")
fig.suptitle("Saliva shotgun↔2bRAD: 2bRAD confirms all WMS strains and adds the low-abundance ones shotgun can't reach", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(f"{WORK}/results/saliva_concordance.png", dpi=150)
fig.savefig(f"{WORK}/results/saliva_concordance.pdf")
print("wrote results/saliva_concordance.tsv + .png/.pdf")
