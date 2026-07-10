#!/usr/bin/env python3
"""Mock host-contamination series: native 2bRAD vs in-silico-digested shotgun (WMS) on the
ATCC MSA-1002 20-species even mock. Precision stays 1.0 for both; at 99% human DNA native 2bRAD
keeps full 20/20 recall while shotgun drops to 12/20, because 2bRAD preserves ~10x more usable
BcgI markers under heavy host contamination."""
import os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

WORK = "/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad"
rows = [l.rstrip("\n").split("\t") for l in open(f"{WORK}/results/mock_hostcontam.tsv")][1:]
d = [dict(dt=r[0], host=int(r[1]), markers=int(r[3]), recall=float(r[7]) * 20, prec=float(r[6])) for r in rows]

def agg(dt):
    xs = sorted(set(x["host"] for x in d if x["dt"] == dt))
    rec = [np.mean([x["recall"] for x in d if x["dt"] == dt and x["host"] == h]) for h in xs]
    mk  = [np.mean([x["markers"] for x in d if x["dt"] == dt and x["host"] == h]) for h in xs]
    return xs, rec, mk

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
C = {"2bRAD": "#1f77b4", "WMS": "#d62728"}
lab = {"2bRAD": "native BcgI 2bRAD", "WMS": "shotgun WMS (in-silico BcgI)"}
HOSTS = [0, 90, 99]                      # categorical, evenly spaced
xpos = {h: i for i, h in enumerate(HOSTS)}
for dt in ("2bRAD", "WMS"):
    xs, rec, mk = agg(dt)
    X = [xpos[h] for h in xs]
    ax1.plot(X, rec, "-o", color=C[dt], lw=2, ms=9, label=lab[dt])
    ax2.plot(X, mk, "-o", color=C[dt], lw=2, ms=9, label=lab[dt])
for ax in (ax1, ax2):
    ax.set_xticks(list(xpos.values())); ax.set_xticklabels([f"{h}%" for h in HOSTS])
    ax.set_xlim(-0.3, 2.3); ax.set_xlabel("human DNA contamination")
ax1.set_ylim(0, 21.5); ax1.set_yticks(range(0, 21, 5))
ax1.set_ylabel("species recall (of 20)")
ax1.set_title("Recall: 2bRAD holds 20/20 at 99% host;\nshotgun drops to 12/20", fontsize=11)
ax1.axhspan(19.5, 21.5, color="#2ca02c", alpha=0.07)
ax1.annotate("precision = 1.00 at every point\n(no false positives, either method)",
             (0.5, 0.06), xycoords="axes fraction", ha="center", fontsize=8.5, color="#555")
ax1.grid(alpha=0.25); ax1.legend(fontsize=9, loc="center left")

ax2.set_yscale("log")
ax2.set_ylabel("usable BcgI markers in sample (log)")
ax2.set_title("Mechanism: 2bRAD preserves ~10× more\nusable markers under host load", fontsize=11)
ax2.grid(alpha=0.25, which="both"); ax2.legend(fontsize=9, loc="lower left")
fig.suptitle("ATCC MSA-1002 mock, host-contamination series — 2bRAD vs shotgun (Strain2bScan)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(f"{WORK}/results/mock_hostcontam.png", dpi=150)
fig.savefig(f"{WORK}/results/mock_hostcontam.pdf")
print("wrote mock_hostcontam.png + .pdf")
