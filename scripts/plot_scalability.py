#!/usr/bin/env python3
"""Fig 2 - scalability: species gradient (flat), sample gradient (linear), build scaling,
and multi-species accuracy vs the Layer-1 species gate."""
import os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def rows(path): return [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]
# multispecies_scaling.tsv has two sections (species gradient, sample gradient)
sp_n, sp_t, sa_n, sa_t = [], [], [], []
for f in rows("results/multispecies_scaling.tsv"):
    if f[0] == "n_species_dbs" or f[0] == "n_samples": continue
    if len(f) == 2:      # species gradient
        sp_n.append(int(f[0])); sp_t.append(float(f[1]))
    elif len(f) == 3:    # sample gradient
        sa_n.append(int(f[0])); sa_t.append(float(f[1]))
# build scaling from parallel_and_build_scaling.tsv (section B: genomes/time/rss/method)
bg, bt = [], []
for f in rows("results/parallel_and_build_scaling.tsv"):
    if f[0].isdigit() and len(f) >= 2:
        bg.append(int(f[0])); bt.append(float(f[1]))
# gate sweep
gg, gp, gr = [], [], []
for f in rows("results/multispecies_accuracy_gate.tsv"):
    if f[0] in ("species_gate",): continue
    try: gg.append(int(f[0])); gp.append(float(f[1])); gr.append(float(f[2]))
    except ValueError: pass

fig, ax = plt.subplots(2, 2, figsize=(10.5, 8))
ax[0,0].plot(sp_n, sp_t, "o-", color="#9467bd"); ax[0,0].set_ylim(0, max(sp_t)*1.6)
ax[0,0].set_xlabel("# species DBs matched"); ax[0,0].set_ylabel("time / sample (s)")
ax[0,0].set_title("Species-count gradient — flat\n(digest once, match many)")
ax[0,1].plot(sa_n, sa_t, "s-", color="#2ca02c")
ax[0,1].set_xlabel("# samples"); ax[0,1].set_ylabel("total time (s)")
ax[0,1].set_title(f"Sample-count gradient — linear (~{sa_t[-1]/sa_n[-1]:.2f} s/sample)")
ax[1,0].plot(bg, bt, "^-", color="#1f77b4")
ax[1,0].set_xlabel("# genomes"); ax[1,0].set_ylabel("build time (s)")
ax[1,0].set_title("DB-build scaling (MinHash >96 genomes)")
ax[1,1].plot(gg, gp, "o-", label="species precision", color="#1f77b4")
ax[1,1].plot(gg, gr, "s--", label="species recall", color="#2ca02c")
ax[1,1].set_xscale("log"); ax[1,1].set_xlabel("Layer-1 species gate (min species-specific markers)")
ax[1,1].set_ylabel("accuracy"); ax[1,1].set_ylim(0, 1.05); ax[1,1].legend(loc="lower right")
ax[1,1].set_title("Multi-species accuracy vs species gate\n(recall stays 1.0)")
for a in ax.flat: a.grid(alpha=0.3)
fig.tight_layout(); os.makedirs("figures", exist_ok=True)
fig.savefig("figures/scalability.png", dpi=150); fig.savefig("figures/scalability.pdf")
print("wrote figures/scalability.png + .pdf")
