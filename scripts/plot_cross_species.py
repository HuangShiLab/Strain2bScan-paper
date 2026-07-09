import csv, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, numpy as np
src = next(p for p in ["results/cross_species.tsv", "xspec/cross_species.tsv"] if os.path.exists(p))
r = list(csv.DictReader(open(src), delimiter="\t"))
sp = [x["species"].replace("_", " ").replace("Cutibacterium", "C.").replace("Staphylococcus", "S.") for x in r]
P = [float(x["precision"]) for x in r]; R = [float(x["recall"]) for x in r]
BC = [float(x["bray_curtis"]) for x in r]; ncl = [int(x["n_clusters"]) for x in r]; ng = [int(x["n_genomes"]) for x in r]
x = np.arange(len(sp)); w = 0.27
fig, ax = plt.subplots(figsize=(7.6, 4.6))
ax.bar(x - w, P, w, label="precision", color="#1f77b4")
ax.bar(x, R, w, label="recall", color="#2ca02c")
ax.bar(x + w, [1 - b for b in BC], w, label="abundance accuracy (1 − Bray–Curtis)", color="#ff7f0e")
ax.set_xticks(x); ax.set_xticklabels(sp); ax.set_ylabel("accuracy"); ax.set_ylim(0, 1.14)
for i in x:
    ax.annotate(f"{ncl[i]}/{ng[i]} clusters", (i, 1.05), ha="center", fontsize=8, color="#555")
ax.set_title("Cross-species strain profiling (strand-fixed): precision 1.0 across species\n"
             "(the earlier resolvability gradient / over-detection was a digestion artifact)")
ax.legend(loc="lower center", fontsize=8, ncol=3)
fig.tight_layout(); os.makedirs("figures", exist_ok=True)
fig.savefig("figures/cross_species.png", dpi=150); fig.savefig("figures/cross_species.pdf")
print("wrote figures/cross_species.png + .pdf")
