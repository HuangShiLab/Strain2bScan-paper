import csv, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, numpy as np
src=next(p for p in ["results/cross_species.tsv","xspec/cross_species.tsv"] if os.path.exists(p))
r=list(csv.DictReader(open(src),delimiter="\t"))
sp=[x["species"].replace("_"," ").replace("Cutibacterium","C.").replace("Staphylococcus","S.") for x in r]
P=[float(x["precision"]) for x in r]; R=[float(x["recall"]) for x in r]
cf=[float(x["cluster_frac"]) for x in r]; ncl=[int(x["n_clusters"]) for x in r]; ng=[int(x["n_genomes"]) for x in r]
x=np.arange(len(sp)); w=0.35
fig,ax=plt.subplots(figsize=(7,4.4))
ax.bar(x-w/2,P,w,label="precision",color="#1f77b4")
ax.bar(x+w/2,R,w,label="recall",color="#2ca02c")
ax.set_xticks(x); ax.set_xticklabels(sp); ax.set_ylabel("accuracy"); ax.set_ylim(0,1.05)
ax2=ax.twinx()
ax2.plot(x,cf,"D-",color="#d62728",label="resolvability (clusters/genomes)")
ax2.set_ylabel("clusters / genomes",color="#d62728"); ax2.set_ylim(0,1.0)
for i in x: ax2.annotate(f"{ncl[i]}/{ng[i]}",(i,cf[i]),textcoords="offset points",xytext=(0,7),fontsize=8,color="#d62728",ha="center")
ax.set_title("Cross-species strain profiling: accuracy tracks intra-species resolvability\n(14 enzymes, simulated strain mixtures)")
ax.legend(loc="upper right",fontsize=8); ax2.legend(loc="lower left",fontsize=8)
fig.tight_layout(); os.makedirs("figures",exist_ok=True)
fig.savefig("figures/cross_species.png",dpi=150); fig.savefig("figures/cross_species.pdf")
print("wrote figures/cross_species.png + .pdf")
