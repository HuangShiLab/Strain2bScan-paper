import csv, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src=next(p for p in ["results/enzyme_sweep.tsv","enzsweep/enzyme_sweep.tsv"] if os.path.exists(p))
r=list(csv.DictReader(open(src),delimiter="\t"))
ne=[int(x["n_enzymes"]) for x in r]; ss=[int(x["strain_specific"]) for x in r]; ncl=[int(x["n_clusters"]) for x in r]
P=[float(x["precision"]) for x in r]; R=[float(x["recall"]) for x in r]; BC=[float(x["bray_curtis"]) for x in r]
bt=[float(x["build_s"]) for x in r]
fig,(a1,a2)=plt.subplots(1,2,figsize=(11,4.3))
# Panel A: resolution + marker yield
a1.plot(ne,ncl,"o-",color="#9467bd",label="# clusters (resolution)")
a1.axhline(64,ls=":",color="grey",lw=1); a1.text(13,62,"64 genomes",fontsize=7,color="grey",ha="right")
a1.set_xscale("log",base=2); a1.set_xticks(ne); a1.set_xticklabels(ne)
a1.set_xlabel("# enzymes"); a1.set_ylabel("# clusters",color="#9467bd"); a1.set_ylim(0,70)
a1b=a1.twinx(); a1b.plot(ne,ss,"s--",color="#8c564b",label="strain-specific markers")
a1b.set_ylabel("strain-specific markers",color="#8c564b")
a1.set_title("Resolution & marker yield vs #enzymes")
a1.axvspan(0.9,3,color="orange",alpha=0.10); a1.text(1.4,8,"under-\nresolved",fontsize=7,color="#cc7000",ha="center")
# Panel B: accuracy + cost
a2.plot(ne,P,"o-",color="#1f77b4",label="precision")
a2.plot(ne,R,"s-",color="#2ca02c",label="recall")
a2.plot(ne,BC,"^--",color="#d62728",label="Bray-Curtis (err)")
a2.set_xscale("log",base=2); a2.set_xticks(ne); a2.set_xticklabels(ne)
a2.set_xlabel("# enzymes"); a2.set_ylabel("accuracy metric"); a2.set_ylim(0,1.05)
a2.axvspan(0.9,3,color="orange",alpha=0.10); a2.text(1.4,0.5,"degenerate\n(coarse / under-marked)",fontsize=7,color="#cc7000",ha="center")
a2.annotate("sweet spot",(8,1.0),xytext=(8,0.78),fontsize=8,ha="center",arrowprops=dict(arrowstyle="->"))
a2b=a2.twinx(); a2b.plot(ne,bt,"d:",color="#7f7f7f",label="build time (s)"); a2b.set_ylabel("build time (s)",color="#7f7f7f")
a2.set_title("Accuracy & cost vs #enzymes (C. acnes)")
a1.legend(loc="center right",fontsize=7); a2.legend(loc="center left",fontsize=7)
fig.tight_layout(); os.makedirs("figures",exist_ok=True)
fig.savefig("figures/enzyme_sweep.png",dpi=150); fig.savefig("figures/enzyme_sweep.pdf")
print("wrote figures/enzyme_sweep.png + .pdf")
