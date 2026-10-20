#!/usr/bin/env python3
"""Figure 3 — enzyme-count sweep on the matched 225-sample simulation benchmark."""
import csv,collections,statistics
import matplotlib;matplotlib.use("Agg")
import matplotlib.pyplot as plt
species=list(csv.DictReader(open('results/enzyme_sweep_225/by_species_enzyme.tsv'),delimiter='\t'))
dbs=list(csv.DictReader(open('results/enzyme_sweep_225/database_provenance.tsv'),delimiter='\t'))
counts=['1','2','4','8','14']
agg={}
for c in counts:
 rs=[r for r in species if r['enzyme_count']==c]
 ee=[r for r in dbs if r['enzyme_count']==c]
 agg[c]={
  'markers':[int(r['strain_specific_markers']) for r in ee],
  'clusters':[int(r['n_clusters']) for r in ee],
 }
for c in counts:
 rs=[r for r in species if r['enzyme_count']==c]
 agg[c]['p']=[statistics.median(float(r[f'median_{m}']) for r in rs) for m in ('precision','recall','f1')]
fig,ax=plt.subplots(1,2,figsize=(11.5,4.4))
x=range(len(counts))
for i,(m,label,color) in enumerate(zip(('precision','recall','f1'),('Precision','Recall','F1'),('#1f77b4','#2ca02c','#ff7f0e'))):
 ax[0].plot(x,[agg[c]['p'][i] for c in counts],marker='o',lw=2,color=color,label=label)
ax[0].set_xticks(x);ax[0].set_xticklabels([f'{c}\nenzyme{"s" if int(c)>1 else ""}' for c in counts])
ax[0].set_ylim(0,1.08);ax[0].set_ylabel('median score across species');ax[0].set_title('A  Accuracy across enzyme ladder')
ax[0].legend(ncol=3,fontsize=8);ax[0].grid(axis='y',alpha=.25)
markers=[statistics.median(agg[c]['markers']) for c in counts]
clusters=[statistics.median(agg[c]['clusters']) for c in counts]
ax[1].plot(x,markers,'o-',lw=2,color='#6a3d9a',label='strain-specific markers')
ax[1].set_xticks(x);ax[1].set_xticklabels([f'{c}\nenzyme{"s" if int(c)>1 else ""}' for c in counts])
ax[1].set_yscale('log');ax[1].set_ylabel('median strain-specific markers');ax[1].set_title('B  Marker yield')
ax2=ax[1].twinx();ax2.plot(x,clusters,'s--',lw=1.5,color='#555555',label='clusters')
ax2.set_ylabel('median clusters',color='#555555');ax2.set_ylim(0,max(clusters)*1.35)
lines1,labels1=ax[1].get_legend_handles_labels();lines2,labels2=ax2.get_legend_handles_labels()
ax[1].legend(lines1+lines2,labels1+labels2,fontsize=8,loc='upper left')
fig.tight_layout()
fig.savefig('figures/source_panels/enzyme_sweep.png',dpi=180,bbox_inches='tight')
fig.savefig('figures/source_panels/enzyme_sweep.pdf',bbox_inches='tight')
print('wrote enzyme sweep figure from 225-sample rerun')
