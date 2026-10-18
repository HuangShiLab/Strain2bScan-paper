#!/usr/bin/env python3
import csv,collections,statistics
from pathlib import Path
rows=list(csv.DictReader(open('results/instrain_sim225/persample.tsv'),delimiter='\t'))
out=Path('results/instrain_sim225'); out.mkdir(exist_ok=True)
fields=['group','n','median_precision','median_recall','median_f1','mean_precision','mean_recall','mean_f1']
def write(path,groups):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader()
  for k,rs in sorted(groups.items()):
   vals=[tuple(float(r[m]) for m in ('precision','recall','f1')) for r in rs]
   meds=[statistics.median(v[i] for v in vals) for i in range(3)]
   means=[sum(v[i] for v in vals)/len(vals) for i in range(3)]
   w.writerow({'group':k,'n':len(rs),**{f'median_{m}':meds[i] for i,m in enumerate(('precision','recall','f1'))},**{f'mean_{m}':means[i] for i,m in enumerate(('precision','recall','f1'))}})
bydepth=collections.defaultdict(list); byspecies=collections.defaultdict(list)
for r in rows:
 bydepth[r['sample'].rsplit('_d',1)[-1]].append(r); byspecies[r['species']].append(r)
write(out/'by_depth.tsv',bydepth); write(out/'by_species.tsv',byspecies)
print('wrote',len(rows),'rows')
