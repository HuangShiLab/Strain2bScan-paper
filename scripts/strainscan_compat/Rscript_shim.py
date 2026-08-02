#!/usr/bin/env python3
# Minimal Rscript shim for StrainScan's tem_hier.R.
# Place this on PATH as 'Rscript' when R is not available.
import sys, re
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

if len(sys.argv) < 2 or not sys.argv[1].endswith('.R'):
    sys.stderr.write('Rscript shim: only tem_hier.R is supported\n')
    sys.exit(1)

rfile = sys.argv[1]
with open(rfile) as f:
    code = f.read()

m = re.search(r'read\.table\("([^"]+)"', code)
if not m:
    sys.stderr.write('Rscript shim: could not parse input file\n')
    sys.exit(1)
matrix_file = m.group(1)

m = re.search(r'hclust\(d,\s*method="([^"]+)"\)', code)
method = m.group(1) if m else 'average'

m = re.search(r'cutree\(hc,\s*h=([0-9.]+)\)', code)
if not m:
    sys.stderr.write('Rscript shim: could not parse cutoff\n')
    sys.exit(1)
cutoff = float(m.group(1))

df = pd.read_csv(matrix_file, sep='\t', index_col=0)
labels = df.index.tolist()
dist = df.values
for i in range(len(dist)):
    dist[i, i] = 0.0
cd = squareform(dist, checks=False)
Z = linkage(cd, method=method)
cs = fcluster(Z, t=cutoff, criterion='distance')
membership = {labels[i]: int(cs[i]) for i in range(len(labels))}
items = sorted(membership.items(), key=lambda kv: (kv[1], kv[0]))
groups = {}
for name, clst in items:
    groups.setdefault(clst, []).append(name)
out = []
for clst, names in groups.items():
    out.append(' '.join(names))
    out.append(' '.join([str(clst)] * len(names)))
print('\n'.join(out))
