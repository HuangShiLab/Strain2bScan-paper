#!/bin/bash
set -euo pipefail
WORK=/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/work/monkey_panel_20261023
BIN=/lustre1/g/aos_shihuang/Strain2bScan/target/release/strain2bscan
ls -l "$WORK"/profile_shotgun/login_manifest_*.csv
sleep 5
pids=()
for i in 00 01 02 03 04 05 06 07; do
  STRAIN2BSCAN_THREADS=8 "$BIN" batch \
    --dbs "$WORK/dbs" \
    --manifest "$WORK/profile_shotgun/login_manifest_${i}.csv" \
    --out "$WORK/profile_shotgun/login_part${i}.tsv" \
    --enzyme BcgI \
    --min-species-markers 20 \
    --min-species-detect 2 \
    --min-support 2 \
    --min-coverage 0.01 \
    --min-consistency 0.5 \
    > "$WORK/profile_logs/login_shotgun_${i}.out" 2> "$WORK/profile_logs/login_shotgun_${i}.err" &
  pids+=($!)
done
rc=0
for p in "${pids[@]}"; do wait "$p" || rc=1; done
python3 - <<'PY'
import pandas as pd
from pathlib import Path
w=Path('/lustre1/g/aos_shihuang/projects/monkey_strain2bscan/work/monkey_panel_20261023')
frames=[]
for p in sorted((w/'profile_shotgun').glob('login_part*.tsv')):
    if p.stat().st_size:
        frames.append(pd.read_csv(p,sep='\t'))
pd.concat(frames,ignore_index=True).to_csv(w/'profile_shotgun/merged.tsv',sep='\t',index=False)
print('merged rows',sum(len(x) for x in frames))
PY
exit $rc
