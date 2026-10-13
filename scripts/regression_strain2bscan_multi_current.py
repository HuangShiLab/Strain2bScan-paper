#!/usr/bin/env python3
"""Regenerate the current Strain2bScan 12-community multi-species benchmark."""
from __future__ import annotations

import argparse, csv, hashlib, json, os, statistics, subprocess
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
DEFAULT_SOFTWARE = PAPER.parent / "Strain2bScan"


def read_tsv(path):
    with path.open() as fh: return list(csv.DictReader(fh, delimiter="\t"))


def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as fh:
        for b in iter(lambda: fh.read(1024*1024), b''): h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--binary',type=Path,default=DEFAULT_SOFTWARE/'target/release/strain2bscan')
    ap.add_argument('--source-repo',type=Path,default=DEFAULT_SOFTWARE)
    ap.add_argument('--dbs',type=Path,default=PAPER/'work/software_regression_current/db')
    ap.add_argument('--work',type=Path,default=PAPER/'work/software_regression_current/predictions')
    ap.add_argument('--threads',type=int,default=16)
    args=ap.parse_args()
    if not args.binary.exists(): raise SystemExit(f'binary not found: {args.binary}')
    manifest=[r for r in read_tsv(PAPER/'work/strainscan_headtohead_rerun/run_manifest.tsv') if r['kind']=='multi']
    if len(manifest)!=12: raise SystemExit(f'expected 12 multi runs, got {len(manifest)}')
    if not args.dbs.exists(): raise SystemExit(f'database directory not found: {args.dbs}')
    species=sorted(p.name for p in args.dbs.iterdir() if (p/'db.tsv').exists())
    if not species: raise SystemExit('no current species databases found')
    db_links=PAPER/'work/software_regression_current/multi_dbs'; db_links.mkdir(parents=True,exist_ok=True)
    for sp in species:
        link=db_links/f'{sp}.tsv'
        target=(args.dbs/sp/'db.tsv').resolve()
        if not target.exists(): raise SystemExit(f'missing database: {target}')
        if link.is_symlink() or link.exists(): link.unlink()
        link.symlink_to(target)
    env={**os.environ,'STRAIN2BSCAN_THREADS':str(args.threads)}
    rows=[]
    for n,row in enumerate(manifest,1):
        pred=args.work/row['sample']/'multi_pred.tsv';pred.parent.mkdir(parents=True,exist_ok=True)
        fifo=pred.parent/'reads.fq'
        if fifo.exists(): fifo.unlink()
        os.mkfifo(fifo)
        proc=subprocess.Popen([str(args.binary),'multi-profile','--dbs',str(db_links),'--reads',str(fifo),'--enzyme','all','--out',str(pred)],env=env)
        with fifo.open('wb') as fifo_handle:
            writer=subprocess.Popen(['gzip','-dc',row['r1'],row['r2']],stdout=fifo_handle,stderr=subprocess.DEVNULL)
        rc=proc.wait();writer.wait()
        if rc: raise RuntimeError(f'{row["sample"]} failed: {rc}')
        calls=set()
        for line in pred.read_text().splitlines():
            if not line or line.startswith('#'): continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)>=9: calls.add(f'{fields[0].strip()}|{fields[1].strip()}')
        truth=set()
        for line in Path(row['truth']).read_text().splitlines():
            if line.startswith('#'): continue
            fields=line.rstrip('\n').split('\t');truth.add(f'{fields[0]}|{fields[2]}')
        tp=len(truth&calls);fp=len(calls-truth);fn=len(truth-calls)
        pr=tp/(tp+fp) if tp+fp else 0.;re=tp/(tp+fn) if tp+fn else 0.;f1=2*pr*re/(pr+re) if pr+re else 0.
        rows.append({'sample':row['sample'],'n_truth_clusters':len(truth),'n_calls':len(calls),'tp':tp,'fp':fp,'fn':fn,'s2b_precision':pr,'s2b_recall':re,'s2b_f1':f1})
        print(f'[{n:02d}/12] {row["sample"]} P={pr:.3f} R={re:.3f} F1={f1:.3f}',flush=True)
    out=PAPER/'results/strainscan_rerun/current_multi_persample.tsv'
    with out.open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    from collections import defaultdict
    groups=defaultdict(list)
    for r in rows:groups[r['sample'].split('_')[1]].append(r)
    aggregate=PAPER/'results/strainscan_rerun/current_multi_by_depth.tsv'
    with aggregate.open('w',newline='') as fh:
        fields=['group','n','s2b_precision','s2b_recall','s2b_f1'];w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for depth in ('low','med','high'):
            rr=groups[depth];w.writerow({'group':depth,'n':len(rr),'s2b_precision':statistics.mean(float(x['s2b_precision']) for x in rr),'s2b_recall':statistics.mean(float(x['s2b_recall']) for x in rr),'s2b_f1':statistics.mean(float(x['s2b_f1']) for x in rr)})
    commit=subprocess.check_output(['git','-C',str(args.source_repo),'rev-parse','HEAD'],text=True).strip()
    provenance={'source_commit':commit,'binary':str(args.binary),'binary_sha256':sha256(args.binary),'runs':len(rows),'database_hashes':{p.stem.replace('.tsv',''):sha256(p) for p in db_links.glob('*.tsv')},'single_source':'results/strainscan_rerun/current_single_persample.tsv','multi_source':str(out),'aggregate':str(aggregate)}
    (PAPER/'results/strainscan_rerun/current_multi_software_regression_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(f'wrote {out} and {aggregate}')

if __name__=='__main__':raise SystemExit(main())
