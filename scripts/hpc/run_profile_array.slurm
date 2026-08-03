#!/bin/bash
#SBATCH --job-name=s2bs_mock_retest
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=0-03:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=0-42%20

set -euo pipefail

BIN=/lustre1/g/aos_shihuang/Strain2bScan/target/release/strain2bscan
ROOT=/lustre1/g/aos_shihuang/LU/Strain2bScan-raw-data
MANIFEST=$ROOT/wms_analysis_manifest.csv

# pick row (1-based header, task 0 -> line 2)
ROW=$((SLURM_ARRAY_TASK_ID + 2))
IFS=, read -r sample kind db enzyme reads1 reads2 outdir < <(sed -n "${ROW}p" $MANIFEST)

DB=$ROOT/$db
OUTD=$ROOT/wms_analysis/$outdir
mkdir -p $OUTD
PRED=$OUTD/${sample}.pred
LOG=$OUTD/${sample}.log

echo "TASK $SLURM_ARRAY_TASK_ID ROW $ROW: sample=$sample kind=$kind db=$db enzyme=$enzyme outdir=$outdir"

TMP=""
cleanup() { [ -n "$TMP" ] && rm -f "$TMP"; }
trap cleanup EXIT

if [ -n "$reads2" ]; then
    TMP=$(mktemp --suffix=.fq.gz)
    cat "$ROOT/$reads1" "$ROOT/$reads2" > "$TMP"
    READS=$TMP
else
    READS=$ROOT/$reads1
fi

/usr/bin/time -v -o "$LOG.time"     "$BIN" profile --db "$DB" --reads "$READS" --enzyme "$enzyme" --out "$PRED" > "$LOG" 2>&1

echo "DONE: $PRED"
ls -l "$PRED"
