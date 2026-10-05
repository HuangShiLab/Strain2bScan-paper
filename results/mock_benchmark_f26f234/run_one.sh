#!/bin/bash
# Run one row of wms_analysis_manifest.csv for one mode.
# Usage: run_one.sh <mode> <row_index>
# mode = default | tracegap | layers
set -euo pipefail

BIN=/Users/macstudio/Downloads/Strain2bScan/target/release/strain2bscan
ROOT=/Users/macstudio/Downloads/Strain2bScan-raw-data
MANIFEST=/Users/macstudio/Downloads/Strain2bScan-paper/work/mock_retest/Strain2bScan-raw-data/wms_analysis_manifest.csv
OUT_ROOT=/Users/macstudio/Downloads/Strain2bScan-paper/work/mock_retest/retest_f26f234

MODE=$1
IDX=$2
ROW=$((IDX + 2))

IFS=, read -r sample kind db enzyme reads1 reads2 outdir < <(sed -n "${ROW}p" "$MANIFEST")

DB=$ROOT/$db
EXTRA=""

# layers only for WMS 164/120 panels (MSA1002/MSA1003)
if [ "$MODE" = "layers" ]; then
    if [ "$kind" != "WMS" ] || { [ "$outdir" != "wms164" ] && [ "$outdir" != "shot120" ]; }; then
        exit 0
    fi
    OUTD=$OUT_ROOT/Strain2bScan-port-results/mock/${outdir}_port_layers
    EXTRA="--layer1 cst --layer2 enet"
elif [ "$MODE" = "tracegap" ]; then
    OUTD=$OUT_ROOT/wms_analysis_tracegap/$outdir
    EXTRA="--trace-gap 10 --trace-floor 1e-4"
else
    OUTD=$OUT_ROOT/wms_analysis/$outdir
    EXTRA=""
fi

mkdir -p "$OUTD"
PRED=$OUTD/${sample}.pred
LOG=$OUTD/${sample}.log

if [ -f "$PRED" ] && [ -s "$PRED" ]; then
    echo "[skip] $MODE $sample"
    exit 0
fi

TMPDIR=$(mktemp -d -t s2bs_${MODE}_${sample}.XXXXXX)
TMP=""
cleanup() { [ -n "$TMP" ] && rm -f "$TMP"; rm -rf "$TMPDIR"; }
trap cleanup EXIT

if [ -n "$reads2" ]; then
    TMP=$TMPDIR/merged.fq.gz
    cat "$ROOT/$reads1" "$ROOT/$reads2" > "$TMP"
else
    TMP=$TMPDIR/reads.fq.gz
    ln -s "$ROOT/$reads1" "$TMP"
fi

echo "[run] $MODE $sample $(date +%H:%M:%S)"
"$BIN" profile --db "$DB" --reads "$TMP" --enzyme "$enzyme" $EXTRA --out "$PRED" > "$LOG" 2>&1
echo "[done] $MODE $sample $(date +%H:%M:%S)"
