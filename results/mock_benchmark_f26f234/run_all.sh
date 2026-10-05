#!/bin/bash
# Re-run all mock samples for the Fig6/Fig12 benchmark locally.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
RUN_ONE=$SCRIPT_DIR/run_one.sh
MANIFEST=/Users/macstudio/Downloads/Strain2bScan-paper/work/mock_retest/Strain2bScan-raw-data/wms_analysis_manifest.csv
TOTAL=$(( $(wc -l < "$MANIFEST") - 1 ))
LAST=$(( TOTAL - 1 ))
JOBS=${1:-4}

mkdir -p "$SCRIPT_DIR/wms_analysis" "$SCRIPT_DIR/wms_analysis_tracegap"

echo "Starting retest: $TOTAL samples, up to $JOBS concurrent jobs per mode"
echo "Default mode..."
seq 0 $LAST | xargs -P "$JOBS" -n 1 "$RUN_ONE" default

echo "Tracegap mode..."
seq 0 $LAST | xargs -P "$JOBS" -n 1 "$RUN_ONE" tracegap

echo "Layers mode..."
seq 0 $LAST | xargs -P "$JOBS" -n 1 "$RUN_ONE" layers

echo "All done."
