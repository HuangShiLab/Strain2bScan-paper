#!/bin/bash
set -uo pipefail
WORK=/private/tmp/claude-501/-Users-macstudio-Downloads-YangJiazhen/091459b5-4e03-49b6-8502-3f2acf59ff13/scratchpad
BIN=/Users/macstudio/Downloads/Strain2bScan/target/release/strain2bscan
GEN=$WORK/genomes_oral; DBS=$WORK/dbs/oral; mkdir -p "$DBS"
export STRAIN2BSCAN_THREADS=16
t0=$SECONDS; nb=0
for d in "$GEN"/*/; do
  sp=$(basename "$d")
  ls "$d"/*.fna >/dev/null 2>&1 || { echo "skip $sp"; continue; }
  ng=$(ls "$d"/*.fna | wc -l | tr -d ' ')
  if $BIN cluster --genomes "$d" --enzyme BcgI --out "$DBS/$sp.tsv" --similarity 0.95 >/dev/null 2>&1; then
    nb=$((nb+1)); echo "built $sp ($ng genomes)"
  else echo "FAIL $sp"; fi
done
echo "built $nb oral DBs in $((SECONDS-t0))s"
