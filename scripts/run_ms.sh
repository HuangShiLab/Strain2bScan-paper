#!/usr/bin/env bash
set -euo pipefail
cd ${WORKDIR:-.}
BIN=${STRAIN2BSCAN_BIN:-strain2bscan}
SET="CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"
GEN=multispecies/genomes; DBS=multispecies/dbs; mkdir -p "$DBS"
export STRAIN2BSCAN_THREADS=16

echo "### 1. BUILD per-species cluster DBs (16 threads) ###"
t0=$SECONDS; nb=0
for sp in $(ls "$GEN"); do
  [ -n "$(ls "$GEN/$sp"/*.fna 2>/dev/null)" ] || continue
  $BIN cluster --genomes "$GEN/$sp" --enzyme "$SET" --out "$DBS/$sp.tsv" --similarity 0.95 >/dev/null 2>&1 && nb=$((nb+1))
done
echo "built $nb species DBs in $((SECONDS-t0))s"

echo "### 2. SPECIES-COUNT gradient: profile 1 sample vs N species DBs (up to the full panel) ###"
ALLSP=($(ls "$DBS"/*.tsv | grep -v members | sed 's#.*/##; s#.tsv##'))
TOTSP=${#ALLSP[@]}
printf '%-10s %-10s\n' "n_species" "time_s"
for N in 10 20 30 40 50 $TOTSP; do
  [ $N -le $TOTSP ] || continue
  d="multispecies/dbs_$N"; rm -rf "$d"; mkdir -p "$d"
  for sp in "${ALLSP[@]:0:$N}"; do ln -sf "$(pwd)/$DBS/$sp.tsv" "$d/$sp.tsv"; done
  s=$SECONDS
  $BIN multi-profile --dbs "$d" --reads multispecies/samples/sample_00.fq --enzyme "$SET" >/dev/null 2>&1
  printf '%-10s %-10s\n' "$N" "$((SECONDS-s))"
done

echo "### 3. SAMPLE-COUNT gradient: profile N samples vs ALL species DBs (cycle the sample pool) ###"
POOL=($(ls multispecies/samples/*.fq)); P=${#POOL[@]}
printf '%-10s %-10s\n' "n_samples" "time_s"
for N in 10 50 100 200 500; do
  s=$SECONDS
  for ((i=0;i<N;i++)); do
    $BIN multi-profile --dbs "$DBS" --reads "${POOL[$((i % P))]}" --enzyme "$SET" >/dev/null 2>&1
  done
  printf '%-10s %-10s\n' "$N" "$((SECONDS-s))"
done
echo "### DONE ###"
