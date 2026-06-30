#!/usr/bin/env bash
set -euo pipefail
cd ${WORKDIR:-.}
BIN=${STRAIN2BSCAN_BIN:-strain2bscan}
G=acnes/genomes; R=acnes/reads; T=acnes/truth; O=acnes/results
# 14 well-defined enzymes (exclude degenerate HaeIV id12, Hin4I id14)
SET14="CspCI,AloI,BsaXI,BaeI,BcgI,CjeI,PpiI,PsrI,BplI,FalI,Bsp24I,CjePI,AlfI,BslFI"

run_mode () {  # $1=label  $2=enzyme-set
  local lab="$1" enz="$2"
  echo "############################################################"
  echo "### MODE=$lab  enzymes=$enz"
  echo "############################################################"
  local db="$O/${lab}.db.tsv"
  $BIN cluster --genomes "$G" --enzyme "$enz" --out "$db" --similarity 0.95 2>&1 | tail -8
  local mem="$O/${lab}.db.members.tsv"
  printf "sample\tTP\tFP\tFN\tprecision\trecall\tF1\tL1\tBrayCurtis\n" > "$O/${lab}.summary.tsv"
  for i in 1 2 3 4 5; do
    $BIN profile --db "$db" --reads "$R/sample${i}.fq" --out "$O/${lab}.s${i}.pred.tsv" >/dev/null 2>&1 || true
    # remap truth genomes -> clusters
    awk -F'\t' 'FNR==NR{if($1!~/^#/)cl[$1]=$2;next} $1~/^#/{next}{c=($1 in cl)?cl[$1]:$1; ab[c]+=$2} END{for(k in ab)printf "%s\t%s\n",k,ab[k]}' \
      "$mem" "$T/sample${i}.truth.tsv" > "$O/${lab}.s${i}.truth.clusters.tsv"
    line=$($BIN evaluate --pred "$O/${lab}.s${i}.pred.tsv" --truth "$O/${lab}.s${i}.truth.clusters.tsv" --present 0.01)
    echo "  sample${i}: $line"
    echo "$line" | sed -E "s/^/sample${i}\t/; s/[A-Za-z_-]+=//g; s/  */\t/g" >> "$O/${lab}.summary.tsv"
  done
}

run_mode "bcgi" "BcgI"
echo
run_mode "multi14" "$SET14"
echo "ALL DONE"
