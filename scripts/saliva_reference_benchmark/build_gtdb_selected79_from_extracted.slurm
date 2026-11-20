#!/bin/bash
#SBATCH --job-name=s2bs_gtdb79b
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/projects/strain2bscan_reference_benchmark_20261120/logs/gtdb_selected79_build_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/projects/strain2bscan_reference_benchmark_20261120/logs/gtdb_selected79_build_%j.err
set -euo pipefail
P=/lustre1/g/aos_shihuang/projects/strain2bscan_reference_benchmark_20261120
S2BS=/lustre1/g/aos_shihuang/Strain2bScan/target/release/strain2bscan
export STRAIN2BSCAN_THREADS=16
tail -n +2 "$P/results/gtdb_selected79_panel_manifest.tsv" | while IFS=$'\t' read -r safe species dir n; do
  if [ "$n" -lt 2 ]; then
    echo "skip singleton species: $safe $species n=$n"
    continue
  fi
  echo "building $safe ($species): $n genomes"
  "$S2BS" cluster \
    --genomes "$dir" \
    --enzyme BcgI \
    --similarity 0.95 \
    --out "$P/gtdb_bcgI_dbs/${safe}.tsv"
done
"$S2BS" batch \
  --dbs "$P/gtdb_bcgI_dbs" \
  --manifest "$P/saliva_32_manifest.csv" \
  --enzyme BcgI \
  --min-species-markers 50 \
  --min-species-detect 5 \
  --min-support 8 \
  --min-coverage 0.1 \
  --min-consistency 0.5 \
  --out "$P/results/saliva_GTDB_selected79_bcgI.tsv" \
  > "$P/results/saliva_GTDB_selected79_bcgI.stdout.txt"
