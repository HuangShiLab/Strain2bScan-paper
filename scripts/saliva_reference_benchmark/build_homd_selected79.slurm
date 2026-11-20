#!/bin/bash
#SBATCH --job-name=s2bs_homd79
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/projects/strain2bscan_reference_benchmark_20261120/logs/homd_selected79_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/projects/strain2bscan_reference_benchmark_20261120/logs/homd_selected79_%j.err
set -euo pipefail
P=/lustre1/g/aos_shihuang/projects/strain2bscan_reference_benchmark_20261120
S2BS=/lustre1/g/aos_shihuang/Strain2bScan/target/release/strain2bscan
GENOMES=/lustre1/g/aos_shihuang/databases/HROM/data/genome_catalog/HROM_nonredundant_genomes
export STRAIN2BSCAN_THREADS=16
manifest="$P/results/homd_selected79_panel_manifest.tsv"
echo -e "safe_species\throm_cluster\tspecies_definition\tgtdb_species\tgenome_dir\tn_genome_files" > "$manifest"
tail -n +2 "$P/HOMD_selected79_clusters.tsv" | while IFS=$'\t' read -r cluster species_definition gtdb exact; do
  dir="$GENOMES/$cluster"
  if [ ! -d "$dir" ]; then echo "missing genome dir: $dir"; continue; fi
  n=$(find "$dir" -maxdepth 1 -type f -name '*.fna' | wc -l)
  safe="HOMD__${cluster}__$(printf '%s' "$gtdb" | sed 's/^s__//;s/[^A-Za-z0-9_.]/_/g')"
  echo -e "${safe}\t${cluster}\t${species_definition}\t${gtdb}\t${dir}\t${n}" | tee -a "$manifest"
  if [ "$n" -lt 2 ]; then
    echo "skip singleton HROM species: $cluster n=$n"
    continue
  fi
  "$S2BS" cluster \
    --genomes "$dir" \
    --enzyme BcgI \
    --similarity 0.95 \
    --out "$P/homd_bcgI_dbs/${safe}.tsv"
done
"$S2BS" batch \
  --dbs "$P/homd_bcgI_dbs" \
  --manifest "$P/saliva_32_manifest.csv" \
  --enzyme BcgI \
  --min-species-markers 50 \
  --min-species-detect 5 \
  --min-support 8 \
  --min-coverage 0.1 \
  --min-consistency 0.5 \
  --out "$P/results/saliva_HOMD_selected79_bcgI.tsv" \
  > "$P/results/saliva_HOMD_selected79_bcgI.stdout.txt"
