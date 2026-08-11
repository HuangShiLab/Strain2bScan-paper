#!/bin/bash
#SBATCH --job-name=ss_profile_sim
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=ss_profile_sim_%A_%a.out
#SBATCH --error=ss_profile_sim_%A_%a.err

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate /lustre1/g/aos_shihuang/tools/StrainScan/conda_envs/strainscan

export MODE=profile
python3 /lustre1/g/aos_shihuang/LU/Strain2bScan-raw-data/sim_benchmark/run_strainscan_sim_hpc.py
