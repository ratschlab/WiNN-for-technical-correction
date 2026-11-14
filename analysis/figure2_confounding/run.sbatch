#!/bin/bash
#SBATCH --job-name=winn_f2conf
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-190%40
set -euo pipefail
set +u
source /cluster/apps/local/env2lmod.sh >/dev/null 2>&1 || true
module purge >/dev/null 2>&1 || true
module load stack/2024-06 gcc/12.2.0 r/4.5.1
set -u
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
export R_LIBS_USER=/cluster/scratch/ttanna/winn016_release_20261110/Rlib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript /cluster/scratch/ttanna/winn_figure2_confounding_20261114/code/run.R \
 /cluster/scratch/ttanna/winn016_release_20261110 \
 /cluster/scratch/ttanna/winn_figure2_confounding_20261114/results \
 "$SLURM_ARRAY_TASK_ID" /cluster/scratch/ttanna/winn_figure2_confounding_20261114/code
