#!/bin/bash
#SBATCH --job-name=winn_f2conf
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8

set -euo pipefail
set +u
source /cluster/apps/local/env2lmod.sh >/dev/null 2>&1 || true
module purge >/dev/null 2>&1 || true
module load stack/2024-06 gcc/12.2.0 r/4.5.1
set -u
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
export R_LIBS_USER=/cluster/scratch/ttanna/winn016_release_20261110/Rlib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export F2_ROOT=/cluster/scratch/ttanna/winn_figure2_confounding_20261114
# Resume the same fixed grid. One allocation avoids per-task scheduler delays.
seq 1 190 | xargs -P 8 -I TASK bash -c '
  task=$(printf "%03d" TASK)
  if [ ! -f "$F2_ROOT/results/task_${task}/complete.json" ]; then
    Rscript "$F2_ROOT/code/run.R" /cluster/scratch/ttanna/winn016_release_20261110 \
      "$F2_ROOT/results" TASK "$F2_ROOT/code" \
      > "$F2_ROOT/logs/resume_${task}.out" 2> "$F2_ROOT/logs/resume_${task}.err"
  fi
'
