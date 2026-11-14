#!/bin/bash
#SBATCH --job-name=winn017_gate_validation
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-40%40
set -euo pipefail
set +u
source /cluster/apps/local/env2lmod.sh >/dev/null 2>&1 || true
module purge >/dev/null 2>&1 || true
module load stack/2024-06 gcc/12.2.0 r/4.5.1
set -u
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${WINN_GATE_VALIDATION_ROOT:?Set the isolated validation directory}"
Rscript "${WINN_GATE_VALIDATION_ROOT}/scripts/run_gate_validation.R" \
  "${WINN_GATE_VALIDATION_ROOT}/package" "${WINN_GATE_VALIDATION_ROOT}/output" \
  "${SLURM_ARRAY_TASK_ID}" 40
