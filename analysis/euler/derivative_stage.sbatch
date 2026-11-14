#!/bin/bash
#SBATCH --job-name=winn017_derivative_stage
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1

set -euo pipefail
set +u
source /cluster/apps/local/env2lmod.sh >/dev/null 2>&1 || true
module purge >/dev/null 2>&1 || true
module load stack/2024-06 gcc/12.2.0 r/4.5.1
set -u
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
: "${WINN_RELEASE_ROOT:?Set WINN_RELEASE_ROOT}"
export R_LIBS_USER="${WINN_RELEASE_ROOT}/Rlib"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export DERIVATIVE_CORES="${SLURM_CPUS_PER_TASK:-1}"
export DERIVATIVE_PREVIOUS_ROOT="${WINN_RELEASE_ROOT}/previous_derivative016"
cd "${WINN_RELEASE_ROOT}"
Rscript analysis/prepare_derivative_rerun_inputs.R "${WINN_RELEASE_ROOT}"
