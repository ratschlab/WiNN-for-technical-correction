#!/bin/bash
#SBATCH --job-name=winn017_derivative_checks
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
set -euo pipefail
set +u
source /cluster/apps/local/env2lmod.sh >/dev/null 2>&1 || true
module purge >/dev/null 2>&1 || true
module load stack/2024-06 gcc/12.2.0 r/4.5.1
set -u
: "${WINN_RELEASE_ROOT:?Set WINN_RELEASE_ROOT}"
export R_LIBS_USER="${WINN_RELEASE_ROOT}/Rlib"
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${WINN_RELEASE_ROOT}"
Rscript analysis/test_ablation_stage_impacts.R
Rscript analysis/test_derivative_outputs.R
