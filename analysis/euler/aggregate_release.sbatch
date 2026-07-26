#!/bin/bash
#SBATCH --job-name=winn014_aggregate
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=48G
#SBATCH --cpus-per-task=1

set -euo pipefail
source ~/.bashrc >/tmp/winn_bashrc_${SLURM_JOB_ID}.log 2>&1
source ~/.bash_profile >/tmp/winn_profile_${SLURM_JOB_ID}.log 2>&1
module load r/4.5.1
: "${WINN_RELEASE_ROOT:?Set WINN_RELEASE_ROOT before submitting the job}"
export WINN_CANONICAL_SOURCE_ROOT="${WINN_CANONICAL_SOURCE_ROOT:-${WINN_RELEASE_ROOT}}"
export R_LIBS_USER="${WINN_R_LIBS_USER:-${WINN_RELEASE_ROOT}/Rlib}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

Rscript "${WINN_RELEASE_ROOT}/analysis/finalize_reuse_audit.R"
Rscript "${WINN_RELEASE_ROOT}/analysis/aggregate_release.R"
Rscript "${WINN_RELEASE_ROOT}/analysis/build_execution_ledger.R"
Rscript "${WINN_RELEASE_ROOT}/analysis/audit_reference_separation.R"
Rscript "${WINN_RELEASE_ROOT}/analysis/validate_release.R"
