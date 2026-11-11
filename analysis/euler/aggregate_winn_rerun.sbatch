#!/bin/bash
#SBATCH --job-name=winn016_aggregate_rerun
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=48G
#SBATCH --cpus-per-task=1

set -euo pipefail
# Clean Lmod initialisation; personal shell profiles are not sourced so that no
# user LD_LIBRARY_PATH or module defaults leak into the frozen environment.
set +u
source /cluster/apps/local/env2lmod.sh >/dev/null 2>&1 || true
module purge >/dev/null 2>&1 || true
module load stack/2024-06 gcc/12.2.0 r/4.5.1
set -u
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
: "${WINN_RELEASE_ROOT:?Set WINN_RELEASE_ROOT before submitting the job}"
export WINN_CANONICAL_SOURCE_ROOT="${WINN_CANONICAL_SOURCE_ROOT:-${WINN_RELEASE_ROOT}}"
export R_LIBS_USER="${WINN_R_LIBS_USER:-${WINN_RELEASE_ROOT}/Rlib}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# WiNN-only re-run: aggregate the re-run methods (WINN_RERUN_METHODS). The
# competitor reuse audit and full-release validation assume all nine methods were
# executed here, so they are replaced by analysis/splice_winn_rerun.R, which merges
# the frozen competitor rows and validates the combined tables.
: "${WINN_RERUN_METHODS:?WINN_RERUN_METHODS must be set for the WiNN-only aggregate}"
Rscript "${WINN_RELEASE_ROOT}/analysis/aggregate_release.R"
set +e
Rscript "${WINN_RELEASE_ROOT}/analysis/build_execution_ledger.R" > "${WINN_RELEASE_ROOT}/results/final/validation/build_execution_ledger.log" 2>&1
echo "build_execution_ledger exit $?"
Rscript "${WINN_RELEASE_ROOT}/analysis/audit_reference_separation.R" > "${WINN_RELEASE_ROOT}/results/final/validation/audit_reference_separation.log" 2>&1
echo "audit_reference_separation exit $?"
set -e
