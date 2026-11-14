#!/bin/bash
#SBATCH --job-name=winn017_derivative_empirical
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-5%5

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
datasets=(waveica batchcorr_set1 clinical_fiams mtbls79 sacurine)
dataset="${datasets[$((SLURM_ARRAY_TASK_ID-1))]}"
case "${dataset}" in
 waveica|batchcorr_set1|clinical_fiams) export DERIVATIVE_COMPUTE_METHODS=winn_auto_qc,winn_auto_batch_qc,winn_fixed_no_qc ;;
 *) unset DERIVATIVE_COMPUTE_METHODS ;;
esac
Rscript analysis/run_derivative_biology.R --dataset="${dataset}"
# Restore the exact historical serialization after checking unchanged estimates.
python3 analysis/restore_derivative_comparator_rows.py --new "results/final/derivative/${dataset}" --previous "${DERIVATIVE_PREVIOUS_ROOT}/${dataset}"
