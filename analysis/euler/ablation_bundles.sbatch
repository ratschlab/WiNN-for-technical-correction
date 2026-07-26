#!/bin/bash
#SBATCH --job-name=winn014_ablation
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-6%6

set -euo pipefail
source ~/.bashrc >/tmp/winn_bashrc_${SLURM_JOB_ID}.log 2>&1
source ~/.bash_profile >/tmp/winn_profile_${SLURM_JOB_ID}.log 2>&1
module load r/4.5.1

: "${WINN_RELEASE_ROOT:?Set WINN_RELEASE_ROOT before submitting the job}"
export WINN_CANONICAL_SOURCE_ROOT="${WINN_CANONICAL_SOURCE_ROOT:-${WINN_RELEASE_ROOT}}"
export R_LIBS_USER="${WINN_R_LIBS_USER:-${WINN_RELEASE_ROOT}/Rlib}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p "${WINN_RELEASE_ROOT}/scheduler/logs" "${WINN_RELEASE_ROOT}/scheduler/status"

datasets=(simulation mtbls79 clinical_fiams batchcorr_set1 sacurine waveica)
dataset=${datasets[$((SLURM_ARRAY_TASK_ID - 1))]}
status_path="${WINN_RELEASE_ROOT}/scheduler/status/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.tsv"
printf 'task_id\tjob_id\tarray_index\tstarted\tstatus\nABLATION_%s\t%s\t%s\t%s\trunning\n' \
  "${dataset}" "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$(date --iso-8601=seconds)" > "${status_path}"

if Rscript "${WINN_RELEASE_ROOT}/analysis/run_ablation_bundle.R" --dataset="${dataset}"; then
  printf 'ABLATION_%s\t%s\t%s\t%s\tcompleted\n' "${dataset}" "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$(date --iso-8601=seconds)" >> "${status_path}"
else
  code=$?
  printf 'ABLATION_%s\t%s\t%s\t%s\tfailed_%s\n' "${dataset}" "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$(date --iso-8601=seconds)" "${code}" >> "${status_path}"
  exit "${code}"
fi
