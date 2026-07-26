#!/bin/bash
#SBATCH --job-name=winn014_primary
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=24G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-4,6-13,15-22,24-31,33-40,42-49,51-54%24

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

task_id=$(printf 'PRIMARY_%03d' "${SLURM_ARRAY_TASK_ID}")
status_path="${WINN_RELEASE_ROOT}/scheduler/status/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.tsv"
printf 'task_id\tjob_id\tarray_index\tstarted\tstatus\n%s\t%s\t%s\t%s\trunning\n' \
  "${task_id}" "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$(date --iso-8601=seconds)" > "${status_path}"

if Rscript "${WINN_RELEASE_ROOT}/analysis/run_primary_method.R" --task-id="${task_id}"; then
  printf '%s\t%s\t%s\t%s\tcompleted\n' "${task_id}" "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$(date --iso-8601=seconds)" >> "${status_path}"
else
  code=$?
  printf '%s\t%s\t%s\t%s\tfailed_%s\n' "${task_id}" "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$(date --iso-8601=seconds)" "${code}" >> "${status_path}"
  exit "${code}"
fi
