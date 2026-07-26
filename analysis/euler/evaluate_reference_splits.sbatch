#!/bin/bash
#SBATCH --job-name=winn014_ref_eval
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=48G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-60%24

set -euo pipefail
source ~/.bashrc >/tmp/winn_bashrc_${SLURM_JOB_ID}.log 2>&1
source ~/.bash_profile >/tmp/winn_profile_${SLURM_JOB_ID}.log 2>&1
module load r/4.5.1
: "${WINN_RELEASE_ROOT:?Set WINN_RELEASE_ROOT before submitting the job}"
export WINN_CANONICAL_SOURCE_ROOT="${WINN_CANONICAL_SOURCE_ROOT:-${WINN_RELEASE_ROOT}}"
export R_LIBS_USER="${WINN_R_LIBS_USER:-${WINN_RELEASE_ROOT}/Rlib}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "${WINN_RELEASE_ROOT}/scheduler/logs" "${WINN_RELEASE_ROOT}/scheduler/status"

datasets=(simulation mtbls79 clinical_fiams batchcorr_set1 sacurine waveica)
dataset_index=$(( (SLURM_ARRAY_TASK_ID - 1) / 10 ))
split_number=$(( (SLURM_ARRAY_TASK_ID - 1) % 10 + 1 ))
dataset=${datasets[$dataset_index]}
split_id=$(printf 'SPLIT%02d' "${split_number}")
task_id="REFEVAL_${dataset}_${split_id}"
status_path="${WINN_RELEASE_ROOT}/scheduler/status/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.tsv"
printf 'task_id\tjob_id\tarray_index\tstarted\tstatus\n%s\t%s\t%s\t%s\trunning\n' "${task_id}" "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$(date --iso-8601=seconds)" > "${status_path}"
force_args=()
if [[ "${EVAL_FORCE:-0}" == "1" ]]; then force_args+=(--force); fi
if Rscript "${WINN_RELEASE_ROOT}/analysis/evaluate_reference_split.R" --dataset="${dataset}" --split-id="${split_id}" "${force_args[@]}"; then
  printf '%s\t%s\t%s\t%s\tcompleted\n' "${task_id}" "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$(date --iso-8601=seconds)" >> "${status_path}"
else
  code=$?; printf '%s\t%s\t%s\t%s\tfailed_%s\n' "${task_id}" "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "$(date --iso-8601=seconds)" "${code}" >> "${status_path}"; exit "${code}"
fi
