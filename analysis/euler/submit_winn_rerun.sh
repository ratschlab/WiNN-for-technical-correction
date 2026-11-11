#!/bin/bash
# Re-run only WiNN (all three modes, ablations, partial confounding) plus the
# no-correction Raw reference, and evaluate only those methods. Competitor
# (ComBat, QC-RLSC, QC-RFSC, TIGER, SERRF) outputs are NOT recomputed; their
# frozen aggregate rows are spliced back in afterwards by analysis/splice_winn_rerun.R.
set -euo pipefail

: "${WINN_RELEASE_ROOT:?Set WINN_RELEASE_ROOT to the prepared release directory}"
export WINN_CANONICAL_SOURCE_ROOT="${WINN_CANONICAL_SOURCE_ROOT:-${WINN_RELEASE_ROOT}}"
export WINN_RERUN_METHODS="${WINN_RERUN_METHODS:-raw,winn_auto_qc,winn_auto_batch_qc,winn_fixed_no_qc}"
script_root="${WINN_RELEASE_ROOT}/analysis/euler"
config="${WINN_RELEASE_ROOT}/analysis/config"

# Array indices (task number suffix) of the selected methods in a manifest.
indices() {
  local manifest=$1 prefix=$2
  Rscript -e '
    a <- commandArgs(TRUE); m <- read.csv(a[1], stringsAsFactors = FALSE)
    keep <- strsplit(a[3], ",", fixed = TRUE)[[1]]
    id <- as.integer(sub(paste0("^", a[2]), "", m$task_id[m$method_id %in% keep]))
    cat(paste(sort(id), collapse = ","))' "$manifest" "$prefix" "$WINN_RERUN_METHODS"
}
primary_idx=$(indices "${config}/primary_run_matrix.csv" "PRIMARY_")
simulation_idx=$(indices "${config}/simulation_task_manifest.csv" "SIM_")
reference_idx=$(indices "${config}/reference_split_task_manifest.csv" "REF_")
printf 'primary tasks: %s\nsimulation tasks: %s\nreference tasks: %s\n' \
  "$(tr ',' '\n' <<< "$primary_idx" | wc -l)" "$(tr ',' '\n' <<< "$simulation_idx" | wc -l)" \
  "$(tr ',' '\n' <<< "$reference_idx" | wc -l)"

submit() {
  local response
  response=$(sbatch --parsable --export=ALL "$@")
  printf '%s\n' "${response%%;*}"
}

primary=$(submit --array="${primary_idx}%24" "${script_root}/primary_fast.sbatch")
simulation=$(submit --array="${simulation_idx}%30" "${script_root}/simulation_methods_fast.sbatch")
reference=$(submit --array="${reference_idx}%60" "${script_root}/reference_splits_fast.sbatch")
ablations=$(submit "${script_root}/ablation_bundles.sbatch")
simulation_ablations=$(submit "${script_root}/simulation_ablation_bundles.sbatch")
partial_confounding=$(submit "${script_root}/partial_confounding.sbatch")

primary_evaluation=$(submit --dependency="afterok:${primary}" \
  --export=ALL,EVAL_FAMILY=primary "${script_root}/evaluate_dataset_family.sbatch")
ablation_evaluation=$(submit --dependency="afterok:${ablations}" \
  --export=ALL,EVAL_FAMILY=ablations "${script_root}/evaluate_dataset_family.sbatch")
simulation_evaluation=$(submit --dependency="afterok:${simulation}:${simulation_ablations}" \
  "${script_root}/evaluate_simulation_seeds.sbatch")
reference_evaluation=$(submit --dependency="afterok:${reference}" \
  "${script_root}/evaluate_reference_splits.sbatch")
aggregate=$(submit \
  --dependency="afterok:${primary_evaluation}:${ablation_evaluation}:${simulation_evaluation}:${reference_evaluation}:${partial_confounding}" \
  "${script_root}/aggregate_winn_rerun.sbatch")

printf '%-22s %s\n' \
  primary "$primary" simulation "$simulation" reference "$reference" \
  ablations "$ablations" simulation_ablations "$simulation_ablations" \
  partial_confounding "$partial_confounding" primary_evaluation "$primary_evaluation" \
  ablation_evaluation "$ablation_evaluation" simulation_evaluation "$simulation_evaluation" \
  reference_evaluation "$reference_evaluation" aggregate "$aggregate"
