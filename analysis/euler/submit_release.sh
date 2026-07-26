#!/bin/bash

set -euo pipefail

: "${WINN_RELEASE_ROOT:?Set WINN_RELEASE_ROOT to the prepared release directory}"
export WINN_CANONICAL_SOURCE_ROOT="${WINN_CANONICAL_SOURCE_ROOT:-${WINN_RELEASE_ROOT}}"
script_root="${WINN_RELEASE_ROOT}/analysis/euler"

submit() {
  local response
  response=$(sbatch --parsable "$@")
  printf '%s\n' "${response%%;*}"
}

primary_fast=$(submit "${script_root}/primary_fast.sbatch")
primary_tiger=$(submit "${script_root}/primary_tiger.sbatch")
ablations=$(submit "${script_root}/ablation_bundles.sbatch")
simulation_fast=$(submit "${script_root}/simulation_methods_fast.sbatch")
simulation_tiger=$(submit "${script_root}/simulation_methods_tiger.sbatch")
simulation_ablations=$(submit "${script_root}/simulation_ablation_bundles.sbatch")
reference_fast=$(submit "${script_root}/reference_splits_fast.sbatch")
reference_tiger=$(submit "${script_root}/reference_splits_tiger.sbatch")
partial_confounding=$(submit "${script_root}/partial_confounding.sbatch")

primary_evaluation=$(submit \
  --dependency="afterok:${primary_fast}:${primary_tiger}" \
  --export=ALL,EVAL_FAMILY=primary \
  "${script_root}/evaluate_dataset_family.sbatch")
ablation_evaluation=$(submit \
  --dependency="afterok:${ablations}" \
  --export=ALL,EVAL_FAMILY=ablations \
  "${script_root}/evaluate_dataset_family.sbatch")
simulation_evaluation=$(submit \
  --dependency="afterok:${simulation_fast}:${simulation_tiger}:${simulation_ablations}" \
  "${script_root}/evaluate_simulation_seeds.sbatch")
reference_evaluation=$(submit \
  --dependency="afterok:${reference_fast}:${reference_tiger}" \
  "${script_root}/evaluate_reference_splits.sbatch")

aggregate=$(submit \
  --dependency="afterok:${primary_evaluation}:${ablation_evaluation}:${simulation_evaluation}:${reference_evaluation}:${partial_confounding}" \
  "${script_root}/aggregate_release.sbatch")

printf '%-27s %s\n' \
  primary_fast "${primary_fast}" \
  primary_tiger "${primary_tiger}" \
  ablations "${ablations}" \
  simulation_fast "${simulation_fast}" \
  simulation_tiger "${simulation_tiger}" \
  simulation_ablations "${simulation_ablations}" \
  reference_fast "${reference_fast}" \
  reference_tiger "${reference_tiger}" \
  partial_confounding "${partial_confounding}" \
  primary_evaluation "${primary_evaluation}" \
  ablation_evaluation "${ablation_evaluation}" \
  simulation_evaluation "${simulation_evaluation}" \
  reference_evaluation "${reference_evaluation}" \
  aggregate_and_validation "${aggregate}"
