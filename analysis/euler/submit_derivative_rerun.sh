#!/bin/bash
# Schedule evaluations of frozen matrices after a WiNN-only correction rerun.
set -euo pipefail
: "${WINN_RELEASE_ROOT:?Set the isolated release directory}"
if [ "$#" -ne 2 ]; then
  printf 'Usage: %s PRIMARY_ARRAY_JOB_ID CORRECTION_AGGREGATE_JOB_ID\n' "$0" >&2
  exit 2
fi
primary_job=$1
aggregate_job=$2
[[ "$primary_job" =~ ^[0-9]+$ && "$aggregate_job" =~ ^[0-9]+$ ]] || exit 2
cd "$WINN_RELEASE_ROOT"
for path in previous_final016/validation/release_validation.csv previous_derivative016 results/frozen_endpoint_free_inputs/provenance/matrix_checksums_pre.csv; do
  test -e "$path" || { printf 'Missing frozen prior input: %s\n' "$path" >&2; exit 1; }
done
submit() {
  local response
  response=$(sbatch --parsable --export=ALL "$@")
  printf '%s\n' "${response%%;*}"
}
stage=$(submit --dependency="afterok:$primary_job" analysis/euler/derivative_stage.sbatch)
empirical=$(submit --dependency="afterok:$stage" analysis/euler/derivative_empirical.sbatch)
lobo=$(submit --dependency="afterok:$stage" analysis/euler/derivative_lobo.sbatch)
simulation=$(submit --dependency="afterok:$aggregate_job:$stage" analysis/euler/derivative_simulation.sbatch)
final=$(submit --dependency="afterok:$simulation:$empirical:$lobo" analysis/euler/derivative_aggregate.sbatch)
printf '%-22s %s\n' stage "$stage" empirical "$empirical" lobo "$lobo" simulation "$simulation" derivative_aggregate "$final"
