#!/bin/bash
#SBATCH --job-name=winn017_derivative_simulation
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=1

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
export WINN_RERUN_METHODS=raw,winn_auto_qc,winn_auto_batch_qc,winn_fixed_no_qc
Rscript analysis/refresh_winn_selection_provenance.R "${WINN_RELEASE_ROOT}"
mkdir -p results/winn_only_aggregate
cp -a results/final/tables results/winn_only_aggregate/
cp results/final/aggregate_manifest.json results/winn_only_aggregate/
Rscript analysis/splice_winn_rerun.R --new=results/winn_only_aggregate --old=previous_final016 --out=results/spliced_aggregate
cp -a results/spliced_aggregate/tables/. results/final/tables/
cp results/spliced_aggregate/validation/winn_rerun_splice_log.csv results/final/validation/
Rscript analysis/repair_ablation_derivatives.R
Rscript analysis/build_execution_ledger.R
Rscript analysis/audit_reference_separation.R
Rscript analysis/validate_release.R
cp results/final/validation/release_validation.csv results/final/validation/rerun_release_validation.csv
Rscript analysis/combine_winn_rerun_validation.R --rerun=results/final/validation/rerun_release_validation.csv --previous=previous_final016/validation/release_validation.csv --out=results/final/validation/release_validation.csv
Rscript analysis/inventory_pre_derivative_tables.R "${WINN_RELEASE_ROOT}"
unset DERIVATIVE_COMPUTE_METHODS
Rscript analysis/run_derivative_biology.R --dataset=simulation
python3 analysis/restore_derivative_comparator_rows.py --new results/final/derivative/simulation --previous "${DERIVATIVE_PREVIOUS_ROOT}/simulation"
