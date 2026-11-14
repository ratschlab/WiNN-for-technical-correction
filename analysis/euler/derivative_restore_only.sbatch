#!/bin/bash
#SBATCH --job-name=winn017_derivative_restore
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1
set -euo pipefail
: "${WINN_RELEASE_ROOT:?Set the isolated release root}"
: "${DERIVATIVE_DATASET:?Set the completed derivative dataset to validate}"
case "$DERIVATIVE_DATASET" in waveica|mtbls79|sacurine|batchcorr_set1|clinical_fiams|simulation) ;; *) exit 2 ;; esac
cd "$WINN_RELEASE_ROOT"
python3 analysis/restore_derivative_comparator_rows.py \
  --new "results/final/derivative/$DERIVATIVE_DATASET" \
  --previous "previous_derivative016/$DERIVATIVE_DATASET"
