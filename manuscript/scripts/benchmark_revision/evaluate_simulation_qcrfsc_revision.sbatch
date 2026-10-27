#!/bin/bash
#SBATCH --job-name=winn-review-rfsc-sim-eval
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-10%10
#SBATCH --output=/cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027/logs/sim_rfsc_eval_%A_%a.out
#SBATCH --error=/cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027/logs/sim_rfsc_eval_%A_%a.err
source ~/.bash_profile >/dev/null 2>&1
set -euo pipefail
module load r/4.5.1
export R_LIBS_USER=/cluster/project/platt/tanmay/metabolomics_winn_robustness/Rlib_4.5
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027
printf -v seed_id 'SIM%02d' "$SLURM_ARRAY_TASK_ID"
Rscript evaluate_simulation_qcrfsc_revision.R --seed-id="$seed_id"
