#!/bin/bash
#SBATCH --job-name=winn-review-rfsc-eval
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-5%5
#SBATCH --output=/cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027/logs/rfsc_eval_%A_%a.out
#SBATCH --error=/cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027/logs/rfsc_eval_%A_%a.err
source ~/.bash_profile >/dev/null 2>&1
set -euo pipefail
module load r/4.5.1
export R_LIBS_USER=/cluster/project/platt/tanmay/metabolomics_winn_robustness/Rlib_4.5
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027
keys=(mtbls79 clinical_fiams batchcorr_set1 sacurine waveica)
Rscript evaluate_empirical_qcrfsc_revision.R --dataset="${keys[$((SLURM_ARRAY_TASK_ID-1))]}"
