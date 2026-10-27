#!/bin/bash
#SBATCH --job-name=winn-review-rfsc-smoke
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027/logs/smoke_%j.out
#SBATCH --error=/cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027/logs/smoke_%j.err
source ~/.bash_profile >/dev/null 2>&1
set -euo pipefail
module load r/4.5.1
export R_LIBS_USER=/cluster/project/platt/tanmay/metabolomics_winn_robustness/Rlib_4.5
cd /cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027
Rscript run_qcrfsc_revision.R --dataset=sacurine --smoke
