#!/bin/bash
#SBATCH --job-name=winn-review-rfsc-test
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027/logs/adapter_test_%j.out
#SBATCH --error=/cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027/logs/adapter_test_%j.err
source ~/.bash_profile >/dev/null 2>&1
set -euo pipefail
module load r/4.5.1
export R_LIBS_USER=/cluster/project/platt/tanmay/metabolomics_winn_robustness/Rlib_4.5
cd /cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027
Rscript test_qcrfsc_review_adapter.R
