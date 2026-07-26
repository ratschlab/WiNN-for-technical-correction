# WiNN benchmark analyses

## Overview

This repository contains the analysis code and compact result tables for the WiNN metabolomics technical-correction study. The benchmark covers a truth-known simulation, MTBLS79, an internal pair-aware FIA-MS cohort, BatchCorrMetabolomics Set 1, Sacurine, and the WaveICA adenocarcinoma plasma dataset.

Nine method variants are evaluated in a fixed order: Raw, ComBat, QC-RLSC, QC-RFSC, TIGER, SERRF, WiNN auto with supplied batches and training QCs, WiNN auto-batch with training QCs, and fixed/default WiNN with supplied batches but no QC identities. The fixed/default mode is the prespecified primary WiNN configuration. The automatic modes are reported separately.

Competitor parameters are fixed from published or package-native defaults, native training procedures, or training-QC-only criteria. Held-out QCs, biological labels, replicate identities, and final evaluation metrics are not available during parameter selection. The exact choices and their provenance are recorded in `analysis/config/endpoint_free_selection_manifest.csv`.

## Repository structure

- `analysis/`: task runners, evaluators, aggregation code, validation checks, and Euler job arrays.
- `analysis/config/`: the frozen run matrix, method settings, QC splits, task manifests, and metric definitions.
- `config/holdouts/`: fixed primary held-out-QC definitions for the public datasets.
- `scripts/`: public-data download and preprocessing code plus shared method and metric helpers.
- `data/simulated/canonical/SIM01/`: the canonical simulation used by both the primary comparison and the main ablation.
- `reports/`: suitability audits for the added public datasets.
- `results/final/`: compact final tables, figures, figure source data, and validation summaries.
- `package/`: the exact WiNN 0.1.4 source archive and its validation record.

The production flow is:

```text
download and preprocess data
        -> execute method tasks
        -> evaluate saved matrices
        -> aggregate dataset tables
        -> render figures and reports
```

## Datasets

The five public or synthetic datasets can be prepared from this repository. The internal FIA-MS cohort is not public; only aggregate results are distributed. Its correction unit is the log-scale mean of each adjacent two-injection pair, which prevents the injection design from being interpreted as lag-1 drift. The generic pair-collapse implementation is included, but the repository contains no clinical matrix, identifiers, feature names, row-level results, or split assignments. Clinical tasks therefore require the restricted inputs and assignments held by the study team; the public split audit reports only aggregate counts and hashes.

Dataset sources, assay details, suitability limitations, and preprocessing decisions are documented in [docs/DATASETS.md](docs/DATASETS.md).

## Software requirements

The production runs used R 4.5.1 and WiNN 0.1.4. Restore the recorded R environment and install the bundled package archive into the project library:

```bash
Rscript -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv")'
Rscript -e 'renv::restore(prompt=FALSE)'
mkdir -p Rlib
R CMD INSTALL --library=Rlib package/winn_0.1.4.tar.gz
```

The package archive has SHA-256
`71a0964cee2778b2e5789d20621147e074c7945e813cf76af2ceeb696104aae1`.
The source commit, build information, smoke tests, and package checks are in `package/frozen_package.json` and `package/VALIDATION.md`.

## Downloading and preprocessing public data

These commands download processed matrices and supporting metadata, not raw spectra:

```bash
python3 scripts/download_mtbls79.py
Rscript scripts/preprocess_mtbls79_public_data.R

python3 scripts/download_batchcorr_set1.py
Rscript scripts/preprocess_batchcorr_set1.R

python3 scripts/download_sacurine.py
Rscript scripts/preprocess_sacurine.R

python3 scripts/download_waveica_adenocarcinoma.py
Rscript scripts/preprocess_waveica_adenocarcinoma.R
```

The downloaders are idempotent and write checksums and source manifests. The preprocessing scripts preserve sample order and feature identifiers, validate batch/run-order/QC fields, apply the declared missing-data rules, and reject misaligned matrices.

Generate the remaining prespecified simulation realizations with:

```bash
make simulations
```

The committed `SIM01` bundle is verified and reused; `SIM02` through `SIM10` are generated from `analysis/config/simulation_seed_ledger.csv`. Generated matrices are checked against the frozen identities in `analysis/config/simulation_bundle_hashes.csv`. The cross-platform gate hashes matrices after rounding to nine decimal places because elementary floating-point operations can differ by a few units in the last place between R platforms; exact reference-platform object hashes are retained as provenance diagnostics.

## Reproducing the analyses

The committed final tables and figures can be checked and rendered without recomputing method matrices:

```bash
make validate-results
make report
```

A lightweight package-and-input smoke test is available after installing the project library:

```bash
make smoke
```

Full correction runs require the processed inputs and method dependencies. Set the repository as both the release and canonical-data root. The public workflow recomputes every method because large correction matrices are not distributed:

```bash
export WINN_RELEASE_ROOT="$PWD"
export WINN_CANONICAL_SOURCE_ROOT="$PWD"
export R_LIBS_USER="$PWD/Rlib"

Rscript analysis/run_primary_method.R --task-id=PRIMARY_001
Rscript analysis/evaluate_dataset_family.R --dataset=simulation --family=primary
```

`analysis/config/primary_run_matrix.csv` defines the 54 primary tasks. The other manifests define 10 simulation seeds, 10 reference splits per dataset, nine WiNN ablation variants, and 160 WiNN-only partial-confounding scenarios. Each task writes to its own directory and can be safely retried. Aggregation refuses to run until the declared output counts are complete.

`analysis/config/competitor_reuse_audit.tsv` distinguishes the provenance of the released Euler results from this clean-checkout behavior. The released analysis reused only checksum-validated endpoint-free artifacts; the public workflow refits those methods because the large matrices themselves are not redistributed.

Corrected matrices are required to remain on the nonnegative intensity scale. If a regression-based competitor extrapolates below zero, the shared wrapper floors those values at zero before training-QC scoring or evaluation and records the affected count and original minimum. Non-finite outputs are rejected.

## Running on Euler

The scheduler scripts in `analysis/euler/` use Slurm arrays. On Euler, load the shell profile and R 4.5.1, define the release and data roots, and submit all independent families and their evaluation dependencies with:

```bash
source ~/.bashrc
source ~/.bash_profile
module load r/4.5.1
export WINN_RELEASE_ROOT="$PWD"
export WINN_CANONICAL_SOURCE_ROOT="$PWD"
bash analysis/euler/submit_release.sh
```

The submission command prints every Slurm job ID. Correction, ablation, simulation, reference-split, and confounding arrays are launched together; evaluations wait only for their own inputs, and aggregation waits for all evaluation families. Individual successful tasks are idempotent, so a failed array element can be resubmitted without overwriting verified outputs.

The full task matrix is compute-intensive because TIGER uses its native 4 x 4 ensemble with 500 trees. Fast methods generally finish in minutes; TIGER determines the wall time. Reference stability consists of exactly 10 splits x 9 methods x 6 datasets (540 tasks). The partial-confounding grid contains 10 seeds x 16 scenarios and runs only Raw and fixed/default WiNN.

## Outputs

`results/final/tables/` contains the per-dataset and combined primary results, cumulative-stage impact tables, selective-versus-forced contrasts, ten-seed simulation summaries, ten-split reference summaries, runtime/retention data, and partial-confounding results. Training-QC candidate rankings are retained for every competitor selected by that procedure; the current selected settings and training-QC scores are saved for each WiNN automatic mode. Every figure in `results/final/figures/` has a numerical source table in `results/final/figures/source_data/`. When an exact competitor matrix was reused, cache-loading time is retained as execution provenance but is not reported as method runtime; unrepeated fit times are therefore `NA` rather than artificially near zero.

The metric crosswalk in `analysis/config/metric_crosswalk.csv` states the dataset scope, analysis unit, feature panel, missing-value handling, direction, implementation, and whether a metric was available during selection. Batch weighted-PC R-squared is always calculated with batch treated as a categorical factor. True metabolite ICC(A,1) and the between/within feature repeatability ratio are kept as separate metrics rather than sharing one label.

## Data availability

MTBLS79, BatchCorrMetabolomics Set 1, Sacurine, and WaveICA are obtained from the public sources listed in [docs/DATASETS.md](docs/DATASETS.md). The clinical cohort is subject to its original access restrictions and is not redistributed. Its public outputs are aggregate statistics only.

## Citation

Please cite the WiNN manuscript and the original source publication for each dataset used. Package citation information is available with `citation("winn")` after installation.

## License

The analysis code is distributed under the license in [LICENSE](LICENSE). Source datasets retain their original licenses and terms.

## Contact

Questions about the benchmark or access to restricted study material should be directed to the manuscript authors.
