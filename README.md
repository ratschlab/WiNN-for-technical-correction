# WiNN manuscript and reproducible analyses

This revision contains the code, numerical results, figure source data, compiled figures and matching main/supplementary LaTeX sources for the WiNN metabolomics manuscript. A subsequent documentation commit records this immutable analysis revision in those sources; the analysis and numerical results are frozen before that citation update. The current publication bundle is in [`manuscript/`](manuscript/); its numerical results and complete file hashes are checked with `make validate-results`. WiNN is available on [CRAN](https://CRAN.R-project.org/package=winn). The benchmark retains the frozen WiNN 0.1.4 implementation; the package availability statement refers to CRAN version 0.1.5.

The primary configuration is called **WiNN**. The supplementary configurations use QC-based tuning or inferred acquisition segments and are identified where a comparison requires that distinction. Stable internal method identifiers are preserved so that existing analysis records remain joinable.

## Current results and provenance

The current release replaces the earlier result tree at commit [`21de7c448bfa5545d27c5a930161fb0bd8e40d85`](https://github.com/ratschlab/WiNN-for-technical-correction/tree/21de7c448bfa5545d27c5a930161fb0bd8e40d85). The earlier snapshot remains available through that immutable revision.

- The actual QC-RFSC wrapper in `scripts/benchmark_helpers.R` assigns neutral identifiers to every injection and marks only authorized training references as QCs. Returned relative abundances are restored using raw training-QC feature medians. Features not returned by the native implementation retain their input abundance.
- All ten simulation fits and independent endpoint evaluations, and all five empirical fits and point-endpoint evaluations, completed successfully. Their completion records, configuration selections, training-index assertions and aggregate outputs are in `manuscript/generated/review_validation/euler_completed/`. The integration script passes 225 checks, including exact simulation input hashes and reconstruction of the previous outputs.
- Current simulation results are in `results/final/tables/` and manuscript Table S4. The corrected empirical point endpoints are in Table S26 and a complete 84-row CSV. Earlier-adapter empirical confidence intervals are not reused, so these corrected points are separate from empirical figures requiring paired uncertainty.
- The publication bundle includes identity retrieval, null-gate calibration, rare/coordinated-response stress experiments, ANCOVA/mean-only comparisons, design audits and their code. `manuscript/generated/review_validation/README.md` describes their estimands and provenance.
- All main and supplementary figure panels, numerical source tables and assembled PDF bundles are included under `manuscript/figures/`. The matching sources are `manuscript/main.tex` and `manuscript/supplementary.tex`; a subsequent documentation commit adds their immutable analysis citation.

`publication_manifest.csv` records the byte count and SHA-256 of every public repository file other than itself. The original task manifests and aggregation code under `analysis/` remain available as historical execution infrastructure. Current revision commands below are the publication workflow; the original aggregator does not reproduce the added analyses or revised empirical uncertainty.

## Environment and public inputs

Restore the recorded R environment and install the frozen benchmark package:

```sh
Rscript -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv")'
Rscript -e 'renv::restore(prompt=FALSE)'
mkdir -p Rlib
R CMD INSTALL --library=Rlib package/winn_0.1.4.tar.gz
export R_LIBS_USER="$PWD/Rlib"
```

The production RFSC reruns require statTarget 1.38.0 and refuse other versions unless an explicit smoke-test flag is supplied. The frozen package SHA-256 and environment records are in `package/`; the adapter regression test was also checked locally with statTarget 1.40.0. Figure generation requires the R packages listed by the scripts and Python with `pypdf`; Python 3.10 or later is recommended for the figure utilities.

The source datasets and preparation decisions are documented in [`docs/DATASETS.md`](docs/DATASETS.md). Download and preprocess the four public empirical datasets and generate the ten canonical simulation seeds:

```sh
make download
make preprocess
make simulations
Rscript analysis/prepare_revision_inputs.R --dataset=public
```

The last command creates `results/revision_inputs/input_<dataset>.rds` from the public prepared matrices, authorized training references, frozen internal folds and candidate registry. SIM01 and all four empirical matrices prepared through this path were verified against the frozen exact object hashes. The generated bundles are ignored by Git. Clinical inputs and their identity-bearing fold manifest require authorized local files and are never included in this repository.

## Reproduce the revised analyses

Run the following from the repository root. Full RFSC execution is compute intensive; each task has independent outputs and a completion record.

```sh
# Test the actual public wrapper with the real statTarget training-index logic.
Rscript analysis/test_qcrfsc_public_wrapper.R

# One empirical dataset; repeat for mtbls79, batchcorr_set1, sacurine and waveica.
Rscript manuscript/scripts/benchmark_revision/run_qcrfsc_revision.R \
  --dataset=sacurine --input-root=results/revision_inputs --output-root=results/revised_rfsc
Rscript manuscript/scripts/benchmark_revision/evaluate_empirical_qcrfsc_revision.R \
  --dataset=sacurine --input-root=results/revision_inputs \
  --matrix-root=results/revised_rfsc --output-root=results/revised_rfsc_evaluation

# One canonical seed; repeat for SIM01 through SIM10.
Rscript manuscript/scripts/benchmark_revision/run_simulation_qcrfsc_revision.R \
  --seed-id=SIM01 --input-root=results/revision_inputs --output-root=results/revised_simulation_rfsc
Rscript manuscript/scripts/benchmark_revision/evaluate_simulation_qcrfsc_revision.R \
  --seed-id=SIM01 --input-root=results/revision_inputs \
  --matrix-root=results/revised_simulation_rfsc --output-root=results/revised_simulation_evaluation

# Gate and biological-response stress experiments use the frozen source snapshot.
Rscript manuscript/scripts/targeted_validation/run_targeted_biology_audit.R
Rscript manuscript/scripts/targeted_validation/verify_null_gate_calibration.R
```

The regression baseline and biological endpoint commands are documented in [`manuscript/scripts/nearest_baseline/README.md`](manuscript/scripts/nearest_baseline/README.md). The RFSC revision protocol and integration checks are documented in [`manuscript/scripts/benchmark_revision/README.md`](manuscript/scripts/benchmark_revision/README.md). Fresh evaluations should be written into separate output directories; do not substitute new matrix results into existing uncertainty intervals.

For Euler, connect with `ssh euler-ext`, run `source ~/.bash_profile` before Git operations, and load the recorded R environment. The revision `.sbatch` files document the actual job commands and resource requests. Their project paths are execution provenance and must be adapted to another installation.

## Reproduce the supplementary tables

The complete minimal historical inputs for S4–S20 are supplied as 15 compressed aggregate CSVs in `results/final/revision_inputs/supplement_20260727/`. The original public base did not include the derivative portion of this export. No individual clinical data are included. Regeneration applies the validated RFSC replacements and then the current retrieval, baseline and support-table outputs:

```sh
make tables
```

The exact commands and S4–S31 source mappings are in [`manuscript/generated/supplementary_table_data/README.md`](manuscript/generated/supplementary_table_data/README.md). An isolated clean-input regeneration passed all 46 comparisons: 18 parsed CSV tables and exact bytes for all 28 TeX tables. The archived source files are historical inputs; current publication values are in the manuscript outputs.

## Verify and render the publication bundle

```sh
make validate-results
make report
```

Validation uses Python's standard library and checks all file hashes, simulation means, corrected RFSC endpoints, selected settings, aggregate clinical-data fields and required source/figure assets. `make report` regenerates current main and supplementary figure panels and assembles the vector PDF bundles from their committed numerical inputs; it does not refit correction methods. Re-rendering figures can change PDF metadata and file hashes, so validate the committed bundle before rendering and compare the numerical sources separately after intentional regeneration.

Compile `manuscript/main.tex` and `manuscript/supplementary.tex` in a LaTeX environment supporting the supplied Bioinformatics template. The main paper and supplement use relative paths to the included assets.

## Data availability

MTBLS79, BatchCorrMetabolomics Set 1, Sacurine and WaveICA are obtained from the public sources in [`docs/DATASETS.md`](docs/DATASETS.md). The internal clinical cohort data will be released with a forthcoming manuscript describing that study. This repository currently distributes aggregate clinical results and generic analysis code; it contains no clinical matrices, biological identifiers, identity-level summaries or private split assignments.

## Citation and license

Please cite the WiNN manuscript and the original publication for each dataset used. Use `citation("winn")` for the package citation. Analysis code is distributed under [LICENSE](LICENSE); vendored code and source datasets retain their stated licenses. The unchanged Wehrens implementation and frozen WiNN source snapshot include their license notices.
