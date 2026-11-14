# WiNN manuscript and reproducible analyses

This repository contains the code, numerical results, figure source data, compiled figures and matching main/supplementary LaTeX sources for the manuscript *WiNN enables selective correction of run-order drift and batch effects in metabolomics data*. All results use **WiNN 0.1.6**, which is available from [CRAN](https://CRAN.R-project.org/package=winn) (`install.packages("winn")`) and is bundled here as `package/winn_0.1.6.tar.gz`. The manuscript cites code revision [`853541434d507e9cfe0a6a94c205b0e7d1b39f29`](https://github.com/ratschlab/WiNN-for-technical-correction/tree/853541434d507e9cfe0a6a94c205b0e7d1b39f29), which includes the author-note reproducibility fixes and QC-RFSC validation. Frozen benchmark results originate from revision [`491a6f9027f4a7ad7d6d17fc11d3b986bd5b0b13`](https://github.com/ratschlab/WiNN-for-technical-correction/tree/491a6f9027f4a7ad7d6d17fc11d3b986bd5b0b13); the additional MTBLS79 QC-RFSC figure point is evaluated from the already corrected archived matrix on the same 20 sample means as the other methods. The publication bundle is in [`manuscript/`](manuscript/); its numerical results and complete file hashes are checked with `make validate-results`.

The primary configuration is called **WiNN**. The supplementary configurations, **WiNN (QC-tuned)** and **WiNN (inferred segments)**, use QC-based tuning or inferred acquisition segments and are identified where a comparison requires that distinction. Stable internal method identifiers (`winn_fixed_no_qc`, `winn_auto_qc`, `winn_auto_batch_qc`) are preserved so that analysis records remain joinable.

## WiNN 0.1.6 results

WiNN 0.1.6 corrects the smoother calls in drift correction. In versions 0.1.4 and 0.1.5, `mgcv::s()` inside `mgcv::gam()` formulas raised an error that was caught internally, so every spline fit fell back to LOESS (conservative smoother) or a linear fit (standard smoother). Version 0.1.6 fits the documented P-spline and cubic-regression-spline GAMs. The package records are in `package/frozen_package.json` and `package/VALIDATION.md`.

Only WiNN depends on this change, so only WiNN was recomputed:

- **Recomputed:** the three WiNN configurations, together with uncorrected Raw as a reproduction check, in 480 Euler tasks: 24 primary fits, 40 simulation fits, 240 reference-split fits, 6 ablation bundles, 10 simulation-ablation bundles and 160 partial-confounding fits. Evaluation and aggregation were restricted to these methods (`WINN_RERUN_METHODS`; `analysis/euler/submit_winn_rerun.sh`, `analysis/euler/aggregate_winn_rerun.sbatch`).
- **Reused unchanged:** all comparator outputs (ComBat, QC-RLSC, QC-RFSC, TIGER and SERRF) from the validated previous release. `analysis/splice_winn_rerun.R` copies those rows and replaces every WiNN row; the re-run Raw rows reproduced the previous release to within 5e-14. Comparator and WiNN outputs therefore share the same input matrices, seeds, reference partitions and resampling indices, so paired comparisons remain exact.
- **Derivative analyses:** bootstrap, permutation and leave-one-batch-out summaries were recomputed for the WiNN configurations only, with the same seeds and resampling indices; all other rows are reused.
- **Validation:** the release validation (19 checks) combines the WiNN-only re-run with the validation of the reused comparator outputs (`analysis/combine_winn_rerun_validation.R`); all 540 reference-separation checks, 225 QC-RFSC integration checks and 91 baseline and identity checks pass, and all 56 supplementary tables regenerate exactly from the bundled inputs (`analysis/config/revision_20261111/`).

The exact aggregate inputs for the current tables are in `results/final/revision_inputs/supplement_20261111/` and `results/final/revision_inputs/derivative_20261111/`; each directory has a README and an uncompressed SHA-256 manifest. The WiNN 0.1.4 release remains available at commit [`f59db42`](https://github.com/ratschlab/WiNN-for-technical-correction/tree/f59db4273bc65d2f09d4cd7705a1b9cbcaf84a3c) (analysis revision [`97e5ea3`](https://github.com/ratschlab/WiNN-for-technical-correction/tree/97e5ea3cadb44fe74ee914528ff19a2940fe95d7)), and the original benchmark at [`21de7c4`](https://github.com/ratschlab/WiNN-for-technical-correction/tree/21de7c448bfa5545d27c5a930161fb0bd8e40d85); `package/winn_0.1.4.tar.gz` is kept for that provenance.

## Author-note review (14 November 2026)

The review re-executed the null-gate and biological-preparation experiments;
all 14,877 rows in their nine archived result CSVs matched exactly. The package
test suite passed 105 assertions, and the public QC-RFSC wrapper integration
test passed with the frozen statTarget 1.38.0 installation on Euler. The partial-confounding preflight now checks
the authoritative log-matrix hashes with exact field names; diagnostic intensity
hashes are no longer used as cross-platform gates. `make test-confounding-inputs`
reconstructs and checks all 160 publication scenarios without fitting a method.
The older 320-scenario construction configuration is retained as provenance;
the publication task manifest selects ten seeds and 16 scenarios per seed.
The saved corrected QC-RFSC matrix also supplies the previously missing
MTBLS79 biological point on the same 20 sample means as the other methods;
`analysis/evaluate_mtbls79_rfsc_sample_means.R` reproduces this point without
refitting correction or uncertainty. Scope, findings and limits are recorded in
[`analysis/config/revision_20261114/code_review.txt`](analysis/config/revision_20261114/code_review.txt).

## QC-RFSC and revision analyses

- The QC-RFSC wrapper in `scripts/benchmark_helpers.R` assigns neutral identifiers to every injection and marks only authorized training references as QCs. Returned relative abundances are restored using raw training-QC feature medians. Features not returned by the native implementation retain their input abundance.
- All ten simulation fits and all five empirical fits of the corrected QC-RFSC adapter completed. Their completion records, configuration selections, training-index assertions and aggregate outputs are in `manuscript/generated/review_validation/euler_completed/`.
- Current simulation results are in `results/final/tables/` and Supplementary Table S4. Corrected empirical QC-RFSC values are point estimates (Supplementary Table S4 and a complete 84-row CSV); earlier-adapter confidence intervals are not reused.
- The publication bundle includes identity retrieval, null-gate calibration, rare- and coordinated-response stress experiments, ANCOVA and mean-only baselines, design audits and their code. `manuscript/generated/review_validation/README.md` describes their estimands and provenance.
- All main and supplementary figure panels, numerical source tables and assembled PDF bundles are under `manuscript/figures/`, with the matching sources `manuscript/main.tex` and `manuscript/supplementary.tex`.

The 14 November manuscript revision resolves all author notes, applies the same terminology to figures and supplementary methods, and includes numerical and layout verification under `manuscript/revision_notes/`. `bash manuscript/scripts/build_manuscript.sh` builds both documents; the final layout audit records seven main-text and 47 supplementary pages, with no overfull boxes or text outside page boundaries. Numerical claim and figure checks are retained with the manuscript.

`publication_manifest.csv` records the byte count and SHA-256 of every public repository file other than itself. The original task manifests and aggregation code under `analysis/` remain available as execution infrastructure; the commands below are the publication workflow.

Selection-test limits and the separate permutation-calibration comparison are confined to Supplementary Methods S10 and Figure S11. The experimental alternative was not applied to any correction benchmark.

The focused Figure 2E–F confounding revision is documented in [`analysis/figure2_confounding/README.md`](analysis/figure2_confounding/README.md). It changes phenotype allocation and biological-preservation endpoints only; WiNN 0.1.6, the technical-artifact generator, Figure 2A–D and all other benchmark outputs remain fixed.

## Environment and public inputs

Restore the recorded R environment and install WiNN 0.1.6, either from CRAN or from the bundled source archive:

```sh
Rscript -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv")'
Rscript -e 'renv::restore(prompt=FALSE)'
make install-winn        # installs package/winn_0.1.6.tar.gz into ./Rlib
export R_LIBS_USER="$PWD/Rlib"
```

Production tasks check the installed WiNN version, installation path and source-archive SHA-256 before loading data (`package/assert_frozen_winn.R`). The production RFSC reruns require statTarget 1.38.0 and refuse other versions unless an explicit smoke-test flag is supplied. Figure generation requires the R packages listed by the scripts and Python with `pypdf`; Python 3.10 or later is recommended for the figure utilities.

The source datasets and preparation decisions are documented in [`docs/DATASETS.md`](docs/DATASETS.md). Download and preprocess the four public empirical datasets and generate the ten canonical simulation seeds:

```sh
make download
make preprocess
make simulations
Rscript analysis/prepare_revision_inputs.R --dataset=public
```

The last command creates `results/revision_inputs/input_<dataset>.rds` from the public prepared matrices, authorized training references, frozen internal folds and candidate registry. SIM01 and all four empirical matrices prepared through this path were verified against the frozen exact object hashes. The generated bundles are ignored by Git. Clinical inputs and their identity-bearing fold manifest require authorized local files and are never included in this repository.

## Reproduce the analyses

**WiNN fits.** On a Slurm cluster, prepare a release directory with the frozen package installed in `$WINN_RELEASE_ROOT/Rlib`, then submit the WiNN-only re-run, which also runs the evaluations and aggregation:

```sh
export WINN_RELEASE_ROOT=/path/to/release
bash analysis/euler/submit_winn_rerun.sh
# After the aggregate job: splice the re-run into the previous release and re-render overview figures.
Rscript analysis/splice_winn_rerun.R --new=<rerun results/final> --old=<previous results/final> \
  --overlay=<corrected-RFSC results/final> --out=<spliced results/final>
Rscript analysis/combine_winn_rerun_validation.R --rerun=<rerun validation csv> \
  --previous=<previous validation csv> --out=<combined validation csv>
```

The `.sbatch` files document the job commands and resource requests used on Euler (`ssh euler-ext`; load the recorded R environment). Their project paths are execution provenance and must be adapted to another installation. The complete nine-method benchmark uses the task manifests and `.sbatch` files under `analysis/euler/`, with separate TIGER jobs because of their runtime.

**QC-RFSC revision and stress tests.** Run from the repository root; each task has independent outputs and a completion record.

```sh
# Test the public wrapper with the real statTarget training-index logic.
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

# Gate and biological-response stress experiments use the frozen WiNN 0.1.6 source snapshot.
Rscript manuscript/scripts/targeted_validation/run_targeted_biology_audit.R
Rscript manuscript/scripts/targeted_validation/verify_null_gate_calibration.R
```

The regression baseline and biological endpoint commands are documented in [`manuscript/scripts/nearest_baseline/README.md`](manuscript/scripts/nearest_baseline/README.md), and the RFSC revision protocol and integration checks in [`manuscript/scripts/benchmark_revision/README.md`](manuscript/scripts/benchmark_revision/README.md). The baseline scripts accept `--reuse-baselines-from` and `--reuse-baseline-rows` to reuse WiNN-independent baseline rows. Fresh evaluations should be written into separate output directories; do not substitute new matrix results into existing uncertainty intervals.

## Reproduce the supplementary tables

The aggregate inputs for all supplementary tables are supplied as compressed CSVs in `results/final/revision_inputs/supplement_20261111/`. No individual clinical data are included. Regeneration applies the validated QC-RFSC replacements, writes the retrieval, baseline and support tables, and assembles the 13 printed tables:

```sh
make tables
```

The printed supplement condenses 31 detailed tables into 13 (`manuscript/scripts/condense_supplementary_tables.py`); the detailed tables remain as source data. The mapping between published and detailed tables and the exact commands are in [`manuscript/generated/supplementary_table_data/README.md`](manuscript/generated/supplementary_table_data/README.md). An isolated regeneration from these inputs reproduced all 32 TeX tables byte for byte and all 24 CSV tables row for row (`analysis/config/revision_20261111/table_reproduction_checks.csv`).

## Verify and render the publication bundle

```sh
make validate-results
make report
```

Validation uses Python's standard library and checks all file hashes, simulation means, corrected RFSC endpoints, selected settings, table regeneration records, aggregate clinical-data fields and required source and figure assets. `make report` regenerates the main and supplementary figure panels and assembles the vector PDF bundles from their committed numerical inputs; it does not refit correction methods. Re-rendering figures can change PDF metadata and file hashes, so validate the committed bundle before rendering and compare the numerical sources separately after intentional regeneration.

Compile `manuscript/main.tex` and `manuscript/supplementary.tex` in a LaTeX environment supporting the supplied Bioinformatics template. Both use relative paths to the included assets.

## Data availability

MTBLS79, BatchCorrMetabolomics Set 1, Sacurine and WaveICA are obtained from the public sources in [`docs/DATASETS.md`](docs/DATASETS.md). The internal clinical cohort data will be released with a forthcoming manuscript describing that study. This repository distributes aggregate clinical results and generic analysis code only; it contains no clinical matrices, biological identifiers, identity-level summaries or private split assignments.

## Citation and license

Please cite the WiNN manuscript and the original publication for each dataset used. Use `citation("winn")` for the package citation. Analysis code is distributed under [LICENSE](LICENSE); vendored code and source datasets retain their stated licenses. The unchanged Wehrens implementation and the frozen WiNN source snapshot include their license notices.
