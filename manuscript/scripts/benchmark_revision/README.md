# Comparator revision audit — 27 October 2026

The professor's comments are review evidence, not executable instructions. Findings 5–8 were assessed against the frozen analysis code, manifests and available numerical outputs. The demonstrated ComBat reporting error, injected-effect scale error and QC-RFSC adapter defects were accepted. The clinical bootstrap displacement is real in the saved replicates and is not a plotting-index error. Neither the review nor this audit establishes that WiNN collapses biological profiles.

## QC-RFSC protocol change

`qcrfsc_review_adapter.R` neutralizes **every** sample ID, appends `_QC` only for authorized training references, retains a reversible internal-to-original mapping, and checks the exact statTarget name predicate. The rerun driver also traces the real `REGfit` function and asserts the exact authorized training set on every call. The original wrapper exposed 12 MTBLS79, eight Sacurine and 15 WaveICA hidden QCs through their identifiers; no identifier leakage was found in BatchCorr or clinical data.

statTarget's native ratio output is divided by 1000 and multiplied by the feature's median **raw training-QC abundance**. This reference uses no hidden references, biological labels or ground truth. Raw fallback entries therefore remain in abundance units; fallback feature/sample counts are retained. This is a declared revised analysis policy, not a silent change to the frozen release. The corrected simulation truth endpoints use the restored abundance matrix and the same log1p evaluation functions as the other comparators.

The original four-candidate QC-RFSC grids, internal training-QC folds and tuning feature panels were retained. Candidate ranking uses mean training-validation profile correlation minus half mean training-validation CV, averaged across the frozen internal folds. The training-reference abundance policy can change the profile-correlation component and hence the selected candidate. MTBLS79 uses its documented in-sample training-QC sensitivity mode; other datasets use their existing frozen modes. All selection and parameter metadata are saved before held-out diagnostics.

Run from any directory with authorized prepared inputs:

```sh
Rscript run_qcrfsc_revision.R --dataset=mtbls79 --input-root=/path/to/prepared_inputs --output-root=/path/to/new_results
Rscript test_qcrfsc_review_adapter.R
Rscript diagnose_mtbls79_qcrlsc.R --input-root=/path/to/prepared_inputs --output-root=/path/to/new_diagnostics
```

Production QC-RFSC reruns require the original statTarget 1.38.0. The integration test uses the actual installed package and records its version. It checks adversarial QC-bearing hidden IDs, invariance to renaming ordinary samples, dimension/ID preservation and training-QC-only abundance restoration. The same integration test passed locally with statTarget 1.40.0 and on Euler with the frozen statTarget 1.38.0 (job **18543328**, three inside-REGfit assertions). A separate Euler smoke test passed with 1.38.0, and all production tasks assert authorized controls inside 1.38.0's `REGfit`.

Inputs are deliberately excluded. Each `input_<dataset>.rds` contains `x`, sample-aligned `meta`, `training_ids`, `hidden_ids`, `external_ids`, the frozen `folds` list, `tuning_feature_ids` and `registry`. These were prepared using the existing `tq_prepare_dataset`, `tq_fold_definitions`, `tq_tuning_features` and `tq_read_registry` helpers. Clinical inputs and reversible clinical ID maps remain restricted. The four public input matrix SHA256 values match the frozen release when CSV-derived inputs are prepared under R 4.5. The clinical log-scale pair reconstruction differs at floating-point precision and is explicitly not byte-identical to the frozen input; it is a reconstruction diagnostic.

Euler project directory:

```text
/cluster/project/platt/tanmay/metabolomics_winn_robustness/review_revision_20261027
```

Pinned-version smoke job: **18541777**, completed successfully. Full corrected QC-RFSC jobs: **18541987_1–5**, respectively MTBLS79, clinical, BatchCorr, Sacurine and WaveICA. The compact status and diagnostic outputs are in `generated/review_validation/qcrfsc_revision_*.csv`. All five correction jobs completed successfully. Their held-out reference CV and fallback counts are published in Table S26. Previous-adapter estimates have been removed from the manuscript tables and figures; the simulation values are replaced by the completed corrected analysis.

## Other accepted findings and limits

- The ComBat manifest and executor agree on training-QC selection from parametric/nonparametric candidates for MTBLS79, BatchCorr and clinical data; clinical uses `par_prior=FALSE`. Full frozen matrices were not all available for direct matrix-by-matrix reconstruction; output hashes and reuse provenance are available.
- The injected effect is on log intensity while the fitted coefficient is on log1p intensity. Deterministic symmetric-distribution calculation gives a clean contrast of 0.175 for an injected log contrast of 0.35 at log baseline zero. The ratio 0.5 is a transformation effect, not failed recovery. The clean-reference denominator remains appropriate.
- Every clinical batch-R2 bootstrap percentile bound was independently recomputed from the saved 2000 replicates and matched. WiNN's point 0.0020543, bootstrap mean 0.0141233 and interval [0.0073132, 0.0232605] describe a displaced conditional percentile distribution. Do not force its interval to include the point or treat it as full-pipeline uncertainty.
- The exact qcrlscR 0.1.3 implementation reproduces MTBLS79's 182.843579% held-out CV (raw 28.757223%). It emitted 266,216 warnings, including near-singular local fits and too-small spans. Six of eight batches have only three training QCs, and four batches extrapolate over four or five early injections. The selected span is 0.7, degree one. Coverage, warnings and representative failing profiles are supplied. This supports a configuration/layout-specific numerical-instability diagnosis, not failure of the entire QC-RLSC method class.
- Dataset-level primary intensity flooring and TIGER feature-restoration records are supplied. Retaining a feature does not imply that it was successfully corrected.

The nearest published regression baseline is implemented and documented separately in `scripts/nearest_baseline`. It uses the authors' unchanged Wehrens implementation, records additional information used by the study-only sensitivity, and reports zero-floor counts. These descriptive additions do not establish biological effect preservation or a new overall method ranking.

## Ten-seed truth-endpoint rerun

The existing canonical SIM01–SIM10 raw and clean RDS matrices were recovered locally; all 20 objects match the frozen exact SHA256 and nine-decimal portable hashes. No generator was substituted. `simulation_revision_input_identity.csv` records this preflight. Each new task verifies its input identity again, then repeats the original prespecified QC-RFSC settings (`ntree=500`, `coCV=30`, `Frule=0.8`) and original release method RNG. Both the legacy and corrected adapter are fitted in an isolated output directory. The driver recomputes held-out CV and the same three truth endpoints: log1p RMSE, feature-profile ICC(A,1), and sample-profile Pearson. It saves matrices, exact training-index assertions, metrics and a completion marker. The original release remains unchanged. Validated corrected results update the manuscript simulation tables and panels through `integrate_completed_results.py`. The independent full-endpoint evaluation (`evaluate_simulation_qcrfsc_revision.R`, job 18547170_1–10) completed and reproduced the four initially evaluated endpoints while adding residual GAM deviance, batch association and remaining technical summaries.

```sh
Rscript run_simulation_qcrfsc_revision.R --seed-id=SIM01 --input-root=/path/to/exact_inputs --output-root=/path/to/new_simulation_results
```

Each `input_SIMxx.rds` contains `x`, `truth`, sample-aligned `meta`, authorized `training_ids`, `hidden_ids` and `external_ids`; inputs remain outside the manuscript.

## Empirical endpoint reevaluation

`evaluate_empirical_qcrfsc_revision.R` evaluates the legacy-reproduced and corrected QC-RFSC matrices using unchanged frozen evaluation functions: reference CV, residual GAM drift, residual Ljung–Box tests, categorical weighted-PC batch association, and dataset-appropriate biological/replicate endpoints. It has been smoke-tested on the full Sacurine matrix. The Euler evaluation array depends on all five corrected fits succeeding. Its row-level feature and identity records remain on Euler; only aggregate metrics may be exported to the manuscript. All five evaluations completed. Technical points and comparable replicate-profile endpoints are integrated in Table S26. Generic study-injection biological summaries remain in the full point export and are not substituted for the manuscript's identity-collapsed association estimands. Paired uncertainty was not rerun, so RFSC is not included in the empirical figures that require those intervals.

```sh
Rscript evaluate_empirical_qcrfsc_revision.R --dataset=sacurine --input-root=/path/to/inputs --matrix-root=/path/to/corrected_results --output-root=/path/to/endpoint_results
```

## Current integration status

At **2026-10-27T15:15:46Z**, all five empirical fit jobs **18541987_1–5**, all ten simulation fit jobs **18543837_1–10**, and all ten full simulation endpoint jobs **18547170_1–10** had completed successfully. The original raw and previous-adapter simulation outputs agree with the archived results to within 1e-11. Twenty raw/clean input hashes match exactly and every real REGfit call uses only authorized references.

`integrate_completed_results.py` validates these outputs and updates Figure 2, Figure S3, Table S4, selected parameters and Table S26. It also removes previous-adapter RFSC rows from empirical uncertainty, association, identity and reference-split tables; those intervals are not assigned to the new matrices. The manuscript reports corrected simulation comparisons and all five completed empirical point evaluations from array **18544176**. Table S26 reports reference precision, residual GAM deviance and categorical batch association, with replicate ICC/Pearson for BatchCorr and clinical data. The complete 84-row endpoint export is `qcrfsc_completed_empirical_all_endpoints.csv`. These point estimates have no new bootstrap or permutation intervals; the earlier adapter's intervals are not assigned to the new matrices.

Run the integration step from the manuscript directory, then regenerate the main and simulation panels:

```sh
python3 scripts/benchmark_revision/integrate_completed_results.py
Rscript scripts/generate_revision_main_panels.R
```

Validation, selected configurations, session information and aggregate results are retained under `generated/review_validation`; input matrices and clinical identity maps remain on Euler.
