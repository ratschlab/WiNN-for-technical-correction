# QC-RFSC implementation and uncertainty provenance

## Current chronological implementation — 22 November 2026

The primary adapter preserves chronological injection rank for every candidate-fold and final fit. A code-level review established that statTarget 1.38.0's `REGfit` fits feature abundance against training-QC **column positions**, so moving a reference to a sequence boundary changes its predictor. The former boundary-QC wrapper changed 11 clinical coordinates (including seven study units; maximum displacement eight ranks). The acquisition-geometry audit found no such movement in any full or tuning fit in the four public datasets.

`qcrfsc_review_adapter.R` makes isolated copies of the installed `shiftCor` and `REGfit` functions in a private environment. It removes exactly four statements: the two endpoint-QC index assignments and the two corresponding stop checks. The native forest fitting and prediction body is unchanged. Forests predict observations outside the training-QC range using terminal-leaf endpoint values; no reference is relocated or fabricated. Checks inside every copied `REGfit` require the full chronological column sequence, `y=1:n`, and exactly the authorized masked training references. No installed namespace is modified. `test_qcrfsc_review_adapter.R` confirms exact numerical agreement with native statTarget when boundary QCs are present, valid chronological fitting with both boundary QCs withheld, invariance to adversarial hidden-QC names, training-reference abundance restoration and unchanged native function hashes.

Native statTarget preprocessing treats nonpositive values as missing, applies `FilterMV` and performs KNN imputation using the full input matrix. Held-out references are excluded from forest training, training-QC abundance scaling and parameter selection; they are not excluded from this native preprocessing. The aggregate input audit found 48,820 zero entries in WaveICA (2,742 in training QCs, 587 in held-out references and 45,491 in study injections), with no negative or preexisting missing values. Native `FilterMV` removes nine WaveICA features, leaving 47,912 missing entries to be imputed (2,635 in training QCs, 573 in held-out references and 44,704 in study injections). The other four empirical datasets and all ten simulation inputs contain no zero, negative or nonfinite entries, and have no missing entries after native filtering. `qcrfsc_preprocessing_input_audit.csv` records these counts and input hashes without sample identifiers; `scripts/benchmark_revision/audit_qcrfsc_preprocessing_inputs.R` reproduces the audit.

The five empirical datasets are refitted using the same four-candidate grids, tuning feature panels, internal folds and per-fit seeds. Independent candidate-fold fits run in separate forked processes with explicit original seeds; the final fit also resets its original seed. Full public fits must reproduce the previous numeric matrices exactly, checked independently; the clinical fit uses the corrected coordinates. Ten simulation fits retain their prespecified parameters, original raw/clean matrices and seeds. Current production commands are:

```sh
sbatch scripts/benchmark_revision/run_chronological_empirical_qcrfsc.sbatch
sbatch --dependency=aftercorr:<empirical-array-id> scripts/benchmark_revision/run_chronological_empirical_followup.sbatch
sbatch scripts/benchmark_revision/run_simulation_qcrfsc_chronology.sbatch
sbatch --dependency=aftercorr:<simulation-array-id> scripts/benchmark_revision/evaluate_simulation_qcrfsc_chronology.sbatch
```

The empirical follow-up regenerates corrected-matrix point endpoints and all 2,000-draw conditional bootstrap summaries. Integration is separated by dataset type and does not erase other methods or previously completed uncertainty:

```sh
python3 scripts/benchmark_revision/integrate_chronological_empirical.py
python3 scripts/benchmark_revision/integrate_chronology_simulation.py
python3 scripts/benchmark_revision/integrate_qcrfsc_bootstrap.py
python3 scripts/benchmark_revision/validate_qcrfsc_bootstrap_pairing.py --archive-root=/path/to/frozen/derivative
```

Point-estimate data, parameter manifests, matrices' hashes, chronology/authorized-reference assertions, native function bodies, adapter tests, geometry audit, exact public-matrix comparisons and full aggregate bootstrap series are retained under `generated/review_validation`. The complete raw 2,000-draw replay reproduces all archived empirical PC and replicate series (maximum absolute discrepancy below 5e-12). Each raw and QC-RFSC draw list has the same SHA256, independently verifying pairing beyond shared seed numbers. Query-weighted retrieval draws additionally reproduce the current matching archive.

The residual Ljung–Box summary now explicitly excludes unavailable tests and returns NA if none are valid. All five preceding empirical RFSC evaluations had zero unavailable adjusted tests, so this bounded evaluator correction leaves those reported fractions unchanged; the aggregate audit is supplied. Conditional bootstrap uncertainty fixes the fitted correction and parameter choice, as for the other methods. Reference-split refitting, association-count permutation and leave-one-batch-out analyses remain separate sensitivity analyses restricted to their named eight outputs; no old RFSC intervals are assigned to new matrices.

The superseded pre-chronology uncertainty outputs are retained only under `generated/review_validation/archive_boundary_relocation_20261122/`, explicitly marked as audit artifacts. Current publication inputs are the chronology-safe files generated by the integration chain above.

## Historical audit record


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

`evaluate_empirical_qcrfsc_revision.R` evaluates the legacy-reproduced and corrected QC-RFSC matrices using unchanged frozen evaluation functions: reference CV, residual GAM drift, residual Ljung–Box tests, categorical weighted-PC batch association, and dataset-appropriate biological/replicate endpoints. It has been smoke-tested on the full Sacurine matrix. The Euler evaluation array depends on all five corrected fits succeeding. Its row-level feature and identity records remain on Euler; only aggregate metrics may be exported to the manuscript. All five evaluations completed. Technical points and comparable replicate-profile endpoints are integrated in Table S26. Generic study-injection biological summaries remain in the full point export and are not substituted for the manuscript's identity-collapsed association estimands. Conditional empirical uncertainty is now complete: the 22 November 2026 analysis below evaluates the corrected matrices with the same 2,000 paired biological-unit draws as the other methods.

```sh
Rscript evaluate_empirical_qcrfsc_revision.R --dataset=sacurine --input-root=/path/to/inputs --matrix-root=/path/to/corrected_results --output-root=/path/to/endpoint_results
```

## Current integration status

At **2026-10-27T15:15:46Z**, all five empirical fit jobs **18541987_1–5**, all ten simulation fit jobs **18543837_1–10**, and all ten full simulation endpoint jobs **18547170_1–10** had completed successfully. The original raw and previous-adapter simulation outputs agree with the archived results to within 1e-11. Twenty raw/clean input hashes match exactly and every real REGfit call uses only authorized references.

`integrate_completed_results.py` validates these outputs and updates Figure 2, Figure S3, Table S4, selected parameters and Table S26. It also removes previous-adapter RFSC rows from empirical uncertainty, association, identity and reference-split tables; those intervals are not assigned to the new matrices. The manuscript reports corrected simulation comparisons and all five completed empirical point evaluations from array **18544176**. Table S26 reports reference precision, residual GAM deviance and categorical batch association, with replicate ICC/Pearson for BatchCorr and clinical data. The complete 84-row endpoint export is `qcrfsc_completed_empirical_all_endpoints.csv`. The current conditional bootstrap intervals are computed from these corrected matrices as documented below. Earlier-adapter intervals are never assigned to new matrices. Association-count permutations, leave-one-batch-out sensitivity and reference-split refitting remain separate analyses restricted to their explicitly listed outputs.

Run the integration step from the manuscript directory, then regenerate the main and simulation panels:

```sh
python3 scripts/benchmark_revision/integrate_completed_results.py
Rscript scripts/generate_revision_main_panels.R
```

Validation, selected configurations, session information and aggregate results are retained under `generated/review_validation`; input matrices and clinical identity maps remain on Euler.


## Completed conditional empirical uncertainty — 22 November 2026

`bootstrap_empirical_qcrfsc.R` reads the five completed `corrected_training_qc_abundance.rds` matrices produced with statTarget 1.38.0. It does not refit correction, change selected parameters, or repeat tuning. The same frozen derivative evaluator, seed (`2026072701` plus dataset offsets 0/79/113/199/1594), biological-unit order and 2,000 draws reproduce the manuscript's conditional bootstrap protocol. Each draw recomputes PCA after feature standardization. WaveICA resamples within batch-by-group strata; MTBLS79 resamples extracts within species, carrying injections together and averaging log1p profiles for the species endpoint; Sacurine resamples within batch; BatchCorr resamples plant lines and clinical FIA-MS resamples recorded biological-sample identities.

The endpoint export includes batch association in all five datasets, biological association in WaveICA and MTBLS79, age/BMI/gender association in Sacurine, and replicate Pearson/ICC plus query-weighted sample matching in BatchCorr and the clinical cohort. `complete_qcrfsc_benchmark_points.R` additionally supplies the existing table's point-only MTBLS79 feature ICC and WaveICA ordinary-model association count/standardized effect. Its raw results reproduce the frozen table source. The generic full endpoint export's limma-based WaveICA count and injection-based MTBLS79 biological-PC summary are not substituted for the table's ordinary-model count and extract-mean PC estimand. Matching keeps candidate pools and fitted matrices fixed and uses each sampled identity's eligible-query count, consistent with the published point estimand. Aggregate replicates permit paired WiNN-minus-QC-RFSC differences without exporting sample identities. Standard errors are the standard deviations of the 2,000 bootstrap estimates; percentile bounds use R's default type-7 quantiles. Technical QC-CV and drift bars remain descriptive feature/profile SEMs, calculated from the corrected-matrix SD and denominator exports.

The identifier and scale corrections are unchanged: all sample names are neutralized; only training references receive the `_QC` tag; each actual `REGfit` call asserts the exact authorized reference set. Native ratio values are divided by 1,000 and multiplied by each feature's raw training-QC median. This restores abundance units using training references alone and leaves fallback entries in the same units. The 1.38.0 integration test and inside-REGfit assertions passed; fallback counts are retained in the source tables. These counts identify features omitted from statTarget output (including its feature filtering), not a proven count of nonfinite fits; all final matrices are finite. The selected four-candidate training-QC configurations remain those in `generated/supplementary_table_data/endpoint_free_selection_manifest.csv` and each corrected fit's `selected_parameters.csv`: `ntree=500, coCV=30, Frule=0.8` in MTBLS79, clinical, Sacurine and WaveICA, and `ntree=200, coCV=30, Frule=0.8` in BatchCorr. No new parameter selection was performed for uncertainty. `scripts/benchmark_revision/actual_comparator_parameters.csv` is retained only as the historical configuration input needed to reproduce the legacy-adapter diagnostic; it is not the current QC-RFSC selection manifest.

Reproduction on Euler:

```sh
sbatch bootstrap_empirical_qcrfsc.sbatch
```

After retrieving only aggregate output files into `generated/review_validation/qcrfsc_bootstrap_20261122/<dataset>/`:

```sh
python3 scripts/benchmark_revision/integrate_qcrfsc_bootstrap.py
python3 scripts/condense_supplementary_tables.py
Rscript scripts/generate_revision_supplementary_uncertainty_panels.R
Rscript scripts/generate_revision_main_panels.R
```

The integration validates 2,000 finite replicates per endpoint, independently recomputes SEs and percentile bounds, checks unchanged points, restores QC-RFSC in Figures 3–4 and supplementary uncertainty panels, and appends current-matrix matching intervals and paired differences. Raw data, corrected matrices, bootstrap identity draws and reversible clinical ID maps remain on Euler. Public input/fit/evaluator/draw hashes, aggregate replicates and checks are retained alongside the results.
