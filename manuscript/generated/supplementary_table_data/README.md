# Current supplementary table sources

These CSVs are the current publication sources. Tables S4–S20 begin with aggregate outputs from the original benchmark and its local derivative evaluation; the validated 27 October RFSC integration then replaces simulation RFSC results and removes obsolete empirical RFSC estimates. The previous public analysis base was `21de7c448bfa5545d27c5a930161fb0bd8e40d85`; it did not contain the full derivative export or the new revision analyses. This revision publishes the exact aggregate inputs needed for table assembly, compressed with an uncompressed SHA-256 manifest under `results/final/revision_inputs/supplement_20260727/` in the analysis repository. Clinical inputs there are aggregate dataset/method/metric results only.

From the public analysis repository root, regenerate the current tables without refitting any correction model:

```sh
Rscript manuscript/scripts/generate_supplementary_tables.R --final-root=results/final/revision_inputs/supplement_20260727
Rscript manuscript/scripts/nearest_baseline/generate_validation_tables.R
python3 manuscript/scripts/generate_review2_support_tables.py
```

The first command automatically calls `manuscript/scripts/benchmark_revision/integrate_completed_results.py`. That integration checks the completed Euler outputs, replaces all RFSC simulation points, removes earlier-adapter empirical values/intervals, synchronizes the selected-setting manifest, and writes S26. The next commands write the corrected retrieval and baseline comparisons and the stage/design tables. The `--analysis-root` option remains supported for a historical local analysis tree containing `results/final`; it is not a claim that the old public base contains the added derivative files. `supplementary_table_manifest.csv` covers only the S4–S20 tables produced by the original assembler and records the current revision step.

From an Overleaf manuscript checkout, use the equivalent command with an absolute path to the aggregate archive in the public analysis checkout:

```sh
Rscript scripts/generate_supplementary_tables.R --final-root=/path/to/WiNN-for-technical-correction/results/final/revision_inputs/supplement_20260727
Rscript scripts/nearest_baseline/generate_validation_tables.R
python3 scripts/generate_review2_support_tables.py
```

Revision table mappings:

- **S21:** `endpoint_free_selection_manifest.csv`; all 54 dataset/output settings, with corrected RFSC selected parameters and fit seeds. The printed table focuses on four comparator settings.
- **S22:** `intensity_domain_safeguard.csv` and `tiger_nonfinite_feature_fallback.csv`, restricted to primary fits with a recorded trigger. Corrected RFSC raw-restoration counts are in S26.
- **S23:** `../review_validation/frozen_identity_retrieval_estimates.csv` and `frozen_identity_retrieval_paired_differences.csv`. All current intervals use 2,000 paired biological-identity resamples with query-count weighting aligned to the point estimator. Earlier identity-weighted intervals remain only in explicitly suffixed audit files.
- **S24:** `../review_validation/paired_feature_effects_descriptive.csv`; matched-feature effects are relative to raw input, not known biological truth.
- **S25:** `../review_validation/nearest_baseline_metrics.csv`; reference precision and batch association for the regression baselines.
- **S26:** `table_s26_qcrfsc_empirical.csv` and `../review_validation/qcrfsc_completed_empirical_all_endpoints.csv`; all five completed corrected RFSC empirical point evaluations. There are no new RFSC empirical uncertainty intervals.
- **S27:** `../review_validation/baseline_simulation_summary.csv` and `baseline_simulation_paired_differences.csv`; the exact ten simulation seeds, with paired t intervals.
- **S28:** `../review_validation/baseline_empirical_feature_effects.csv` and `baseline_empirical_lobo_effects.csv`; descriptive feature effects against input and consistency across batches.
- **S29:** `../review_validation/preparation_biology_summary.csv`; the separate MAD, PQN and combined stage experiments.
- **S30:** `../review_validation/prepared_design_group_by_batch.csv` and `prepared_design_run_order_balance.csv`; group counts and acquisition-order association. Quartile detail is in `design_group_acquisition_quartiles.csv` and `design_covariate_acquisition_quartiles.csv`.
- **S31:** `../review_validation/prepared_design_reference_placement.csv`; aggregate counts, positions and spacing of training/evaluation references.

Current figure generation uses included numerical inputs. From the manuscript directory:

```sh
Rscript scripts/generate_revision_main_panels.R
Rscript scripts/generate_revision_supplementary_uncertainty_panels.R
Rscript scripts/generate_revision_s6_panels.R
Rscript scripts/generate_review_fold_panels.R
python3 scripts/assemble_main_figures.py
python3 scripts/assemble_supplementary_figures.py
```

Main/S5 sources are in `../revision_panel_inputs/`. The S6 exact 5,000-draw null distributions are supplied compressed as `../revision_panel_inputs/permutation_null.csv.gz`; S7 uses the included `leave_one_batch_out_bootstrap_summary.csv` and displays individual batch folds. These rendering commands do not perform correction or tune parameters.
