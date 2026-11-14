# Current supplementary table sources

These CSVs are the current publication sources for the permutation-calibrated WiNN development version 0.1.7.9999 (package commit `ab4969792c98bd86208a28a4e83b7a69e204b9fc`, 14 November 2026). Only WiNN outputs were recomputed; raw and comparator rows are preserved from the preceding publication snapshot, including the corrected QC-RFSC results. The table numbers below (S4–S31) are the detailed generator outputs. The published supplement condenses them into 13 tables (see "Published numbering" at the end); tables that are not printed remain available here as source data.

The detailed tables S4–S20 begin with aggregate outputs from the benchmark and its derivative evaluation; the validated 27 October RFSC integration then replaces simulation RFSC results and removes obsolete empirical RFSC estimates. The previous public analysis base was `21de7c448bfa5545d27c5a930161fb0bd8e40d85`; it did not contain the full derivative export or the new revision analyses. This revision publishes the exact aggregate inputs needed for table assembly, compressed with an uncompressed SHA-256 manifest under `results/final/revision_inputs/supplement_20261114/` in the analysis repository (the earlier 0.1.4 and 0.1.6 inputs remain in their dated archives). Clinical inputs there are aggregate dataset/method/metric results only.

From the public analysis repository root, regenerate the current tables without refitting any correction model:

```sh
Rscript manuscript/scripts/generate_supplementary_tables.R --final-root=results/final/revision_inputs/supplement_20261114
Rscript manuscript/scripts/nearest_baseline/generate_validation_tables.R
python3 manuscript/scripts/generate_review2_support_tables.py
Rscript manuscript/scripts/generate_s24_paired_effects_table.R
python3 manuscript/scripts/condense_supplementary_tables.py
```

The first command automatically calls `manuscript/scripts/benchmark_revision/integrate_completed_results.py`. That integration checks the completed Euler outputs, replaces all RFSC simulation points, removes earlier-adapter empirical values/intervals, synchronizes the selected-setting manifest, and writes S26. The next commands write the corrected retrieval and baseline comparisons and the stage/design tables. The `--analysis-root` option remains supported for a historical local analysis tree containing `results/final`; it is not a claim that the old public base contains the added derivative files. `supplementary_table_manifest.csv` covers only the S4–S20 tables produced by the original assembler and records the current revision step.

From an Overleaf manuscript checkout, use the equivalent command with an absolute path to the aggregate archive in the public analysis checkout:

```sh
Rscript scripts/generate_supplementary_tables.R --final-root=/path/to/WiNN-for-technical-correction/results/final/revision_inputs/supplement_20261114
Rscript scripts/nearest_baseline/generate_validation_tables.R
python3 scripts/generate_review2_support_tables.py
Rscript scripts/generate_s24_paired_effects_table.R
python3 scripts/condense_supplementary_tables.py
```

Revision table mappings:

- **S21:** `endpoint_free_selection_manifest.csv`; all 54 dataset/output settings, with corrected RFSC selected parameters and fit seeds. The printed table focuses on four comparator settings.
- **S22:** `intensity_domain_safeguard.csv` and `tiger_nonfinite_feature_fallback.csv`, restricted to primary fits with a recorded trigger. Corrected RFSC raw-restoration counts are in S26.
- **S23:** `../review_validation/frozen_identity_retrieval_estimates.csv` and `frozen_identity_retrieval_paired_differences.csv`. All current intervals use 2,000 paired biological-identity resamples with query-count weighting aligned to the point estimator. Earlier identity-weighted intervals remain only in explicitly suffixed audit files.
- **S24:** `../review_validation/paired_feature_effects_descriptive.csv`; matched-feature effects are relative to raw input, not known biological truth.
- **S25:** `../review_validation/nearest_baseline_metrics.csv`; reference precision and batch association for the regression baselines.
- **S26:** `table_s26_qcrfsc_empirical.csv` and `../review_validation/qcrfsc_completed_empirical_all_endpoints.csv`; all five completed corrected RFSC empirical point evaluations. There are no new RFSC empirical uncertainty intervals. The comparable MTBLS79 biological-group point uses 20 mean log1p sample profiles, reevaluated from the archived corrected matrix in `../review_validation/qcrfsc_mtbls79_sample_means.csv`; its raw endpoints reproduce the main figure, and input hashes are supplied alongside it. Figure 4 and the condensed benchmark table use this sample-level point and the valid 568-sample WaveICA group point.
- **S27:** `../review_validation/baseline_simulation_summary.csv` and `baseline_simulation_paired_differences.csv`; the exact ten simulation seeds, with paired t intervals.
- **S28:** `../review_validation/baseline_empirical_feature_effects.csv` and `baseline_empirical_lobo_effects.csv`; descriptive feature effects against input and consistency across batches.
- **S29:** `../review_validation/preparation_biology_summary.csv`; the separate MAD, PQN and combined stage experiments.
- **S30:** `../review_validation/prepared_design_group_by_batch.csv` and `prepared_design_run_order_balance.csv`; group counts and acquisition-order association. Quartile detail is in `design_group_acquisition_quartiles.csv` and `design_covariate_acquisition_quartiles.csv`.
- **S31:** `../review_validation/prepared_design_reference_placement.csv`; aggregate counts, positions and spacing of training/evaluation references.

Current figure generation uses included numerical inputs. From the manuscript directory:

```sh
Rscript scripts/generate_revision_main_panels.R
Rscript scripts/generate_revision_ablation_panels.R
Rscript scripts/generate_revision_supplementary_uncertainty_panels.R
Rscript scripts/generate_revision_s6_panels.R
Rscript scripts/generate_review_fold_panels.R
Rscript scripts/calibration_revision/plot_gate_validation.R
python3 scripts/assemble_main_figures.py
python3 scripts/assemble_supplementary_figures.py
```

Main/S5 sources are in `../revision_panel_inputs/`. The S6 exact 5,000-draw null distributions are supplied compressed as `../revision_panel_inputs/permutation_null.csv.gz`; S7 uses the included `leave_one_batch_out_bootstrap_summary.csv` and displays individual batch folds. Figure S11 uses the independently validated summaries in `../calibration_revision/`. These rendering commands do not perform correction or tune parameters. The drift gate uses permutation maxT familywise adjustment; batch ANOVA and biological-association tests retain their separate BH adjustments. The residual Ljung–Box endpoint remains the frozen asymptotic/BH diagnostic for comparability, not the new production gate.

## Published numbering (condensed supplement)

`scripts/condense_supplementary_tables.py` merges related detailed tables into `generated/supplementary_tables/condensed_*.tex`. The printed supplement uses:

| Published | Content | Detailed source |
|---|---|---|
| S1 | Benchmark designs and independent units | inline in `supplementary.tex` |
| S2 | Compared outputs and information use | inline in `supplementary.tex` |
| S3 | Selected settings for every dataset and method | S21 + S17 |
| S4 | Complete benchmark results | S4–S9 + S26 (corrected QC-RFSC empirical points) |
| S5 | Runtime, retention and WiNN correction counts | S10 |
| S6 | Selective versus forced correction | S15 (compact) + `figures/source_data/S2D–S2F` |
| S7 | Partial confounding | S13 |
| S8 | Biological associations and permutation nulls | S11 |
| S9 | Matched-feature effects | S24 |
| S10 | Cross-batch sample matching | S23 |
| S11 | Same-information baselines | S27 + S25 + S28 |
| S12 | Preparation and PQN response patterns | S29 |
| S13 | Acquisition design (group balance) | S30 |

Not printed, available as source data: the uncorrected technical summary (now the Raw rows of S4), S12 leave-one-batch-out summaries (shown in Supplementary Figure S7), S14 reference splits (Supplementary Figure S4), S16 cumulative stages (Supplementary Figures S1 and S9), S18/S22 numerical safeguards and S20 non-estimable analyses (summarized in Section S9), S19 release validation, and S31 reference placement.


## Isolated table reproduction check

The included driver deletes each producer's declared outputs in an isolated
copy before running the five table generators, then compares the regenerated
TeX bytes and parsed CSV rows with this publication snapshot. Preserved static
tables and staged inputs are labelled separately in the proof. This checks
table assembly from the supplied aggregates; it does not refit correction models.

From the manuscript directory, choose a new replay directory outside this checkout:

```sh
python3 scripts/calibration_revision/reproduce_publication_tables.py \
  --source=. \
  --archive=/path/to/WiNN-for-technical-correction/results/final/revision_inputs/supplement_20261114 \
  --work=/tmp/winn_table_replay \
  --proof=generated/calibration_revision/table_reproduction_checks.csv
```
