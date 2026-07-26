#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
final_root <- file.path(repo_root, "results", "final")

required <- c(
  "aggregate_manifest.json", "file_manifest.csv",
  "tables/primary_method_metrics_long.csv",
  "tables/combined_primary_benchmark.csv",
  "tables/winn_mode_comparison.csv",
  "tables/feature_retention_and_correction_counts.csv",
  "tables/configuration_selection_manifest.csv",
  "tables/primary_competitor_candidate_scores.csv",
  "tables/primary_winn_selected_parameters.csv",
  "tables/competitor_reuse_audit.csv",
  "tables/primary_run_provenance.csv",
  "tables/runtime_measurement_policy.csv",
  "tables/intensity_domain_safeguard.csv",
  "tables/by_dataset/simulation_method_metrics.csv",
  "tables/by_dataset/mtbls79_method_metrics.csv",
  "tables/by_dataset/clinical_fiams_method_metrics.csv",
  "tables/by_dataset/batchcorr_set1_method_metrics.csv",
  "tables/by_dataset/sacurine_method_metrics.csv",
  "tables/by_dataset/waveica_method_metrics.csv",
  "tables/simulation_seed_summary.csv",
  "tables/reference_split_summary.csv",
  "tables/reference_split_selection.csv",
  "tables/reference_split_selection_stability.csv",
  "tables/reference_split_candidate_scores.csv",
  "tables/partial_confounding_summary.csv",
  "tables/ablation_step_impacts_heldout_qc_cv.csv",
  "tables/ablation_step_impacts_gam_deviance.csv",
  "tables/ablation_step_impacts_batch_r2.csv",
  "tables/ablation_step_impacts_metabolite_icc_a1.csv",
  "tables/ablation_step_impacts_sample_pearson.csv",
  "tables/ablation_step_impacts_feature_repeatability_ratio.csv",
  "figures/figure_manifest.csv",
  "validation/release_validation.csv",
  "validation/task_completeness.csv",
  "validation/execution_ledger.csv",
  "validation/logical_task_ledger.csv",
  "validation/reference_separation_audit.csv",
  "validation/corrected_matrix_domain_audit.csv"
)
missing <- required[!file.exists(file.path(final_root, required))]
if (length(missing)) stop("Missing final result file(s): ", paste(missing, collapse = ", "), call. = FALSE)

primary <- read.csv(file.path(final_root, "tables", "primary_method_metrics_long.csv"), stringsAsFactors = FALSE)
if (length(unique(primary$dataset_key)) != 6L) stop("Primary table does not contain six datasets.", call. = FALSE)
if (length(unique(primary$method_id)) != 9L) stop("Primary table does not contain nine methods.", call. = FALSE)

completion <- read.csv(file.path(final_root, "validation", "task_completeness.csv"), stringsAsFactors = FALSE)
if (!all(completion$complete)) stop("One or more production families are incomplete.", call. = FALSE)
validation <- read.csv(file.path(final_root, "validation", "release_validation.csv"), stringsAsFactors = FALSE)
if (!all(validation$passed)) stop("The saved release validation contains failures.", call. = FALSE)

flooring <- read.csv(
  file.path(final_root, "tables", "intensity_domain_safeguard.csv"),
  stringsAsFactors = FALSE
)
domain <- read.csv(
  file.path(final_root, "validation", "corrected_matrix_domain_audit.csv"),
  stringsAsFactors = FALSE
)
if (nrow(flooring) != 684L || nrow(domain) != 684L ||
    any(domain$negative_values != 0L) || any(domain$nonfinite_values != 0L)) {
  stop("The committed output-domain audit is incomplete or contains invalid matrices.", call. = FALSE)
}

figures <- read.csv(file.path(final_root, "figures", "figure_manifest.csv"), stringsAsFactors = FALSE)
source_ok <- file.exists(file.path(final_root, figures$source_csv))
pdf_ok <- !nzchar(figures$pdf) | file.exists(file.path(final_root, figures$pdf))
png_ok <- !nzchar(figures$png) | file.exists(file.path(final_root, figures$png))
if (!all(source_ok & pdf_ok & png_ok)) stop("A final figure or source table is missing.", call. = FALSE)

manifest <- read.csv(file.path(final_root, "file_manifest.csv"), stringsAsFactors = FALSE)
manifest <- manifest[manifest$path != "file_manifest.csv", , drop = FALSE]
manifest_paths <- file.path(final_root, manifest$path)
if (!all(file.exists(manifest_paths))) stop("The final output manifest lists missing files.", call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("The digest package is required.", call. = FALSE)
observed_sha <- vapply(manifest_paths, digest::digest, character(1), file = TRUE, algo = "sha256")
if (!identical(unname(observed_sha), unname(manifest$sha256))) {
  stop("A committed final output does not match its SHA-256 manifest.", call. = FALSE)
}

cat(
  "Validated final release:", length(unique(primary$dataset_key)), "datasets,",
  length(unique(primary$method_id)), "methods,", nrow(figures),
  "figure-source entries,", nrow(manifest), "verified files.\n"
)
