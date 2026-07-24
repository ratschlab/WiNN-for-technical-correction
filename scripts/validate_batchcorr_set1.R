#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript or R --file.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))

result_dir <- file.path(repo_root, "results", "batchcorr_set1")
processed_dir <- file.path(repo_root, "data", "public", "batchcorr_set1", "processed")
oral_dir <- file.path(repo_root, "data", "public", "oral_cavity")

required_methods <- c(
  "Raw", "ComBat", "QC-RLSC", "QC-RFSC", "TIGER", "SERRF",
  "WINN auto (QC)", "WINN auto-batch (QC)", "WINN default (no QC)"
)
method_files <- c(
  Raw = "raw.rds", ComBat = "combat.rds", `QC-RLSC` = "qc_rlsc.rds",
  `QC-RFSC` = "qc_rfsc.rds", TIGER = "tiger.rds", SERRF = "serrf.rds",
  `WINN auto (QC)` = "winn_auto_qc_.rds",
  `WINN auto-batch (QC)` = "winn_auto_batch_qc_.rds",
  `WINN default (no QC)` = "winn_default_no_qc_.rds"
)

required_files <- c(
  "method_metrics.csv", "method_metrics_long.csv", "method_parameters.csv",
  "runtime.csv", "qc_cv_by_feature.csv", "run_order_gam_by_feature_batch.csv",
  "ljung_box_by_feature_batch.csv", "accession_associations.csv",
  "replicate_agreement_by_accession.csv", "feature_replicate_repeatability.csv",
  "analysis_log.txt", "sessionInfo.txt", "analysis_complete.json"
)

for (path in c(
  file.path(processed_dir, "BatchCorr_set1_metadata.csv"),
  file.path(result_dir, "method_attempts.csv"),
  file.path(result_dir, "hidden_qc_method_protocol.csv"),
  file.path(result_dir, "method_parameters.csv"),
  file.path(result_dir, "method_warnings_messages_errors.csv"),
  file.path(result_dir, "qc_holdout_ids.csv"),
  file.path(result_dir, "qc_training_ids.csv")
)) {
  if (!file.exists(path)) stop("Required validation input is missing: ", path)
}

meta <- read.csv(file.path(processed_dir, "BatchCorr_set1_metadata.csv"), stringsAsFactors = FALSE)
attempts <- read.csv(file.path(result_dir, "method_attempts.csv"), stringsAsFactors = FALSE)
protocol <- read.csv(file.path(result_dir, "hidden_qc_method_protocol.csv"), stringsAsFactors = FALSE)
holdout <- read.csv(file.path(result_dir, "qc_holdout_ids.csv"), stringsAsFactors = FALSE)
training <- read.csv(file.path(result_dir, "qc_training_ids.csv"), stringsAsFactors = FALSE)
parameters <- read.csv(file.path(result_dir, "method_parameters.csv"), stringsAsFactors = FALSE)
method_diagnostics <- read.csv(file.path(result_dir, "method_warnings_messages_errors.csv"), stringsAsFactors = FALSE)
withdrawn_holdout_ids <- c(
  "LCMS sample 21", "LCMS sample 100", "LCMS sample 198", "LCMS sample 249",
  "LCMS sample 324", "LCMS sample 404", "LCMS sample 502", "LCMS sample 582",
  "LCMS sample 642", "LCMS sample 741"
)

read_matrix_checks <- function(label, filename) {
  path <- file.path(result_dir, "method_outputs", filename)
  if (!file.exists(path)) {
    return(data.frame(method = label, exists = FALSE, features = NA_integer_, samples = NA_integer_, aligned = FALSE, finite = FALSE))
  }
  mat <- readRDS(path)
  expected_samples <- meta$sample_id[meta$sample_id %in% colnames(mat)]
  data.frame(
    method = label,
    exists = TRUE,
    features = nrow(mat),
    samples = ncol(mat),
    aligned = identical(colnames(mat), expected_samples),
    finite = is.matrix(mat) && all(is.finite(mat)) && !anyDuplicated(rownames(mat)) && !anyDuplicated(colnames(mat)),
    stringsAsFactors = FALSE
  )
}
matrix_checks <- do.call(rbind, Map(read_matrix_checks, names(method_files), unname(method_files)))

study <- meta[as.logical(meta$is_study_sample), , drop = FALSE]
study <- study[order(study$accession_id, study$run_order), , drop = FALSE]
expected_rep <- ave(study$run_order, study$accession_id, FUN = seq_along)
raw_archive_names <- c("pos1.zip", "pos2.zip", "neg1.zip", "neg2.zip")
raw_archive_paths <- unlist(lapply(raw_archive_names, function(name) {
  list.files(oral_dir, pattern = paste0("^", name, "$"), recursive = TRUE, full.names = TRUE)
}), use.names = FALSE)

checks <- data.frame(
  check = c(
    "suitability_and_qcrlsc_correction_reports_exist",
    "oral_raw_spectra_absent",
    "all_nine_methods_attempted_in_order",
    "all_nine_methods_completed",
    "holdout_has_one_interior_qc_per_batch",
    "fresh_holdout_seed_and_zero_overlap",
    "training_has_at_least_four_qcs_per_batch",
    "heldout_labels_never_supplied",
    "heldout_values_never_used_for_tuning",
    "qcrlsc_uses_corrected_fixed_protocol",
    "qcrlsc_has_no_stale_tuning_or_warnings",
    "qcrlsc_retains_full_matrix",
    "all_method_matrices_saved",
    "all_method_matrices_align_with_metadata",
    "all_method_matrices_finite",
    "sample_and_qc_counts_match_audit",
    "replicate_ids_are_genuine_accession_repeats",
    "all_primary_result_tables_exist",
    "figure_source_data_exist",
    "clean_rendered_notebook_exists"
  ),
  passed = c(
    all(file.exists(file.path(repo_root, "reports", c(
      "oral_cavity_suitability.md", "batchcorr_set1_suitability.md",
      "batchcorr_set1_qc_rlsc_correction.md", "batchcorr_set1_qc_rlsc_correction.csv"
    )))),
    length(raw_archive_paths) == 0L,
    identical(attempts$method, required_methods) && nrow(attempts) == 9L,
    all(attempts$status == "completed") && nrow(attempts) == 9L,
    nrow(holdout) == 10L && length(unique(holdout$batch)) == 10L && all(holdout$within_batch_order > 1L),
    all(holdout$seed == 20260809L) && !any(holdout$sample_id %in% withdrawn_holdout_ids),
    all(table(training$batch) >= 4L) && length(unique(training$batch)) == 10L,
    nrow(protocol) == 9L && all(!as.logical(protocol$supplied_holdout_qc_ids)),
    nrow(protocol) == 9L && all(!as.logical(protocol$tuning_uses_holdout_values)),
    any(
      parameters$method == "QC-RLSC" &
        grepl("analysis_scale=log1p", parameters$parameters, fixed = TRUE) &
        grepl("span=1", parameters$parameters, fixed = TRUE) &
        grepl("degree=1", parameters$parameters, fixed = TRUE) &
        grepl("sample_ids_preserved=TRUE", parameters$parameters, fixed = TRUE)
    ),
    !file.exists(file.path(result_dir, "tuning_candidates", "qc_rlsc.csv")) &&
      all(is.na(method_diagnostics$warnings[method_diagnostics$method == "QC-RLSC"]) |
            method_diagnostics$warnings[method_diagnostics$method == "QC-RLSC"] == "") &&
      all(is.na(method_diagnostics$error[method_diagnostics$method == "QC-RLSC"]) |
            method_diagnostics$error[method_diagnostics$method == "QC-RLSC"] == ""),
    matrix_checks$features[matrix_checks$method == "QC-RLSC"] == 199L &&
      matrix_checks$samples[matrix_checks$method == "QC-RLSC"] == 761L,
    all(matrix_checks$exists),
    all(matrix_checks$aligned) && all(matrix_checks$samples == 761L),
    all(matrix_checks$finite),
    nrow(meta) == 761L && sum(as.logical(meta$is_qc)) == 51L && nrow(study) == 710L && length(unique(meta$batch)) == 10L,
    identical(as.integer(study$biological_replicate), as.integer(expected_rep)) && identical(study$replicate_group, study$accession_id),
    all(file.exists(file.path(result_dir, required_files))),
    length(list.files(file.path(result_dir, "figures", "source_data"), pattern = "\\.csv$")) >= 7L,
    file.exists(file.path(repo_root, "notebooks", "rendered", "batchcorr_set1_comparison.html")) &&
      file.info(file.path(repo_root, "notebooks", "rendered", "batchcorr_set1_comparison.html"))$size > 1e6
  ),
  stringsAsFactors = FALSE
)
checks$details <- c(
  "Oral, BatchCorr, and QC-RLSC correction reports",
  paste("found", length(raw_archive_paths), "of four raw ZIP names"),
  paste(attempts$method, collapse = "; "),
  paste(sum(attempts$status == "completed"), "completed"),
  paste(nrow(holdout), "holdouts across", length(unique(holdout$batch)), "batches"),
  paste("seed", unique(holdout$seed), "; overlap", sum(holdout$sample_id %in% withdrawn_holdout_ids)),
  paste(range(as.integer(table(training$batch))), collapse = "-"),
  paste(sum(as.logical(protocol$supplied_holdout_qc_ids)), "exposed method rows"),
  paste(sum(as.logical(protocol$tuning_uses_holdout_values)), "exposed tuning rows"),
  parameters$parameters[parameters$method == "QC-RLSC"],
  paste("stale tuning file", file.exists(file.path(result_dir, "tuning_candidates", "qc_rlsc.csv")),
        "; warning/error characters",
        sum(nchar(method_diagnostics$warnings[method_diagnostics$method == "QC-RLSC"]), na.rm = TRUE) +
          sum(nchar(method_diagnostics$error[method_diagnostics$method == "QC-RLSC"]), na.rm = TRUE)),
  paste(matrix_checks$features[matrix_checks$method == "QC-RLSC"], "features x",
        matrix_checks$samples[matrix_checks$method == "QC-RLSC"], "injections"),
  paste(sum(matrix_checks$exists), "of", nrow(matrix_checks), "matrices"),
  paste(range(matrix_checks$samples), collapse = "-"),
  paste(sum(matrix_checks$finite), "of", nrow(matrix_checks), "finite matrices"),
  paste(nrow(meta), "total;", sum(as.logical(meta$is_qc)), "QC;", nrow(study), "study"),
  paste(length(unique(study$accession_id)), "accession IDs; no pseudo-replicates"),
  paste(sum(file.exists(file.path(result_dir, required_files))), "of", length(required_files), "files"),
  paste(length(list.files(file.path(result_dir, "figures", "source_data"), pattern = "\\.csv$")), "CSV source tables"),
  paste(file.info(file.path(repo_root, "notebooks", "rendered", "batchcorr_set1_comparison.html"))$size, "bytes")
)

write.csv(checks, file.path(result_dir, "final_validation.csv"), row.names = FALSE, quote = TRUE)
write.csv(matrix_checks, file.path(result_dir, "method_matrix_validation.csv"), row.names = FALSE, quote = TRUE)
jsonlite::write_json(
  list(
    validated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    all_passed = all(checks$passed),
    checks = checks,
    matrices = matrix_checks
  ),
  file.path(result_dir, "final_validation.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  dataframe = "rows"
)

if (!all(checks$passed)) {
  stop("Final validation failed: ", paste(checks$check[!checks$passed], collapse = ", "))
}
message("All ", nrow(checks), " final validation checks passed.")
