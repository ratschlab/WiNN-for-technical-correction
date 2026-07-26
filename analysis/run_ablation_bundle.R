#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
dataset_value <- grep("^--dataset=", args, value = TRUE)
if (length(dataset_value) != 1L) stop("Supply --dataset=<dataset_key>.", call. = FALSE)
dataset_key <- sub("^--dataset=", "", dataset_value)

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root)) stop("WINN_RELEASE_ROOT is not set.", call. = FALSE)
source(file.path(release_root, "analysis", "release_helpers.R"), local = FALSE)

allowed <- read.csv(
  file.path(release_root, "analysis", "config", "dataset_manifest.csv"),
  stringsAsFactors = FALSE
)$dataset_key
if (!dataset_key %in% allowed) stop("Unknown dataset: ", dataset_key, call. = FALSE)

output_dir <- file.path(release_root, "results", "ablations", dataset_key)
dir.create(file.path(output_dir, "variants"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
completion_path <- file.path(output_dir, "complete.json")
if (file.exists(completion_path) && !"--force" %in% args) {
  message(dataset_key, " ablation bundle already completed.")
  quit(save = "no", status = 0L)
}

set.seed(2026072601L + match(dataset_key, allowed))
started <- Sys.time()
prepared <- release_load_dataset(dataset_key)
input_hash <- release_validate_input(prepared, dataset_key)
x <- prepared$x

base_capture <- release_capture(function() winn::winn_ablation(
  x,
  batch = prepared$meta$batch,
  run_order = prepared$meta$run_order,
  control_samples = NULL,
  parameters = "fixed",
  use_outlier_shrinkage = TRUE,
  drift_gate = "selective",
  batch_gate = "selective",
  pqn_mode = "shrink",
  fdr_threshold = 0.05,
  test = "Ljung-Box",
  lag = NULL,
  ljung_box_fitdf = 0L,
  spline_method = "conservative",
  return_intermediates = TRUE,
  return_diagnostics = TRUE
))
if (inherits(base_capture$value, "error")) stop(base_capture$error, call. = FALSE)

forced_drift_capture <- release_capture(function() winn::winn_ablation(
  x,
  batch = prepared$meta$batch,
  run_order = prepared$meta$run_order,
  control_samples = NULL,
  parameters = "fixed",
  use_outlier_shrinkage = TRUE,
  drift_gate = "all",
  batch_gate = "selective",
  pqn_mode = "shrink",
  fdr_threshold = 0.05,
  test = "Ljung-Box",
  lag = NULL,
  ljung_box_fitdf = 0L,
  spline_method = "conservative",
  return_intermediates = TRUE,
  return_diagnostics = TRUE
))
if (inherits(forced_drift_capture$value, "error")) stop(forced_drift_capture$error, call. = FALSE)

base <- base_capture$value
forced_drift <- forced_drift_capture$value
selective_drift_forced_batch <- derive_forced_batch_variant(base, prepared$meta$batch)
forced_drift_forced_batch <- derive_forced_batch_variant(forced_drift, prepared$meta$batch)
variants <- list(
  C0_RAW = empty_raw_result(x),
  C1_OUTLIER = subset_cumulative_result(base, "outlier"),
  C2_SELECTIVE_DRIFT = subset_cumulative_result(base, "drift"),
  C3_SELECTIVE_BATCH = subset_cumulative_result(base, "batch"),
  C4_FULL_FIXED = base,
  G_SS = base,
  G_AS = forced_drift,
  G_SA = selective_drift_forced_batch,
  G_AA = forced_drift_forced_batch
)

fixed_check <- release_run_winn(prepared, auto_batch = FALSE, parameters = "fixed")
fixed_difference <- max(abs(fixed_check$data - variants$C4_FULL_FIXED$data))
if (!is.finite(fixed_difference) || fixed_difference > 1e-8) {
  stop("winn_ablation full fixed output differs from winn fixed mode.", call. = FALSE)
}

definition <- ablation_variant_table()
selection_rows <- list()
matrix_rows <- list()
for (variant_id in definition$variant_id) {
  result <- variants[[variant_id]]
  matrix <- release_validate_output(result$data, prepared, variant_id)
  matrix_path <- file.path(output_dir, "variants", paste0(variant_id, ".rds"))
  diagnostic_path <- file.path(output_dir, "diagnostics", paste0(variant_id, ".rds"))
  release_atomic_save_rds(matrix, matrix_path)
  release_atomic_save_rds(list(
    configuration = result$configuration,
    diagnostics = result$diagnostics,
    stage_runtime_sec = result$stage_runtime_sec,
    total_runtime_sec = result$total_runtime_sec
  ), diagnostic_path)

  drift <- result$diagnostics$drift
  batch <- result$diagnostics$batch
  pqn <- result$diagnostics$pqn
  drift_corrected <- if (is.data.frame(drift) && nrow(drift) && "actual_correction" %in% names(drift)) {
    sum(as.logical(drift$actual_correction), na.rm = TRUE)
  } else 0L
  drift_features <- if (is.data.frame(drift) && nrow(drift) && all(c("feature_id", "actual_correction") %in% names(drift))) {
    length(unique(drift$feature_id[as.logical(drift$actual_correction)]))
  } else 0L
  batch_corrected <- if (is.data.frame(batch) && nrow(batch) && "actual_correction" %in% names(batch)) {
    sum(as.logical(batch$actual_correction), na.rm = TRUE)
  } else 0L
  pqn_changed <- if (is.data.frame(pqn) && nrow(pqn) && "altered" %in% names(pqn)) {
    sum(as.logical(pqn$altered), na.rm = TRUE)
  } else 0L
  selection_rows[[variant_id]] <- data.frame(
    dataset_key = dataset_key,
    variant_id = variant_id,
    drift_corrected_profiles = drift_corrected,
    drift_unique_features_corrected = drift_features,
    batch_corrected_features = batch_corrected,
    pqn_samples_changed = pqn_changed,
    total_runtime_sec = result$total_runtime_sec,
    stringsAsFactors = FALSE
  )
  matrix_rows[[variant_id]] <- data.frame(
    dataset_key = dataset_key,
    variant_id = variant_id,
    n_features = nrow(matrix),
    n_samples = ncol(matrix),
    object_sha256 = release_sha256_object(matrix),
    file_sha256 = release_sha256_file(matrix_path),
    stringsAsFactors = FALSE
  )
}

write.csv(definition, file.path(output_dir, "variant_definitions.csv"), row.names = FALSE, quote = TRUE)
write.csv(dplyr::bind_rows(selection_rows), file.path(output_dir, "stage_selection_summary.csv"), row.names = FALSE, quote = TRUE)
write.csv(dplyr::bind_rows(matrix_rows), file.path(output_dir, "matrix_manifest.csv"), row.names = FALSE, quote = TRUE)
writeLines(
  c(
    "BASE WARNINGS", base_capture$warnings, "", "BASE MESSAGES", base_capture$messages,
    "", "FORCED-DRIFT WARNINGS", forced_drift_capture$warnings,
    "", "FORCED-DRIFT MESSAGES", forced_drift_capture$messages
  ),
  file.path(output_dir, "conditions.log")
)
utils::capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
completed <- Sys.time()
jsonlite::write_json(list(
  schema = "endpoint_free_winn_ablation_bundle_v1",
  dataset_key = dataset_key,
  status = "completed",
  logical_task_count = length(variants),
  input_object_sha256 = input_hash,
  fixed_equivalence_max_abs_difference = fixed_difference,
  package = frozen_winn,
  started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
  completed_at = format(completed, "%Y-%m-%dT%H:%M:%S%z"),
  wall_sec = as.numeric(difftime(completed, started, units = "secs"))
), file.path(output_dir, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
jsonlite::write_json(
  list(dataset_key = dataset_key, status = "completed", logical_task_count = length(variants)),
  completion_path, auto_unbox = TRUE, pretty = TRUE
)
message(dataset_key, " ablation bundle completed with ", length(variants), " variants.")
