#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
seed_value <- grep("^--seed-id=", args, value = TRUE)
if (length(seed_value) != 1L) stop("Supply --seed-id=SIM01 through SIM10.", call. = FALSE)
seed_id <- sub("^--seed-id=", "", seed_value)
if (!grepl("^SIM(0[1-9]|10)$", seed_id)) stop("Only SIM01 through SIM10 are in the release.")

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root)) stop("WINN_RELEASE_ROOT is not set.", call. = FALSE)
source(file.path(release_root, "analysis", "release_helpers.R"), local = FALSE)

output_dir <- file.path(release_root, "results", "simulation_seed_ablations", seed_id)
dir.create(file.path(output_dir, "variants"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
completion_path <- file.path(output_dir, "complete.json")
if (file.exists(completion_path) && !"--force" %in% args) {
  message(seed_id, " ablation bundle already completed.")
  quit(save = "no", status = 0L)
}

prepared <- release_load_simulation(seed_id)
prepared$meta_method <- tq_method_metadata(prepared$meta, prepared$training_ids)
prepared$meta_technical <- tq_technical_metadata(prepared$meta)
input_hash <- release_sha256_object(prepared$x)
set.seed(2026072700L + as.integer(sub("SIM", "", seed_id)))
started <- Sys.time()

base <- winn::winn_ablation(
  prepared$x, batch = prepared$meta$batch, run_order = prepared$meta$run_order,
  control_samples = NULL, parameters = "fixed", use_outlier_shrinkage = TRUE,
  drift_gate = "selective", batch_gate = "selective", pqn_mode = "shrink",
  fdr_threshold = 0.05, test = "Ljung-Box", lag = NULL, ljung_box_fitdf = 0L,
  spline_method = "conservative", return_intermediates = TRUE,
  return_diagnostics = TRUE
)
forced_drift <- winn::winn_ablation(
  prepared$x, batch = prepared$meta$batch, run_order = prepared$meta$run_order,
  control_samples = NULL, parameters = "fixed", use_outlier_shrinkage = TRUE,
  drift_gate = "all", batch_gate = "selective", pqn_mode = "shrink",
  fdr_threshold = 0.05, test = "Ljung-Box", lag = NULL, ljung_box_fitdf = 0L,
  spline_method = "conservative", return_intermediates = TRUE,
  return_diagnostics = TRUE
)
variants <- list(
  C0_RAW = empty_raw_result(prepared$x),
  C1_OUTLIER = subset_cumulative_result(base, "outlier"),
  C2_SELECTIVE_DRIFT = subset_cumulative_result(base, "drift"),
  C3_SELECTIVE_BATCH = subset_cumulative_result(base, "batch"),
  C4_FULL_FIXED = base,
  G_SS = base,
  G_AS = forced_drift,
  G_SA = derive_forced_batch_variant(base, prepared$meta$batch),
  G_AA = derive_forced_batch_variant(forced_drift, prepared$meta$batch)
)

fixed <- release_run_winn(prepared, auto_batch = FALSE, parameters = "fixed")
fixed_difference <- max(abs(fixed$data - base$data))
if (!is.finite(fixed_difference) || fixed_difference > 1e-8) {
  stop(seed_id, ": fixed/ablation equality check failed.", call. = FALSE)
}

matrix_rows <- list()
selection_rows <- list()
for (variant_id in names(variants)) {
  result <- variants[[variant_id]]
  matrix <- release_validate_output(result$data, prepared, variant_id)
  matrix_path <- file.path(output_dir, "variants", paste0(variant_id, ".rds"))
  release_atomic_save_rds(matrix, matrix_path)
  release_atomic_save_rds(list(
    configuration = result$configuration,
    diagnostics = result$diagnostics,
    stage_runtime_sec = result$stage_runtime_sec,
    total_runtime_sec = result$total_runtime_sec
  ), file.path(output_dir, "diagnostics", paste0(variant_id, ".rds")))
  drift <- result$diagnostics$drift
  batch <- result$diagnostics$batch
  selection_rows[[variant_id]] <- data.frame(
    seed_id = seed_id, variant_id = variant_id,
    drift_corrected_profiles = if (is.data.frame(drift) && "actual_correction" %in% names(drift)) sum(drift$actual_correction) else 0L,
    drift_unique_features_corrected = if (is.data.frame(drift) && all(c("feature_id", "actual_correction") %in% names(drift))) length(unique(drift$feature_id[drift$actual_correction])) else 0L,
    batch_corrected_features = if (is.data.frame(batch) && "actual_correction" %in% names(batch)) sum(batch$actual_correction) else 0L,
    total_runtime_sec = result$total_runtime_sec,
    stringsAsFactors = FALSE
  )
  matrix_rows[[variant_id]] <- data.frame(
    seed_id = seed_id, variant_id = variant_id,
    n_features = nrow(matrix), n_samples = ncol(matrix),
    object_sha256 = release_sha256_object(matrix),
    file_sha256 = release_sha256_file(matrix_path), stringsAsFactors = FALSE
  )
}
write.csv(dplyr::bind_rows(matrix_rows), file.path(output_dir, "matrix_manifest.csv"), row.names = FALSE, quote = TRUE)
write.csv(dplyr::bind_rows(selection_rows), file.path(output_dir, "stage_selection_summary.csv"), row.names = FALSE, quote = TRUE)
utils::capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
completed <- Sys.time()
jsonlite::write_json(list(
  schema = "endpoint_free_simulation_ablation_bundle_v1",
  seed_id = seed_id, status = "completed", logical_task_count = length(variants),
  input_object_sha256 = input_hash,
  fixed_equivalence_max_abs_difference = fixed_difference,
  package = frozen_winn,
  wall_sec = as.numeric(difftime(completed, started, units = "secs")),
  completed_at = format(completed, "%Y-%m-%dT%H:%M:%S%z")
), file.path(output_dir, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
jsonlite::write_json(
  list(seed_id = seed_id, status = "completed", logical_task_count = length(variants)),
  completion_path, auto_unbox = TRUE, pretty = TRUE
)
message(seed_id, " simulation ablation bundle completed.")
