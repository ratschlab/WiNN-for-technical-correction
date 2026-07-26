#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
task_value <- grep("^--task-id=", args, value = TRUE)
if (length(task_value) != 1L) stop("Supply --task-id=CONF_NNN.", call. = FALSE)
task_id <- sub("^--task-id=", "", task_value)

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root)) stop("WINN_RELEASE_ROOT is not set.", call. = FALSE)
source(file.path(release_root, "analysis", "release_helpers.R"), local = FALSE)
source(file.path(source_root, "scripts", "robustness", "simulation_core.R"), local = FALSE)
source(file.path(source_root, "scripts", "robustness", "partial_confounding_core.R"), local = FALSE)
source(file.path(source_root, "scripts", "robustness", "canonical_cache.R"), local = FALSE)
source(file.path(source_root, "scripts", "robustness", "partial_confounding_portable_hash.R"), local = FALSE)
source(file.path(source_root, "scripts", "robustness", "partial_confounding_pilot_metrics.R"), local = FALSE)
source(file.path(source_root, "scripts", "robustness", "partial_confounding_full_metrics.R"), local = FALSE)
source(
  file.path(source_root, "scripts", "run_order_drift_helpers.R"),
  local = FALSE
)

manifest <- read.csv(
  file.path(release_root, "analysis", "config", "partial_confounding_task_manifest.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
task <- manifest[manifest$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L) stop("Unknown or duplicate task: ", task_id, call. = FALSE)
scenario_id <- as.character(task$scenario_id)
seed_id <- as.character(task$seed_id)
output_dir <- file.path(release_root, "results", "partial_confounding", scenario_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
completion_path <- file.path(output_dir, "complete.json")
if (file.exists(completion_path) && !"--force" %in% args) {
  message(task_id, " already completed.")
  quit(save = "no", status = 0L)
}

scenario_root <- file.path(source_root, "data", "simulated", "partial_confounding")
bundle_path <- file.path(scenario_root, "bundles", paste0(scenario_id, ".rds"))
if (!file.exists(bundle_path)) stop("Missing frozen scenario bundle: ", bundle_path)
bundle <- readRDS(bundle_path)
if (!identical(as.character(bundle$scenario$scenario_id), scenario_id)) {
  stop("Scenario bundle identity mismatch.", call. = FALSE)
}

started <- Sys.time()
base <- generate_canonical_simulation(
  bundle$generation_seeds, include_artifact_matrices = TRUE
)
allocation <- bundle$evaluation_truth$phenotype_allocation
feature_truth <- bundle$evaluation_truth$feature_truth
constructed <- construct_partial_confounding_matrices(
  base, allocation, feature_truth, include_intensities = TRUE
)
scenario_validation <- validate_partial_confounding_scenario(
  base, allocation, feature_truth, constructed
)
if (is.data.frame(scenario_validation) && "passed" %in% names(scenario_validation) &&
    !all(scenario_validation$passed)) {
  stop("Frozen scenario reconstruction failed validation.", call. = FALSE)
}
x <- constructed$raw_intensity
clean <- constructed$clean_ground_truth
portable_hashes <- read.csv(
  file.path(scenario_root, "config", "portable_matrix_hashes.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
portable_row <- portable_hashes[portable_hashes$scenario_id == scenario_id, , drop = FALSE]
if (nrow(portable_row) != 1L ||
    !identical(
      partial_confounding_portable_matrix_sha256(x),
      portable_row$raw_intensity_quantized_9_sha256[[1L]]
    ) ||
    !identical(
      partial_confounding_portable_matrix_sha256(clean),
      portable_row$clean_ground_truth_quantized_9_sha256[[1L]]
    )) {
  stop("Reconstructed scenario matrices do not match the frozen portable hashes.", call. = FALSE)
}

metadata <- constructed$correction_metadata
metadata$batch <- metadata$plate
metadata$within_batch_order <- metadata$order_in_plate
hidden_seeds <- read.csv(
  if (dir.exists(file.path(scenario_root, "full_grid"))) {
    file.path(scenario_root, "full_grid", "config", "hidden_reference_seeds.csv")
  } else {
    file.path(scenario_root, "config", "hidden_reference_seeds.csv")
  },
  stringsAsFactors = FALSE, check.names = FALSE
)
hidden_row <- hidden_seeds[hidden_seeds$seed_id == seed_id, , drop = FALSE]
if (nrow(hidden_row) != 1L) stop("Missing hidden-reference seed for ", seed_id)
hidden_assignment <- select_canonical_hidden_references(
  base$metadata,
  seed = as.integer(hidden_row$hidden_reference_seed),
  hidden_per_plate = as.integer(hidden_row$hidden_per_plate)
)
hidden_ids <- as.character(hidden_assignment$sample_id[hidden_assignment$is_hidden_reference])
training_ids <- as.character(hidden_assignment$sample_id[!hidden_assignment$is_hidden_reference])
meta_blind <- metadata[, c("sample_id", "batch", "run_order", "within_batch_order"), drop = FALSE]
meta_blind$class <- ifelse(meta_blind$sample_id %in% training_ids, "QC", "Sample")
meta_blind$is_qc <- meta_blind$sample_id %in% training_ids
meta_blind$sample_type <- ifelse(meta_blind$is_qc, "control", "sample")

set.seed(as.integer(bundle$component_seeds[["method_WINN_fixed_no_QC"]]))
capture <- release_capture(function() winn::winn(
  x,
  batch = meta_blind$batch,
  run_order = meta_blind$run_order,
  control_samples = NULL,
  parameters = "fixed",
  ljung_box_fitdf = 0L,
  return_details = TRUE
))
writeLines(
  c("WARNINGS", capture$warnings, "", "MESSAGES", capture$messages, "", "ERROR", capture$error),
  file.path(output_dir, "conditions.log")
)
if (inherits(capture$value, "error")) stop(capture$error, call. = FALSE)
winn_result <- capture$value
winn_matrix <- release_validate_output(winn_result$data, list(x = x), "winn_fixed_no_qc")
method_matrices <- list("Raw" = x, "WINN default (no QC)" = winn_matrix)
runtime <- data.frame(
  method = names(method_matrices),
  runtime_sec = c(0, capture$runtime_sec),
  stringsAsFactors = FALSE
)

metric_dictionary <- partial_full_metric_dictionary()
metrics <- compute_partial_confounding_full_metrics(
  scenario_id = scenario_id,
  seed_id = seed_id,
  method_matrices = method_matrices,
  raw = x,
  clean = clean,
  allocation = allocation,
  metadata_evaluation = metadata,
  hidden_ids = hidden_ids,
  feature_truth = feature_truth,
  runtime = runtime,
  ablation_diagnostics = list(),
  metric_dictionary = metric_dictionary,
  run_dir = output_dir
)

release_atomic_save_rds(x, file.path(output_dir, "raw_matrix.rds"))
release_atomic_save_rds(winn_matrix, file.path(output_dir, "winn_fixed_matrix.rds"))
release_atomic_save_rds(
  unclass(winn_result)[setdiff(names(winn_result), "data")],
  file.path(output_dir, "winn_details.rds")
)
write.csv(metrics$method_metrics, file.path(output_dir, "method_metrics.csv"), row.names = FALSE, quote = TRUE)
write.csv(runtime, file.path(output_dir, "runtime.csv"), row.names = FALSE, quote = TRUE)
write.csv(scenario_validation, file.path(output_dir, "scenario_validation.csv"), row.names = FALSE, quote = TRUE)
write.csv(data.frame(
  sample_id = c(training_ids, hidden_ids),
  assignment = c(rep("training", length(training_ids)), rep("hidden", length(hidden_ids))),
  stringsAsFactors = FALSE
), file.path(output_dir, "reference_assignment.csv"), row.names = FALSE, quote = TRUE)
utils::capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
completed <- Sys.time()
jsonlite::write_json(list(
  schema = "endpoint_free_partial_confounding_winn_only_v1",
  task_id = task_id, scenario_id = scenario_id, seed_id = seed_id,
  status = "completed", methods = names(method_matrices),
  n_features = nrow(x), n_samples = ncol(x),
  n_training_qc = length(training_ids), n_hidden_qc = length(hidden_ids),
  raw_input_sha256 = release_sha256_object(x),
  clean_truth_sha256 = release_sha256_object(clean),
  winn_output_sha256 = release_sha256_object(winn_matrix),
  phenotype_fields_available_to_method = FALSE,
  package = frozen_winn,
  runtime_sec = capture$runtime_sec,
  wall_sec = as.numeric(difftime(completed, started, units = "secs")),
  completed_at = format(completed, "%Y-%m-%dT%H:%M:%S%z")
), file.path(output_dir, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
jsonlite::write_json(
  list(task_id = task_id, scenario_id = scenario_id, status = "completed"),
  completion_path, auto_unbox = TRUE, pretty = TRUE
)
message(task_id, " completed: ", scenario_id)
