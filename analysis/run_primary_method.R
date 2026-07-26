#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
task_value <- grep("^--task-id=", args, value = TRUE)
if (length(task_value) != 1L) stop("Supply --task-id=PRIMARY_NNN.", call. = FALSE)
task_id <- sub("^--task-id=", "", task_value)

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root)) stop("WINN_RELEASE_ROOT is not set.", call. = FALSE)
source(file.path(release_root, "analysis", "release_helpers.R"), local = FALSE)

manifest <- read.csv(
  file.path(release_root, "analysis", "config", "primary_run_matrix.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
task <- manifest[manifest$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L) stop("Unknown or duplicate task: ", task_id, call. = FALSE)
dataset_key <- as.character(task$dataset_key)
method_id <- as.character(task$method_id)
method_label <- as.character(task$method)
output_dir <- file.path(release_root, as.character(task$output_dir))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
completion_path <- file.path(output_dir, "complete.json")
if (file.exists(completion_path) && !"--force" %in% args) {
  message(task_id, " already completed.")
  quit(save = "no", status = 0L)
}

started <- Sys.time()
set.seed(as.integer(task$seed))
prepared <- release_load_dataset(dataset_key)
input_hash <- release_validate_input(prepared, dataset_key)

capture <- release_capture(function() release_run_method(
  dataset_key, method_id, prepared
))
writeLines(
  c(
    "WARNINGS", capture$warnings, "", "MESSAGES", capture$messages,
    "", "ERROR", capture$error
  ),
  file.path(output_dir, "conditions.log")
)
if (inherits(capture$value, "error")) {
  jsonlite::write_json(list(
    task_id = task_id, dataset_key = dataset_key, method_id = method_id,
    status = "failed", error = capture$error,
    package = frozen_winn,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ), file.path(output_dir, "failed.json"), auto_unbox = TRUE, pretty = TRUE)
  stop(capture$error, call. = FALSE)
}

matrix <- release_validate_output(capture$value$data, prepared, method_id)
floor_record <- attr(matrix, "intensity_floor", exact = TRUE)
if (!is.null(floor_record)) {
  capture$value$details$intensity_floor <- floor_record
}
matrix_path <- file.path(output_dir, "corrected_matrix.rds")
details_path <- file.path(output_dir, "method_details.rds")
release_atomic_save_rds(matrix, matrix_path)
release_atomic_save_rds(capture$value$details, details_path)
configuration <- release_selected_parameters(dataset_key, method_id)
write.csv(configuration, file.path(output_dir, "selected_configuration.csv"), row.names = FALSE, quote = TRUE)

completed <- Sys.time()
record <- list(
  schema = "endpoint_free_primary_method_v1",
  task_id = task_id,
  dataset_key = dataset_key,
  method_id = method_id,
  method = method_label,
  status = "completed",
  reuse_decision = release_reuse_decision(dataset_key, method_id),
  seed = as.integer(task$seed),
  runtime_sec = capture$runtime_sec,
  wall_sec = as.numeric(difftime(completed, started, units = "secs")),
  n_features = nrow(matrix),
  n_samples = ncol(matrix),
  input_object_sha256 = input_hash,
  output_object_sha256 = release_sha256_object(matrix),
  matrix_file_sha256 = release_sha256_file(matrix_path),
  details_file_sha256 = release_sha256_file(details_path),
  package = frozen_winn,
  started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
  completed_at = format(completed, "%Y-%m-%dT%H:%M:%S%z")
)
jsonlite::write_json(
  record, file.path(output_dir, "run_manifest.json"),
  auto_unbox = TRUE, pretty = TRUE, na = "null"
)
utils::capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
jsonlite::write_json(
  list(task_id = task_id, status = "completed", output_object_sha256 = record$output_object_sha256),
  completion_path, auto_unbox = TRUE, pretty = TRUE
)
message(task_id, " completed: ", dataset_key, " / ", method_label)
