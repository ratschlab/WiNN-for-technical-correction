#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
task_value <- grep("^--task-id=", args, value = TRUE)
if (length(task_value) != 1L) stop("Supply --task-id=SIM_NNN.", call. = FALSE)
task_id <- sub("^--task-id=", "", task_value)

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root)) stop("WINN_RELEASE_ROOT is not set.", call. = FALSE)
source(file.path(release_root, "analysis", "release_helpers.R"), local = FALSE)

manifest <- read.csv(
  file.path(release_root, "analysis", "config", "simulation_task_manifest.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
task <- manifest[manifest$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L) stop("Unknown or duplicate task: ", task_id, call. = FALSE)
seed_id <- as.character(task$seed_id)
method_id <- as.character(task$method_id)
method_label <- as.character(task$method)
output_dir <- file.path(release_root, as.character(task$output_dir))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
completion_path <- file.path(output_dir, "complete.json")
if (file.exists(completion_path) && !"--force" %in% args) {
  message(task_id, " already completed.")
  quit(save = "no", status = 0L)
}

prepared <- release_load_simulation(seed_id)
prepared$meta_method <- tq_method_metadata(prepared$meta, prepared$training_ids)
prepared$meta_technical <- tq_technical_metadata(prepared$meta)
input_hash <- release_sha256_object(prepared$x)
set.seed(2026072610L + as.integer(sub("SIM", "", seed_id)) * 100L + match(method_id, c(
  "raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf",
  "winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc"
)))
started <- Sys.time()

capture <- release_capture(function() {
  release_run_method("simulation", method_id, prepared)
})
writeLines(
  c("WARNINGS", capture$warnings, "", "MESSAGES", capture$messages, "", "ERROR", capture$error),
  file.path(output_dir, "conditions.log")
)
if (inherits(capture$value, "error")) {
  jsonlite::write_json(list(
    task_id = task_id, seed_id = seed_id, method_id = method_id,
    status = "failed", error = capture$error, package = frozen_winn
  ), file.path(output_dir, "failed.json"), auto_unbox = TRUE, pretty = TRUE)
  stop(capture$error, call. = FALSE)
}

matrix <- release_validate_output(capture$value$data, prepared, method_id)
floor_record <- attr(matrix, "intensity_floor", exact = TRUE)
if (!is.null(floor_record)) {
  capture$value$details$intensity_floor <- floor_record
}
matrix_path <- file.path(output_dir, "corrected_matrix.rds")
release_atomic_save_rds(matrix, matrix_path)
release_atomic_save_rds(capture$value$details, file.path(output_dir, "method_details.rds"))
completed <- Sys.time()
jsonlite::write_json(list(
  schema = "endpoint_free_simulation_method_v1",
  task_id = task_id, seed_id = seed_id, method_id = method_id,
  method = method_label, status = "completed",
  decision = as.character(task$public_reproduction_decision),
  original_release_decision = as.character(task$decision),
  input_object_sha256 = input_hash,
  output_object_sha256 = release_sha256_object(matrix),
  matrix_file_sha256 = release_sha256_file(matrix_path),
  runtime_sec = capture$runtime_sec,
  wall_sec = as.numeric(difftime(completed, started, units = "secs")),
  package = frozen_winn,
  completed_at = format(completed, "%Y-%m-%dT%H:%M:%S%z")
), file.path(output_dir, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
utils::capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
jsonlite::write_json(
  list(task_id = task_id, status = "completed", output_object_sha256 = release_sha256_object(matrix)),
  completion_path, auto_unbox = TRUE, pretty = TRUE
)
message(task_id, " completed: ", seed_id, " / ", method_label)
