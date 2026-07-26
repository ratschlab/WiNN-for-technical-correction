#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
task_value <- grep("^--task-id=", args, value = TRUE)
if (length(task_value) != 1L) stop("Supply --task-id=REF_NNN.", call. = FALSE)
task_id <- sub("^--task-id=", "", task_value)

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root)) stop("WINN_RELEASE_ROOT is not set.", call. = FALSE)
source(file.path(release_root, "analysis", "release_helpers.R"), local = FALSE)

manifest <- read.csv(
  file.path(release_root, "analysis", "config", "reference_split_task_manifest.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
task <- manifest[manifest$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L) stop("Unknown or duplicate task: ", task_id, call. = FALSE)
dataset_key <- as.character(task$dataset_key)
split_id <- as.character(task$split_id)
method_id <- as.character(task$method_id)
method_label <- as.character(task$method)
output_dir <- file.path(release_root, as.character(task$output_dir))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
completion_path <- file.path(output_dir, "complete.json")
if (file.exists(completion_path) && !"--force" %in% args) {
  message(task_id, " already completed.")
  quit(save = "no", status = 0L)
}

prepared <- release_load_dataset(dataset_key)
input_hash <- release_sha256_object(prepared$x)
assignments <- release_reference_assignments(dataset_key)
split_assignment <- assignments[
  assignments$dataset == dataset_key & assignments$seed_id == split_id, , drop = FALSE
]
if (!nrow(split_assignment)) stop("Missing reference assignment for ", dataset_key, " / ", split_id)
training_ids <- as.character(split_assignment$sample_id[split_assignment$assignment == "training"])
hidden_ids <- as.character(split_assignment$sample_id[split_assignment$assignment == "hidden"])
external_ids <- as.character(split_assignment$sample_id[split_assignment$assignment == "external"])
if (!length(external_ids) && length(prepared$external_ids)) external_ids <- prepared$external_ids
if (!all(c(training_ids, hidden_ids, external_ids) %in% colnames(prepared$x)) ||
    length(intersect(training_ids, c(hidden_ids, external_ids))) ||
    length(intersect(hidden_ids, external_ids))) {
  stop("Reference split roles do not align to the canonical matrix.", call. = FALSE)
}
meta_method <- tq_method_metadata(prepared$meta, training_ids)
meta_technical <- tq_technical_metadata(prepared$meta)
write.csv(split_assignment, file.path(output_dir, "reference_assignment.csv"), row.names = FALSE, quote = TRUE)

set.seed(as.integer(split_assignment$seed[[1L]]) + match(method_id, c(
  "raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf",
  "winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc"
)))
started <- Sys.time()

run_empirical_candidate <- function(candidate) {
  if (identical(dataset_key, "batchcorr_set1") && identical(method_id, "qc_rlsc")) {
    return(run_qc_rlsc_id_safe(
      prepared$x, training_ids, meta_method, span = 1, degree = 1L,
      shift_batches = TRUE, min_qcs_per_batch = 4L
    ))
  }
  tq_run_candidate(
    method_label, prepared$x, training_ids, meta_method, meta_technical,
    candidate, dataset = dataset_key
  )
}

candidate_methods <- c("combat", "qc_rlsc", "qc_rfsc", "serrf")
candidate_scores <- data.frame()
selected_details <- list()
capture <- release_capture(function() {
  if (identical(dataset_key, "simulation")) {
    value <- switch(
      method_id,
      raw = prepared$x,
      combat = run_combat(prepared$x, meta_technical, par_prior = TRUE),
      qc_rlsc = run_qc_rlsc_id_safe(
        prepared$x, training_ids, meta_method, span = 0.75, degree = 1L,
        shift_batches = TRUE, min_qcs_per_batch = 4L
      ),
      qc_rfsc = run_qc_rfsc_with_controls(
        prepared$x, training_ids, meta_method, ntree = 500L, coCV = 30, Frule = 0.8
      ),
      tiger = with_no_cluster(run_tiger_all_corrected(
        prepared$x, training_ids, meta_method, use_injection = TRUE,
        mtry_percent = seq(0.2, 0.8, by = 0.2),
        nodesize_percent = seq(0.2, 0.8, by = 0.2),
        ntree = 500L,
        parallel_cores = as.integer(Sys.getenv("WINN_TIGER_CORES", unset = "1"))
      )),
      serrf = run_serrf_all_corrected(prepared$x, training_ids, meta_method, jitter_eps = 0),
      winn_auto_qc = winn::winn(
        prepared$x, batch = prepared$meta$batch, run_order = prepared$meta$run_order,
        control_samples = which(prepared$meta$sample_id %in% training_ids),
        parameters = "auto", ljung_box_fitdf = 0L, return_details = TRUE
      ),
      winn_auto_batch_qc = winn::winn(
        prepared$x, batch = NULL, run_order = prepared$meta$run_order,
        control_samples = which(prepared$meta$sample_id %in% training_ids),
        parameters = "auto", ljung_box_fitdf = 0L, return_details = TRUE
      ),
      winn_fixed_no_qc = winn::winn(
        prepared$x, batch = prepared$meta$batch, run_order = prepared$meta$run_order,
        control_samples = NULL, parameters = "fixed", ljung_box_fitdf = 0L,
        return_details = TRUE
      ),
      stop("Unsupported method: ", method_id)
    )
    if (inherits(value, "winn_result")) {
      selected_details <<- unclass(value)[setdiff(names(value), "data")]
      return(value$data)
    }
    selected_details <<- list(selection = "prespecified_or_native_default")
    return(value)
  }

  if (method_id %in% candidate_methods) {
    registry <- tq_read_registry(source_root, dataset_key, method_label)
    candidate_values <- vector("list", nrow(registry))
    rows <- vector("list", nrow(registry))
    for (index in seq_len(nrow(registry))) {
      candidate <- registry[index, , drop = FALSE]
      fit <- release_capture(function() run_empirical_candidate(candidate))
      if (inherits(fit$value, "error")) {
        rows[[index]] <- data.frame(
          candidate_id = candidate$candidate_id, parameters = candidate$params,
          status = "failed", qc_score = NA_real_, mean_qc_correlation = NA_real_,
          mean_qc_cv = NA_real_, runtime_sec = fit$runtime_sec,
          error = fit$error, stringsAsFactors = FALSE
        )
      } else {
        candidate_values[[index]] <- release_validate_output(fit$value, prepared, method_id)
        score <- training_qc_score(candidate_values[[index]], training_ids)
        rows[[index]] <- data.frame(
          candidate_id = candidate$candidate_id, parameters = candidate$params,
          status = "completed", qc_score = score$qc_score,
          mean_qc_correlation = score$mean_qc_correlation,
          mean_qc_cv = score$mean_qc_cv, runtime_sec = fit$runtime_sec,
          error = "", stringsAsFactors = FALSE
        )
      }
    }
    candidate_scores <<- dplyr::bind_rows(rows)
    valid <- which(candidate_scores$status == "completed" & is.finite(candidate_scores$qc_score))
    if (!length(valid)) stop("No candidate completed with a finite training-QC score.")
    selected_index <- valid[order(-candidate_scores$qc_score[valid], candidate_scores$candidate_id[valid])][1L]
    candidate_scores$selected <- seq_len(nrow(candidate_scores)) == selected_index
    selected_details <<- list(
      selection = "training_qc_only",
      criterion = "mean QC profile correlation minus 0.5 times mean QC CV fraction",
      candidate_id = candidate_scores$candidate_id[selected_index],
      parameters = candidate_scores$parameters[selected_index],
      hidden_qc_values_available = FALSE,
      biological_labels_available = FALSE,
      replicate_identities_available = FALSE
    )
    return(candidate_values[[selected_index]])
  }

  value <- switch(
    method_id,
    raw = prepared$x,
    tiger = with_no_cluster(run_tiger_all_corrected(
      prepared$x, training_ids, meta_method, use_injection = TRUE,
      mtry_percent = seq(0.2, 0.8, by = 0.2),
      nodesize_percent = seq(0.2, 0.8, by = 0.2),
      ntree = 500L,
      parallel_cores = as.integer(Sys.getenv("WINN_TIGER_CORES", unset = "1"))
    )),
    winn_auto_qc = winn::winn(
      prepared$x, batch = prepared$meta$batch, run_order = prepared$meta$run_order,
      control_samples = which(prepared$meta$sample_id %in% training_ids),
      parameters = "auto", ljung_box_fitdf = 0L, return_details = TRUE
    ),
    winn_auto_batch_qc = winn::winn(
      prepared$x, batch = NULL, run_order = prepared$meta$run_order,
      control_samples = which(prepared$meta$sample_id %in% training_ids),
      parameters = "auto", ljung_box_fitdf = 0L, return_details = TRUE
    ),
    winn_fixed_no_qc = winn::winn(
      prepared$x, batch = prepared$meta$batch, run_order = prepared$meta$run_order,
      control_samples = NULL, parameters = "fixed", ljung_box_fitdf = 0L,
      return_details = TRUE
    ),
    stop("Unsupported method: ", method_id)
  )
  if (inherits(value, "winn_result")) {
    selected_details <<- unclass(value)[setdiff(names(value), "data")]
    value$data
  } else {
    selected_details <<- list(selection = "prespecified_or_native_default")
    value
  }
})

writeLines(
  c("WARNINGS", capture$warnings, "", "MESSAGES", capture$messages, "", "ERROR", capture$error),
  file.path(output_dir, "conditions.log")
)
if (inherits(capture$value, "error")) {
  if (nrow(candidate_scores)) write.csv(candidate_scores, file.path(output_dir, "candidate_scores.csv"), row.names = FALSE, quote = TRUE)
  jsonlite::write_json(list(
    task_id = task_id, dataset_key = dataset_key, split_id = split_id,
    method_id = method_id, status = "failed", error = capture$error,
    package = frozen_winn
  ), file.path(output_dir, "failed.json"), auto_unbox = TRUE, pretty = TRUE)
  stop(capture$error, call. = FALSE)
}

matrix <- release_validate_output(capture$value, prepared, method_id)
floor_record <- attr(matrix, "intensity_floor", exact = TRUE)
if (!is.null(floor_record)) {
  selected_details$intensity_floor <- floor_record
}
matrix_path <- file.path(output_dir, "corrected_matrix.rds")
release_atomic_save_rds(matrix, matrix_path)
release_atomic_save_rds(selected_details, file.path(output_dir, "selected_details.rds"))
if (nrow(candidate_scores)) write.csv(candidate_scores, file.path(output_dir, "candidate_scores.csv"), row.names = FALSE, quote = TRUE)
completed <- Sys.time()
jsonlite::write_json(list(
  schema = "endpoint_free_reference_split_method_v1",
  task_id = task_id, dataset_key = dataset_key, split_id = split_id,
  method_id = method_id, method = method_label, status = "completed",
  n_training_qc = length(training_ids), n_hidden_qc = length(hidden_ids),
  n_external_reference = length(external_ids),
  training_qc_ids_sha256 = release_sha256_object(sort(training_ids)),
  hidden_qc_ids_sha256 = release_sha256_object(sort(hidden_ids)),
  input_object_sha256 = input_hash,
  output_object_sha256 = release_sha256_object(matrix),
  matrix_file_sha256 = release_sha256_file(matrix_path),
  selection_rule = as.character(task$selection_rule),
  hidden_qc_values_available_during_selection = FALSE,
  biological_labels_available_during_selection = FALSE,
  replicate_identities_available_during_selection = FALSE,
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
message(task_id, " completed: ", dataset_key, " / ", split_id, " / ", method_label)
