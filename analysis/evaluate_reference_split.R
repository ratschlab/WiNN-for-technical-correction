#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix) {
  value <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (length(value) != 1L) stop("Supply ", prefix, "=...", call. = FALSE)
  sub(paste0("^", prefix, "="), "", value)
}
dataset_key <- arg_value("--dataset")
split_id <- arg_value("--split-id")

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root)) stop("WINN_RELEASE_ROOT is not set.", call. = FALSE)
source(file.path(release_root, "analysis", "release_helpers.R"), local = FALSE)
source(file.path(release_root, "analysis", "evaluation_helpers.R"), local = FALSE)
prepared <- release_load_dataset(dataset_key)
release_validate_input(prepared, dataset_key)

assignments <- release_reference_assignments(dataset_key)
assignment <- assignments[assignments$dataset == dataset_key & assignments$seed_id == split_id, , drop = FALSE]
if (!nrow(assignment)) stop("Missing assignment for ", dataset_key, " / ", split_id)
prepared$training_ids <- as.character(assignment$sample_id[assignment$assignment == "training"])
prepared$hidden_ids <- as.character(assignment$sample_id[assignment$assignment == "hidden"])
# The rotating internal references are the blinded technical endpoint for split
# stability. Fixed external clinical references remain excluded here.
prepared$external_ids <- character()

methods <- read.csv(
  file.path(release_root, "analysis", "config", "method_manifest.csv"),
  stringsAsFactors = FALSE
)
output_dir <- file.path(release_root, "results", "evaluation", "reference_splits", dataset_key, split_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
completion_path <- file.path(output_dir, "complete.json")
if (file.exists(completion_path) && !"--force" %in% args) {
  message(dataset_key, " / ", split_id, " evaluation already completed.")
  quit(save = "no", status = 0L)
}

results <- lapply(seq_len(nrow(methods)), function(index) {
  method_id <- methods$method_id[index]
  directory <- file.path(release_root, "results", "reference_splits", dataset_key, split_id, method_id)
  if (!file.exists(file.path(directory, "complete.json"))) stop("Incomplete reference method: ", directory)
  matrix <- readRDS(file.path(directory, "corrected_matrix.rds"))
  details <- if (file.exists(file.path(directory, "selected_details.rds"))) readRDS(file.path(directory, "selected_details.rds")) else list()
  manifest <- jsonlite::read_json(file.path(directory, "run_manifest.json"), simplifyVector = TRUE)
  evaluate_release_matrix(
    dataset_key, method_id, methods$method[index], matrix, prepared,
    runtime_sec = as.numeric(manifest$runtime_sec), details = details
  )
})

bind_nonempty <- function(field) {
  values <- lapply(results, `[[`, field)
  values <- values[vapply(values, nrow, integer(1)) > 0L]
  if (length(values)) dplyr::bind_rows(values) else data.frame()
}
tables <- list(
  method_metrics = bind_nonempty("metrics"),
  qc_cv_by_feature = bind_nonempty("qc"),
  run_order_gam_by_feature_batch = bind_nonempty("gam"),
  ljung_box_by_feature_batch = bind_nonempty("ljung_box"),
  biological_associations = bind_nonempty("associations"),
  replicate_agreement = bind_nonempty("replicates")
)
for (name in names(tables)) {
  value <- tables[[name]]
  if (nrow(value)) value$split_id <- split_id
  write.csv(value, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE, quote = TRUE, na = "")
}
utils::capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
jsonlite::write_json(list(
  schema = "endpoint_free_reference_split_evaluation_v1",
  dataset_key = dataset_key, split_id = split_id, status = "completed",
  method_count = length(results), hidden_reference_count = length(prepared$hidden_ids),
  package = frozen_winn,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
), file.path(output_dir, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
jsonlite::write_json(
  list(dataset_key = dataset_key, split_id = split_id, status = "completed"),
  completion_path, auto_unbox = TRUE, pretty = TRUE
)
message(dataset_key, " / ", split_id, " reference evaluation completed.")
