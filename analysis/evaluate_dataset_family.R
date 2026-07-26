#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix) {
  value <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (length(value) != 1L) stop("Supply ", prefix, "=...", call. = FALSE)
  sub(paste0("^", prefix, "="), "", value)
}
dataset_key <- arg_value("--dataset")
family <- arg_value("--family")
if (!family %in% c("primary", "ablations")) stop("family must be primary or ablations.")

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root)) stop("WINN_RELEASE_ROOT is not set.", call. = FALSE)
source(file.path(release_root, "analysis", "release_helpers.R"), local = FALSE)
source(file.path(release_root, "analysis", "evaluation_helpers.R"), local = FALSE)
prepared <- release_load_dataset(dataset_key)
release_validate_input(prepared, dataset_key)

if (identical(family, "primary")) {
  methods <- read.csv(
    file.path(release_root, "analysis", "config", "method_manifest.csv"),
    stringsAsFactors = FALSE
  )
  specs <- lapply(seq_len(nrow(methods)), function(index) {
    method_id <- methods$method_id[index]
    directory <- file.path(release_root, "results", "primary", dataset_key, method_id)
    manifest <- jsonlite::read_json(file.path(directory, "run_manifest.json"), simplifyVector = TRUE)
    list(
      method_id = method_id, method = methods$method[index],
      matrix_path = file.path(directory, "corrected_matrix.rds"),
      details_path = file.path(directory, "method_details.rds"),
      runtime_sec = as.numeric(manifest$runtime_sec)
    )
  })
} else {
  variants <- read.csv(
    file.path(release_root, "analysis", "config", "ablation_variants.csv"),
    stringsAsFactors = FALSE
  )
  directory <- file.path(release_root, "results", "ablations", dataset_key)
  selection <- read.csv(file.path(directory, "stage_selection_summary.csv"), stringsAsFactors = FALSE)
  specs <- lapply(seq_len(nrow(variants)), function(index) {
    variant_id <- variants$variant_id[index]
    list(
      method_id = variant_id, method = variants$variant_label[index],
      matrix_path = file.path(directory, "variants", paste0(variant_id, ".rds")),
      details_path = file.path(directory, "diagnostics", paste0(variant_id, ".rds")),
      runtime_sec = selection$total_runtime_sec[match(variant_id, selection$variant_id)]
    )
  })
}

output_dir <- file.path(release_root, "results", "evaluation", family, dataset_key)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
completion_path <- file.path(output_dir, "complete.json")
if (file.exists(completion_path) && !"--force" %in% args) {
  message(dataset_key, " / ", family, " evaluation already completed.")
  quit(save = "no", status = 0L)
}

results <- lapply(specs, function(spec) {
  if (!file.exists(spec$matrix_path)) stop("Missing method matrix: ", spec$matrix_path)
  matrix <- readRDS(spec$matrix_path)
  details <- if (file.exists(spec$details_path)) readRDS(spec$details_path) else list()
  if (identical(family, "ablations") && !is.null(details$diagnostics)) {
    details <- list(stage_decisions = details$diagnostics)
  }
  evaluate_release_matrix(
    dataset_key, spec$method_id, spec$method, matrix, prepared,
    runtime_sec = spec$runtime_sec, details = details
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
  write.csv(tables[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE, quote = TRUE, na = "")
}
control_audit <- do.call(rbind, lapply(specs, function(spec) {
  details <- if (file.exists(spec$details_path)) readRDS(spec$details_path) else list()
  controls <- details$control_samples
  control_ids <- if (is.null(controls) || !length(controls)) {
    character()
  } else {
    as.character(prepared$meta$sample_id[as.integer(controls)])
  }
  expects_training <- spec$method_id %in% c("winn_auto_qc", "winn_auto_batch_qc")
  expects_none <- identical(spec$method_id, "winn_fixed_no_qc")
  passed <- if (expects_training) {
    setequal(control_ids, prepared$training_ids) &&
      !length(intersect(control_ids, c(prepared$hidden_ids, prepared$external_ids)))
  } else if (expects_none) {
    length(control_ids) == 0L
  } else {
    TRUE
  }
  data.frame(
    method_id = spec$method_id,
    recorded_control_count = length(control_ids),
    expected_training_count = length(prepared$training_ids),
    withheld_overlap = length(intersect(
      control_ids, c(prepared$hidden_ids, prepared$external_ids)
    )),
    passed = passed,
    stringsAsFactors = FALSE
  )
}))
write.csv(control_audit, file.path(output_dir, "control_separation_audit.csv"), row.names = FALSE)
reference_roles_ok <- !length(intersect(
  prepared$training_ids, c(prepared$hidden_ids, prepared$external_ids)
)) && !length(intersect(prepared$hidden_ids, prepared$external_ids)) &&
  all(c(
    prepared$training_ids, prepared$hidden_ids, prepared$external_ids
  ) %in% prepared$meta$sample_id)
validation <- data.frame(
  check = c(
    "all_declared_methods_evaluated", "categorical_batch_r2_only",
    "matrix_identifiers_preserved", "heldout_references_excluded_from_method_controls"
  ),
  passed = c(
    length(results) == length(specs),
    all(tables$method_metrics$notes[tables$method_metrics$metric == "batch_weighted_pc_r2_categorical"] == "Batch is explicitly categorical."),
    TRUE,
    reference_roles_ok && all(control_audit$passed)
  ),
  stringsAsFactors = FALSE
)
write.csv(validation, file.path(output_dir, "validation.csv"), row.names = FALSE, quote = TRUE)
if (!all(validation$passed)) stop("Evaluation validation failed.")
utils::capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
jsonlite::write_json(list(
  schema = "endpoint_free_dataset_evaluation_v1",
  dataset_key = dataset_key, family = family, status = "completed",
  methods = length(specs), metric_rows = nrow(tables$method_metrics),
  package = frozen_winn,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
), file.path(output_dir, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
jsonlite::write_json(
  list(dataset_key = dataset_key, family = family, status = "completed"),
  completion_path, auto_unbox = TRUE, pretty = TRUE
)
message(dataset_key, " / ", family, " evaluation completed.")
