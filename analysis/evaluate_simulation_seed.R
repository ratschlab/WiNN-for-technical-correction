#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
seed_value <- grep("^--seed-id=", args, value = TRUE)
if (length(seed_value) != 1L) stop("Supply --seed-id=SIM01 through SIM10.", call. = FALSE)
seed_id <- sub("^--seed-id=", "", seed_value)

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root)) stop("WINN_RELEASE_ROOT is not set.", call. = FALSE)
source(file.path(release_root, "analysis", "release_helpers.R"), local = FALSE)
source(file.path(release_root, "analysis", "evaluation_helpers.R"), local = FALSE)
prepared <- release_load_simulation(seed_id)
prepared$meta_method <- tq_method_metadata(prepared$meta, prepared$training_ids)
prepared$meta_technical <- tq_technical_metadata(prepared$meta)

output_dir <- file.path(release_root, "results", "evaluation", "simulation_seeds", seed_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
completion_path <- file.path(output_dir, "complete.json")
if (file.exists(completion_path) && !"--force" %in% args) {
  message(seed_id, " evaluation already completed.")
  quit(save = "no", status = 0L)
}

evaluate_specs <- function(specs, family) {
  results <- lapply(specs, function(spec) {
    matrix <- readRDS(spec$matrix_path)
    details <- if (file.exists(spec$details_path)) readRDS(spec$details_path) else list()
    if (identical(family, "ablations") && !is.null(details$diagnostics)) {
      details <- list(stage_decisions = details$diagnostics)
    }
    evaluate_release_matrix(
      "simulation", spec$method_id, spec$method, matrix, prepared,
      runtime_sec = spec$runtime_sec, details = details
    )
  })
  bind_nonempty <- function(field) {
    values <- lapply(results, `[[`, field)
    values <- values[vapply(values, nrow, integer(1)) > 0L]
    if (length(values)) dplyr::bind_rows(values) else data.frame()
  }
  family_dir <- file.path(output_dir, family)
  dir.create(family_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- list(
    method_metrics = bind_nonempty("metrics"),
    qc_cv_by_feature = bind_nonempty("qc"),
    run_order_gam_by_feature_batch = bind_nonempty("gam"),
    ljung_box_by_feature_batch = bind_nonempty("ljung_box")
  )
  for (name in names(tables)) {
    value <- tables[[name]]
    if (nrow(value)) value$seed_id <- seed_id
    write.csv(value, file.path(family_dir, paste0(name, ".csv")), row.names = FALSE, quote = TRUE, na = "")
  }
  length(results)
}

methods <- read.csv(file.path(release_root, "analysis", "config", "method_manifest.csv"), stringsAsFactors = FALSE)
method_specs <- lapply(seq_len(nrow(methods)), function(index) {
  directory <- file.path(release_root, "results", "simulation_seeds", seed_id, methods$method_id[index])
  if (!file.exists(file.path(directory, "complete.json"))) stop("Incomplete simulation method: ", directory)
  manifest <- jsonlite::read_json(file.path(directory, "run_manifest.json"), simplifyVector = TRUE)
  list(
    method_id = methods$method_id[index], method = methods$method[index],
    matrix_path = file.path(directory, "corrected_matrix.rds"),
    details_path = file.path(directory, "method_details.rds"),
    runtime_sec = as.numeric(manifest$runtime_sec)
  )
})
method_count <- evaluate_specs(method_specs, "methods")

variants <- read.csv(file.path(release_root, "analysis", "config", "ablation_variants.csv"), stringsAsFactors = FALSE)
ablation_dir <- file.path(release_root, "results", "simulation_seed_ablations", seed_id)
selection <- read.csv(file.path(ablation_dir, "stage_selection_summary.csv"), stringsAsFactors = FALSE)
ablation_specs <- lapply(seq_len(nrow(variants)), function(index) {
  variant_id <- variants$variant_id[index]
  list(
    method_id = variant_id, method = variants$variant_label[index],
    matrix_path = file.path(ablation_dir, "variants", paste0(variant_id, ".rds")),
    details_path = file.path(ablation_dir, "diagnostics", paste0(variant_id, ".rds")),
    runtime_sec = selection$total_runtime_sec[match(variant_id, selection$variant_id)]
  )
})
ablation_count <- evaluate_specs(ablation_specs, "ablations")

utils::capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
jsonlite::write_json(list(
  schema = "endpoint_free_simulation_seed_evaluation_v1",
  seed_id = seed_id, status = "completed",
  method_count = method_count, ablation_count = ablation_count,
  raw_input_sha256 = release_sha256_object(prepared$x),
  package = frozen_winn,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
), file.path(output_dir, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
jsonlite::write_json(list(seed_id = seed_id, status = "completed"), completion_path, auto_unbox = TRUE, pretty = TRUE)
message(seed_id, " method and ablation evaluation completed.")
