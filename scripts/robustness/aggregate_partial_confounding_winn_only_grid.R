#!/usr/bin/env Rscript

# Aggregate the reviewer-scoped partial-confounding experiment.
# The active analysis contains only Raw and fixed QC-free WiNN for ten seeds.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "scripts/robustness/aggregate_partial_confounding_winn_only_grid.R",
    mustWork = TRUE
  )
}
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

project_library <- Sys.getenv("WINN_ROBUSTNESS_R_LIB", unset = "")
if (nzchar(project_library)) {
  .libPaths(c(normalizePath(project_library, mustWork = TRUE), .libPaths()))
}
if (!requireNamespace("jsonlite", quietly = TRUE) ||
    !requireNamespace("digest", quietly = TRUE)) {
  stop("jsonlite and digest are required.", call. = FALSE)
}
source(file.path(script_dir, "canonical_cache.R"), local = FALSE)

result_root <- file.path(
  repo_root, "results", "robustness", "06_partial_confounding"
)
source_config <- file.path(result_root, "full_grid", "config")
grid_root <- file.path(result_root, "winn_only_grid")
run_root <- file.path(grid_root, "runs")
aggregate_dir <- file.path(grid_root, "aggregate")
dir.create(aggregate_dir, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(value, path) {
  utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
}
sha_file <- canonical_sha256_file
expected_methods <- c("Raw", "WINN default (no QC)")
completion_schema <- "partial_confounding_winn_only_scenario_complete_v1"

scenario_order_full <- read.csv(
  file.path(source_config, "scenario_order.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
scenario_order <- scenario_order_full[
  scenario_order_full$seed_id %in% sprintf("CONF%02d", 1:10), , drop = FALSE
]
scenario_order$array_index_original <- scenario_order$array_index
scenario_order$array_index <- seq_len(nrow(scenario_order))
if (nrow(scenario_order) != 160L ||
    anyDuplicated(scenario_order$scenario_id) ||
    !identical(scenario_order$array_index, seq_len(160L))) {
  stop("The ten-seed, sixteen-scenario plan is not exactly 160 rows.",
       call. = FALSE)
}

parameters_full <- read.csv(
  file.path(source_config, "method_parameters.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
method_parameters <- parameters_full[
  match(expected_methods, parameters_full$method), , drop = FALSE
]
if (!identical(method_parameters$method, expected_methods)) {
  stop("Raw/WiNN parameters are absent from the frozen method ledger.",
       call. = FALSE)
}
metric_dictionary <- read.csv(
  file.path(source_config, "metric_dictionary.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(metric_dictionary) != 37L ||
    anyDuplicated(metric_dictionary$metric)) {
  stop("Frozen metric dictionary is invalid.", call. = FALSE)
}

validate_artifacts <- function(run_dir, artifacts) {
  if (anyDuplicated(artifacts$relative_path) ||
      any(grepl("(^/|(^|/)\\.\\.(/|$))", artifacts$relative_path))) {
    return(FALSE)
  }
  paths <- file.path(run_dir, artifacts$relative_path)
  all(file.exists(paths)) &&
    all(as.numeric(file.info(paths)$size) == artifacts$bytes) &&
    identical(
      unname(vapply(paths, sha_file, character(1))),
      unname(artifacts$sha256)
    )
}

validation_rows <- list()
metrics_rows <- list()
manifest_rows <- list()
runtime_rows <- list()
parameter_rows <- list()
degradation_rows <- list()

for (i in seq_len(nrow(scenario_order))) {
  scenario <- scenario_order[i, , drop = FALSE]
  scenario_id <- scenario$scenario_id[[1L]]
  run_dir <- file.path(run_root, scenario_id)
  completion_path <- file.path(run_dir, "analysis_complete.json")
  artifact_path <- file.path(run_dir, "artifact_checksums.csv")
  completion <- if (file.exists(completion_path)) {
    tryCatch(
      jsonlite::read_json(completion_path, simplifyVector = TRUE),
      error = function(e) NULL
    )
  } else NULL
  artifacts <- if (file.exists(artifact_path)) {
    tryCatch(
      read.csv(artifact_path, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
  } else NULL

  method_status_path <- file.path(run_dir, "method_status.csv")
  metric_path <- file.path(run_dir, "method_metrics.csv")
  manifest_path <- file.path(run_dir, "run_manifest.csv")
  runtime_path <- file.path(run_dir, "runtime.csv")
  selected_path <- file.path(run_dir, "selected_parameters.csv")
  degradation_path <- file.path(run_dir, "degradation_inputs.csv")
  required_paths <- c(
    method_status_path, metric_path, manifest_path, runtime_path,
    selected_path, degradation_path
  )
  completion_valid <- !is.null(completion) &&
    identical(completion$schema, completion_schema) &&
    identical(completion$scenario_id, scenario_id) &&
    identical(completion$seed_id, scenario$seed_id[[1L]]) &&
    isTRUE(completion$invariants_passed) &&
    isTRUE(completion$metrics_complete) &&
    isTRUE(completion$full_gam_complete) &&
    isTRUE(completion$detail_rows_complete) &&
    as.integer(completion$method_count) == 2L &&
    as.integer(completion$metric_rows) == 74L
  artifact_valid <- !is.null(artifacts) &&
    all(c("relative_path", "bytes", "sha256") %in% names(artifacts)) &&
    validate_artifacts(run_dir, artifacts)
  files_present <- all(file.exists(required_paths))

  status_valid <- FALSE
  metrics_valid <- FALSE
  selected_valid <- FALSE
  if (files_present) {
    method_status <- read.csv(
      method_status_path, stringsAsFactors = FALSE, check.names = FALSE
    )
    metric_table <- read.csv(
      metric_path, stringsAsFactors = FALSE, check.names = FALSE
    )
    selected <- read.csv(
      selected_path, stringsAsFactors = FALSE, check.names = FALSE
    )
    status_valid <- nrow(method_status) == 2L &&
      identical(method_status$method, expected_methods) &&
      all(grepl("^completed", method_status$status)) &&
      all(method_status$attempted)
    metrics_valid <- nrow(metric_table) == 74L &&
      identical(unique(metric_table$method), expected_methods) &&
      all(vapply(
        split(metric_table$metric, metric_table$method),
        function(x) setequal(x, metric_dictionary$metric), logical(1)
      ))
    selected_valid <- nrow(selected) == 2L &&
      identical(selected$method, expected_methods) &&
      all(!selected$hidden_references_used) &&
      all(!selected$phenotype_used)
  }
  passed <- completion_valid && artifact_valid && files_present &&
    status_valid && metrics_valid && selected_valid
  validation_rows[[scenario_id]] <- data.frame(
    array_index = i, scenario_id = scenario_id,
    seed_id = scenario$seed_id[[1L]],
    completion_valid = completion_valid,
    artifact_checksums_valid = artifact_valid,
    required_files_present = files_present,
    method_status_valid = status_valid,
    method_metrics_valid = metrics_valid,
    selected_parameters_valid = selected_valid,
    passed = passed,
    stringsAsFactors = FALSE
  )
  if (!passed) next

  metric_table$confounding_axis <- scenario$confounding_axis[[1L]]
  metric_table$nominal_strength <- scenario$nominal_strength[[1L]]
  metric_table$realized_batch_cramers_v <-
    scenario$cramers_v_phenotype_batch[[1L]]
  metric_table$realized_order_abs_cor_pooled <-
    scenario$abs_cor_phenotype_within_batch_order[[1L]]
  metric_table$realized_order_weighted_mean_abs_cor <-
    scenario$weighted_mean_abs_within_plate_cor[[1L]]
  metric_table$realized_order_max_abs_cor <-
    scenario$max_abs_within_plate_cor[[1L]]
  metrics_rows[[scenario_id]] <- metric_table
  manifest_rows[[scenario_id]] <- read.csv(
    manifest_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  runtime_rows[[scenario_id]] <- read.csv(
    runtime_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  parameter_rows[[scenario_id]] <- selected
  degradation_rows[[scenario_id]] <- read.csv(
    degradation_path, stringsAsFactors = FALSE, check.names = FALSE
  )
}

validation <- do.call(rbind, validation_rows)
rownames(validation) <- NULL
write_csv(validation, file.path(aggregate_dir, "scenario_validation.csv"))
if (nrow(validation) != 160L || !all(validation$passed)) {
  stop(
    "WiNN-only aggregation stopped: ", sum(!validation$passed),
    " of 160 scenario runs failed validation.", call. = FALSE
  )
}

method_metrics <- do.call(rbind, metrics_rows)
run_manifest <- do.call(rbind, manifest_rows)
runtime <- do.call(rbind, runtime_rows)
selected_parameters <- do.call(rbind, parameter_rows)
degradation_inputs <- do.call(rbind, degradation_rows)
rownames(method_metrics) <- NULL
rownames(run_manifest) <- NULL
rownames(runtime) <- NULL
rownames(selected_parameters) <- NULL
rownames(degradation_inputs) <- NULL
run_manifest_valid <- nrow(run_manifest) == 160L &&
  all(run_manifest$status == "completed") &&
  all(run_manifest$methods_attempted == 2L) &&
  all(run_manifest$methods_completed == 2L) &&
  all(run_manifest$metrics_rows == 74L) &&
  all(run_manifest$invariants_passed) &&
  all(run_manifest$full_gam_complete) &&
  all(run_manifest$detail_rows_complete) &&
  all(run_manifest$hidden_exclusion_passed) &&
  all(run_manifest$phenotype_blinding_passed) &&
  length(unique(run_manifest$code_bundle_sha256)) == 1L &&
  length(unique(run_manifest$winn_commit)) == 1L &&
  length(unique(run_manifest$frozen_parameters_sha256)) == 1L &&
  !anyDuplicated(run_manifest$scenario_id)
if (!run_manifest_valid) {
  stop("Combined run manifests failed scope, blinding, or source-identity checks.",
       call. = FALSE)
}

write_csv(scenario_order, file.path(aggregate_dir, "scenario_order_10_seeds.csv"))
write_csv(method_parameters, file.path(aggregate_dir, "method_parameters.csv"))
write_csv(metric_dictionary, file.path(aggregate_dir, "metric_dictionary.csv"))
write_csv(method_metrics, file.path(aggregate_dir, "method_metrics_all.csv"))
write_csv(run_manifest, file.path(aggregate_dir, "run_manifest_all.csv"))
write_csv(runtime, file.path(aggregate_dir, "runtime_all.csv"))
write_csv(selected_parameters,
          file.path(aggregate_dir, "selected_parameters_all.csv"))
write_csv(degradation_inputs,
          file.path(aggregate_dir, "degradation_inputs_all.csv"))

scenario_design <- unique(method_metrics[, c(
  "scenario_id", "seed_id", "confounding_axis", "nominal_strength",
  "realized_batch_cramers_v", "realized_order_abs_cor_pooled",
  "realized_order_weighted_mean_abs_cor", "realized_order_max_abs_cor"
), drop = FALSE])
if (nrow(scenario_design) != 160L) {
  stop("Scenario-design source table is not exactly 160 rows.", call. = FALSE)
}
write_csv(scenario_design,
          file.path(aggregate_dir, "scenario_design_realized.csv"))

summarize_values <- function(table, group_columns, value_column = "value") {
  groups <- split(table, interaction(table[group_columns], drop = TRUE,
                                     lex.order = TRUE))
  result <- lapply(groups, function(group) {
    values <- group[[value_column]]
    values <- values[is.finite(values)]
    output <- group[1L, group_columns, drop = FALSE]
    output$n <- length(values)
    output$mean <- if (length(values)) mean(values) else NA_real_
    output$sd <- if (length(values) > 1L) stats::sd(values) else NA_real_
    output$median <- if (length(values)) stats::median(values) else NA_real_
    output$q025 <- if (length(values)) {
      unname(stats::quantile(values, 0.025, names = FALSE))
    } else NA_real_
    output$q975 <- if (length(values)) {
      unname(stats::quantile(values, 0.975, names = FALSE))
    } else NA_real_
    output
  })
  output <- do.call(rbind, result)
  rownames(output) <- NULL
  output
}

curve_summary <- summarize_values(
  method_metrics,
  c("confounding_axis", "nominal_strength", "method", "metric",
    "metric_direction", "units")
)
write_csv(curve_summary, file.path(aggregate_dir, "curve_summary.csv"))

raw <- method_metrics[method_metrics$method == "Raw", , drop = FALSE]
winn <- method_metrics[
  method_metrics$method == "WINN default (no QC)", , drop = FALSE
]
key_columns <- c("scenario_id", "seed_id", "metric")
raw_key <- do.call(paste, c(raw[key_columns], sep = "\r"))
winn_key <- do.call(paste, c(winn[key_columns], sep = "\r"))
match_raw <- match(winn_key, raw_key)
if (anyNA(match_raw) || length(unique(match_raw)) != nrow(raw)) {
  stop("Raw and WiNN metric rows do not form exact scenario-metric pairs.",
       call. = FALSE)
}
paired <- winn[, c(
  "scenario_id", "seed_id", "confounding_axis", "nominal_strength",
  "realized_batch_cramers_v", "realized_order_abs_cor_pooled",
  "realized_order_weighted_mean_abs_cor", "realized_order_max_abs_cor",
  "metric", "metric_direction", "units"
), drop = FALSE]
paired$raw_value <- raw$value[match_raw]
paired$winn_value <- winn$value
paired$winn_minus_raw <- paired$winn_value - paired$raw_value
paired$favorable_difference <- ifelse(
  paired$metric_direction == "higher", paired$winn_minus_raw,
  ifelse(paired$metric_direction == "lower", -paired$winn_minus_raw, NA_real_)
)
write_csv(paired, file.path(aggregate_dir, "paired_winn_vs_raw.csv"))

paired_summary <- summarize_values(
  paired,
  c("confounding_axis", "nominal_strength", "metric",
    "metric_direction", "units"),
  value_column = "winn_minus_raw"
)
names(paired_summary)[names(paired_summary) %in%
                        c("mean", "sd", "median", "q025", "q975")] <-
  paste0("winn_minus_raw_", names(paired_summary)[names(paired_summary) %in%
                                                    c("mean", "sd", "median", "q025", "q975")])
write_csv(paired_summary,
          file.path(aggregate_dir, "paired_winn_vs_raw_summary.csv"))

paired_favorable <- paired[is.finite(paired$favorable_difference),
                           , drop = FALSE]
paired_favorable_summary <- summarize_values(
  paired_favorable,
  c("confounding_axis", "nominal_strength", "metric",
    "metric_direction", "units"),
  value_column = "favorable_difference"
)
names(paired_favorable_summary)[names(paired_favorable_summary) %in%
                                 c("mean", "sd", "median", "q025", "q975")] <-
  paste0("favorable_", names(paired_favorable_summary)[
    names(paired_favorable_summary) %in%
      c("mean", "sd", "median", "q025", "q975")
  ])
write_csv(
  paired_favorable_summary,
  file.path(aggregate_dir, "paired_favorable_difference_summary.csv")
)

primary_metrics <- c(
  "median_attenuation_ratio_clean_responsive",
  "effect_rmse_vs_clean_responsive", "true_positive_rate",
  "heldout_qc_cv_mean", "batch_weighted_pc_r2_categorical",
  "phenotype_weighted_pc_r2_categorical",
  "residual_ljung_box_proportion_significant",
  "residual_run_order_gam_mean_deviance", "runtime_sec"
)
primary <- paired[paired$metric %in% primary_metrics, , drop = FALSE]
write_csv(primary, file.path(aggregate_dir, "primary_winn_vs_raw.csv"))

utils::capture.output(sessionInfo(),
                      file = file.path(aggregate_dir, "sessionInfo.txt"))
aggregate_manifest <- list(
  schema = "partial_confounding_winn_only_aggregate_v1",
  scenario_count = 160L,
  seed_count = 10L,
  scenario_count_per_seed = 16L,
  methods = expected_methods,
  corrected_method = "WINN default (no QC)",
  raw_role = "uncorrected within-scenario baseline only",
  metric_count = 37L,
  all_scenarios_valid = TRUE,
  competitors_run = FALSE,
  auto_modes_run = FALSE,
  ablations_run = FALSE,
  source_scenario_order_sha256 = sha_file(
    file.path(source_config, "scenario_order.csv")
  ),
  source_method_parameters_sha256 = sha_file(
    file.path(source_config, "method_parameters.csv")
  ),
  aggregate_script_sha256 = sha_file(script_path),
  run_code_bundle_sha256 = unique(run_manifest$code_bundle_sha256),
  winn_commit = unique(run_manifest$winn_commit),
  frozen_selected_method_parameters_sha256 =
    unique(run_manifest$frozen_parameters_sha256),
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
jsonlite::write_json(
  aggregate_manifest, file.path(aggregate_dir, "aggregate_manifest.json"),
  auto_unbox = TRUE, pretty = TRUE, na = "null"
)

artifact_path <- file.path(aggregate_dir, "artifact_checksums.csv")
files <- list.files(aggregate_dir, recursive = TRUE, full.names = TRUE)
files <- setdiff(files, artifact_path)
files <- files[file.info(files)$isdir %in% FALSE]
artifacts <- data.frame(
  relative_path = substring(files, nchar(aggregate_dir) + 2L),
  bytes = as.numeric(file.info(files)$size),
  sha256 = vapply(files, sha_file, character(1)),
  stringsAsFactors = FALSE
)
write_csv(artifacts, artifact_path)
message("Aggregated and validated 160 Raw/WiNN partial-confounding runs.")
