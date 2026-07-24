#!/usr/bin/env Rscript

# Validate and aggregate the 30 independent canonical simulation realizations.
# By default, aggregation is all-or-nothing. Use --allow-incomplete only for an
# explicitly labelled interim diagnostic; manuscript-facing summaries require
# all 30 validated seeds.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("scripts/robustness/aggregate_canonical_simulation_seeds.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
allow_incomplete <- "--allow-incomplete" %in% args
validate_only <- "--validate-only" %in% args
value_arg <- function(prefix, default = NULL) {
  found <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (!length(found)) return(default)
  if (length(found) != 1L) stop("Supply exactly one ", prefix, " argument.", call. = FALSE)
  sub(paste0("^", prefix, "="), "", found[[1L]])
}

project_library <- Sys.getenv("WINN_ROBUSTNESS_R_LIB", unset = "")
if (nzchar(project_library)) {
  if (!dir.exists(project_library)) {
    stop("WINN_ROBUSTNESS_R_LIB does not exist: ", project_library, call. = FALSE)
  }
  .libPaths(c(normalizePath(project_library, mustWork = TRUE), .libPaths()))
}

required_packages <- c("digest", "jsonlite", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

result_root_arg <- value_arg("--result-root")
result_root <- if (is.null(result_root_arg)) {
  file.path(repo_root, "results", "robustness", "04_simulation_seed_stability")
} else {
  normalizePath(result_root_arg, mustWork = TRUE)
}
config_dir <- file.path(result_root, "config")
runs_dir <- file.path(result_root, "runs")
aggregate_dir <- file.path(result_root, "aggregate")
dir.create(aggregate_dir, recursive = TRUE, showWarnings = FALSE)
validation_output_dir <- if (validate_only) file.path(result_root, "validation") else aggregate_dir
dir.create(validation_output_dir, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(value, path) {
  utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
}
sha_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
read_csv <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
scalar_true <- function(value) {
  length(value) == 1L && !is.na(value) && isTRUE(as.logical(value))
}
collapse_reasons <- function(reasons) {
  reasons <- unique(reasons[nzchar(reasons)])
  if (length(reasons)) paste(reasons, collapse = "; ") else ""
}

study_config_path <- file.path(config_dir, "study_config.json")
parameter_hash_path <- file.path(config_dir, "method_parameters.sha256")
method_order_path <- file.path(config_dir, "method_order.csv")
metric_dictionary_path <- file.path(config_dir, "metric_dictionary.csv")
contrast_path <- file.path(config_dir, "primary_contrasts.csv")
config_manifest_path <- file.path(config_dir, "config_checksums.csv")
required_config <- c(
  study_config_path, parameter_hash_path, method_order_path,
  metric_dictionary_path, contrast_path, config_manifest_path
)
if (!all(file.exists(required_config))) {
  stop(
    "Prepared seed-stability configuration is incomplete: ",
    paste(basename(required_config[!file.exists(required_config)]), collapse = ", "),
    call. = FALSE
  )
}

study_config <- jsonlite::read_json(study_config_path, simplifyVector = TRUE)
seed_ids <- as.character(study_config$seed_ids)
method_order <- read_csv(method_order_path)$method
metric_dictionary <- read_csv(metric_dictionary_path)
contrasts <- read_csv(contrast_path)
expected_parameter_hash <- trimws(readLines(parameter_hash_path, warn = FALSE)[[1L]])
bootstrap_seed <- as.integer(study_config$bootstrap$seed)
bootstrap_replicates <- as.integer(study_config$bootstrap$replicates)
equality_tolerance <- as.numeric(study_config$equality_tolerance)

if (!identical(seed_ids, sprintf("SIM%02d", seq_len(30L)))) {
  stop("The prepared study configuration does not contain SIM01--SIM30 in order.", call. = FALSE)
}
if (length(method_order) != 18L || anyDuplicated(method_order)) {
  stop("The method-order configuration must contain 18 unique analysis units.", call. = FALSE)
}
if (nrow(metric_dictionary) != 28L || anyDuplicated(metric_dictionary$metric)) {
  stop("The metric dictionary must contain 28 unique metrics.", call. = FALSE)
}
if (!all(metric_dictionary$direction %in% c("lower", "higher"))) {
  stop("Every metric must declare direction as lower or higher.", call. = FALSE)
}
if (!all(contrasts$method_a %in% method_order & contrasts$method_b %in% method_order)) {
  stop("At least one primary contrast refers to an unknown analysis unit.", call. = FALSE)
}

config_manifest <- read_csv(config_manifest_path)
config_paths <- file.path(config_dir, config_manifest$relative_path)
config_manifest_valid <- all(file.exists(config_paths)) &&
  all(as.numeric(file.info(config_paths)$size) == config_manifest$bytes) &&
  identical(
    unname(vapply(config_paths, sha_file, character(1))),
    unname(config_manifest$sha256)
  )
if (!config_manifest_valid) {
  stop("Prepared configuration failed its file-size or SHA-256 manifest.", call. = FALSE)
}

required_seed_files <- c(
  "analysis_complete.json", "artifact_checksums.csv", "run_manifest.csv",
  "method_status.csv", "errors.csv", "runtime.csv", "method_metrics.csv",
  "validation_checks.csv", "method_matrix_checks.csv", "matrix_checksums.csv",
  "metric_equivalence_checks.csv", "method_exposure_protocol.csv",
  "input_validation.csv", "selected_parameters.csv"
)

seed_tables <- list()
validation_rows <- list()
failure_rows <- list()

validate_one_seed <- function(seed_id) {
  run_dir <- file.path(runs_dir, seed_id)
  reasons <- character()
  required_paths <- file.path(run_dir, required_seed_files)
  missing <- required_seed_files[!file.exists(required_paths)]
  if (length(missing)) {
    reasons <- c(reasons, paste0("missing required artifact(s): ", paste(missing, collapse = ",")))
  }

  tables <- list()
  safe_table <- function(name) {
    path <- file.path(run_dir, name)
    if (!file.exists(path)) return(NULL)
    tryCatch(read_csv(path), error = function(e) {
      reasons <<- c(reasons, paste0("unreadable ", name, ": ", conditionMessage(e)))
      NULL
    })
  }
  tables$manifest <- safe_table("run_manifest.csv")
  tables$status <- safe_table("method_status.csv")
  tables$errors <- safe_table("errors.csv")
  tables$runtime <- safe_table("runtime.csv")
  tables$metrics <- safe_table("method_metrics.csv")
  tables$invariants <- safe_table("validation_checks.csv")
  tables$matrix_checks <- safe_table("method_matrix_checks.csv")
  tables$matrix_manifest <- safe_table("matrix_checksums.csv")
  tables$metric_equivalences <- safe_table("metric_equivalence_checks.csv")
  tables$exposure <- safe_table("method_exposure_protocol.csv")
  tables$input_validation <- safe_table("input_validation.csv")
  tables$selected_parameters <- safe_table("selected_parameters.csv")

  completion <- NULL
  completion_path <- file.path(run_dir, "analysis_complete.json")
  if (file.exists(completion_path)) {
    completion <- tryCatch(
      jsonlite::read_json(completion_path, simplifyVector = TRUE),
      error = function(e) {
        reasons <<- c(reasons, paste0("unreadable analysis_complete.json: ", conditionMessage(e)))
        NULL
      }
    )
  }

  artifact_ok <- FALSE
  artifact_count <- 0L
  artifact_manifest <- safe_table("artifact_checksums.csv")
  if (!is.null(artifact_manifest)) {
    required_columns <- c("relative_path", "bytes", "sha256")
    if (!all(required_columns %in% names(artifact_manifest))) {
      reasons <- c(reasons, "artifact manifest lacks required columns")
    } else {
      artifact_paths <- file.path(run_dir, artifact_manifest$relative_path)
      artifact_count <- nrow(artifact_manifest)
      artifact_ok <- all(file.exists(artifact_paths)) &&
        all(as.numeric(file.info(artifact_paths)$size) == artifact_manifest$bytes) &&
        identical(
          unname(vapply(artifact_paths, sha_file, character(1))),
          unname(artifact_manifest$sha256)
        )
      if (!artifact_ok) reasons <- c(reasons, "artifact checksum validation failed")
    }
  }

  manifest_ok <- FALSE
  if (!is.null(tables$manifest)) {
    manifest <- tables$manifest
    manifest_ok <- nrow(manifest) == 1L &&
      identical(as.character(manifest$seed_id), seed_id) &&
      identical(as.character(manifest$status), "completed") &&
      identical(as.character(manifest$frozen_parameters_sha256), expected_parameter_hash) &&
      as.integer(manifest$n_features) == 1000L &&
      as.integer(manifest$n_samples) == 500L &&
      as.integer(manifest$n_batches) == 5L &&
      as.integer(manifest$n_qc_total) == 50L &&
      as.integer(manifest$n_hidden) == 10L &&
      as.integer(manifest$n_training) == 40L &&
      as.integer(manifest$methods_attempted) == 18L &&
      as.integer(manifest$methods_completed) == 18L &&
      as.integer(manifest$metrics_rows) == 18L * 28L &&
      scalar_true(manifest$invariants_passed) &&
      scalar_true(manifest$metric_equivalences_passed) &&
      scalar_true(manifest$hidden_exclusion_passed)
    if (!manifest_ok) reasons <- c(reasons, "run manifest failed canonical counts or provenance checks")
  }

  completion_ok <- !is.null(completion) &&
    identical(as.character(completion$seed_id), seed_id) &&
    scalar_true(completion$invariants_passed) &&
    scalar_true(completion$metrics_complete) &&
    as.integer(completion$method_count) == 18L &&
    as.integer(completion$metric_rows) == 18L * 28L
  if (!completion_ok) reasons <- c(reasons, "completion marker is absent or invalid")
  if (completion_ok && !is.null(tables$manifest) &&
      !identical(
        as.character(completion$analysis_context_sha256),
        as.character(tables$manifest$analysis_context_sha256)
      )) {
    completion_ok <- FALSE
    reasons <- c(reasons, "completion and run-manifest context hashes differ")
  }

  status_ok <- FALSE
  if (!is.null(tables$status)) {
    status_ok <- nrow(tables$status) == 18L &&
      identical(as.character(tables$status$method), method_order) &&
      all(tables$status$attempted) &&
      all(grepl("^completed", tables$status$status))
    if (!status_ok) reasons <- c(reasons, "not all 18 analysis units completed in canonical order")
  }

  invariant_ok <- FALSE
  if (!is.null(tables$invariants)) {
    required_checks <- c(
      "all_18_analysis_units_attempted", "all_18_analysis_units_completed",
      "raw_equals_c0", "fixed_equals_c4", "c4_equals_gss",
      "c4_equals_independent_current_winn_fixed",
      "all_matrices_exactly_aligned", "same_hidden_references_all_methods",
      "hidden_references_absent_from_controls_and_tuning",
      "single_output_per_method_for_all_metrics", "canonical_bundle_identity_recorded"
    )
    invariant_ok <- identical(as.character(tables$invariants$check), required_checks) &&
      all(tables$invariants$passed)
    equality_rows <- tables$invariants$check %in% c(
      "raw_equals_c0", "fixed_equals_c4", "c4_equals_gss",
      "c4_equals_independent_current_winn_fixed"
    )
    equality_values <- suppressWarnings(as.numeric(tables$invariants$observed[equality_rows]))
    invariant_ok <- invariant_ok && length(equality_values) == 4L &&
      all(is.finite(equality_values)) && all(equality_values <= equality_tolerance)
    if (!invariant_ok) reasons <- c(reasons, "one or more exact seed invariants failed")
  }

  matrix_ok <- FALSE
  if (!is.null(tables$matrix_checks) && !is.null(tables$matrix_manifest)) {
    matrix_columns <- c(
      "numeric_matrix", "exact_dimensions", "exact_feature_order",
      "exact_sample_order", "all_finite", "log1p_valid"
    )
    matrix_paths <- if ("relative_path" %in% names(tables$matrix_manifest)) {
      file.path(run_dir, tables$matrix_manifest$relative_path)
    } else {
      character()
    }
    matrix_artifacts_ok <- length(matrix_paths) == 18L &&
      all(file.exists(matrix_paths)) &&
      all(as.numeric(file.info(matrix_paths)$size) == tables$matrix_manifest$bytes) &&
      identical(
        unname(vapply(matrix_paths, sha_file, character(1))),
        unname(tables$matrix_manifest$sha256)
      )
    matrix_ok <- nrow(tables$matrix_checks) == 18L &&
      identical(as.character(tables$matrix_checks$method), method_order) &&
      all(matrix_columns %in% names(tables$matrix_checks)) &&
      all(as.matrix(tables$matrix_checks[, matrix_columns, drop = FALSE])) &&
      nrow(tables$matrix_manifest) == 18L &&
      setequal(as.character(tables$matrix_manifest$method), method_order) &&
      all(as.integer(tables$matrix_manifest$n_features) == 1000L) &&
      all(as.integer(tables$matrix_manifest$n_samples) == 500L) &&
      matrix_artifacts_ok
    if (!matrix_ok) reasons <- c(reasons, "matrix dimensions, identities, or alignment checks failed")
  }

  metric_equivalence_ok <- !is.null(tables$metric_equivalences) &&
    nrow(tables$metric_equivalences) == 3L &&
    identical(as.character(tables$metric_equivalences$check), c(
      "raw_equals_c0_metrics", "fixed_equals_c4_metrics", "c4_equals_gss_metrics"
    )) &&
    all(tables$metric_equivalences$metrics_compared == 27L) &&
    all(tables$metric_equivalences$same_na_pattern) &&
    all(tables$metric_equivalences$max_abs_difference <= 1e-12) &&
    all(tables$metric_equivalences$passed)
  if (!metric_equivalence_ok) reasons <- c(reasons, "metric-alias equivalence checks failed")

  exposure_ok <- FALSE
  if (!is.null(tables$exposure)) {
    exposure_ok <- nrow(tables$exposure) == 18L &&
      identical(as.character(tables$exposure$method), method_order) &&
      all(as.integer(tables$exposure$n_hidden_controls_supplied) == 0L) &&
      all(!tables$exposure$hidden_labels_visible) &&
      all(!tables$exposure$hidden_used_for_tuning) &&
      all(tables$exposure$one_output_for_all_endpoints) &&
      length(unique(tables$exposure$hidden_ids_hash)) == 1L
    if (!exposure_ok) reasons <- c(reasons, "hidden-reference label-exposure audit failed")
  }

  input_ok <- !is.null(tables$input_validation) &&
    nrow(tables$input_validation) > 0L && all(tables$input_validation$passed)
  if (!input_ok) reasons <- c(reasons, "input validation did not pass")

  metrics_ok <- FALSE
  if (!is.null(tables$metrics)) {
    metric_keys <- paste(tables$metrics$method, tables$metrics$metric, sep = "::")
    expected_keys <- as.vector(outer(method_order, metric_dictionary$metric, paste, sep = "::"))
    directions <- metric_dictionary$direction[match(tables$metrics$metric, metric_dictionary$metric)]
    metrics_ok <- nrow(tables$metrics) == 18L * 28L &&
      !anyDuplicated(metric_keys) && setequal(metric_keys, expected_keys) &&
      identical(as.character(tables$metrics$seed_id), rep(seed_id, nrow(tables$metrics))) &&
      identical(as.character(tables$metrics$metric_direction), directions)
    if (!metrics_ok) reasons <- c(reasons, "method-metric table is incomplete or inconsistent with its dictionary")
  }

  selected_ok <- !is.null(tables$selected_parameters) &&
    nrow(tables$selected_parameters) == 18L &&
    identical(as.character(tables$selected_parameters$method), method_order) &&
    all(!tables$selected_parameters$hidden_references_used)
  if (!selected_ok) reasons <- c(reasons, "selected-parameter audit is incomplete")

  valid <- !length(reasons) && all(c(
    artifact_ok, manifest_ok, completion_ok, status_ok, invariant_ok,
    matrix_ok, metric_equivalence_ok, exposure_ok, input_ok, metrics_ok, selected_ok
  ))
  validation <- data.frame(
    seed_id = seed_id,
    valid_for_aggregation = valid,
    required_files_present = !length(missing),
    artifact_checksums_valid = artifact_ok,
    run_manifest_valid = manifest_ok,
    completion_marker_valid = completion_ok,
    methods_18_of_18_completed = status_ok,
    exact_invariants_passed = invariant_ok,
    matrices_aligned = matrix_ok,
    metric_aliases_equivalent = metric_equivalence_ok,
    hidden_reference_exclusion_passed = exposure_ok,
    inputs_validated = input_ok,
    metrics_504_rows_complete = metrics_ok,
    parameters_audited = selected_ok,
    artifact_count = artifact_count,
    failure_reason = collapse_reasons(reasons),
    stringsAsFactors = FALSE
  )
  list(validation = validation, tables = tables)
}

for (seed_id in seed_ids) {
  audited <- validate_one_seed(seed_id)
  validation_rows[[seed_id]] <- audited$validation
  seed_tables[[seed_id]] <- audited$tables

  status <- audited$tables$status
  if (is.null(status) || !all(c("method", "status") %in% names(status))) {
    failure_rows[[seed_id]] <- data.frame(
      seed_id = seed_id, method = method_order,
      attempted = FALSE, status = "missing_or_unreadable_seed_run",
      error = audited$validation$failure_reason,
      valid_seed = FALSE, stringsAsFactors = FALSE
    )
  } else {
    index <- match(method_order, status$method)
    failure_rows[[seed_id]] <- data.frame(
      seed_id = seed_id, method = method_order,
      attempted = ifelse(is.na(index), FALSE, as.logical(status$attempted[index])),
      status = ifelse(is.na(index), "missing_method_status", as.character(status$status[index])),
      error = ifelse(is.na(index), "method absent from status table", as.character(status$error[index])),
      valid_seed = audited$validation$valid_for_aggregation,
      stringsAsFactors = FALSE
    )
  }
}

seed_validation <- dplyr::bind_rows(validation_rows)
method_failures <- dplyr::bind_rows(failure_rows)
write_csv(seed_validation, file.path(validation_output_dir, "seed_validation.csv"))
write_csv(method_failures, file.path(validation_output_dir, "method_completion_and_failures.csv"))

valid_seeds <- seed_validation$seed_id[seed_validation$valid_for_aggregation]
complete_30 <- identical(valid_seeds, seed_ids)

validation_report <- c(
  "# Canonical repeated-simulation aggregation validation",
  "",
  paste0("- Expected independent seed units: ", length(seed_ids)),
  paste0("- Validated independent seed units: ", length(valid_seeds)),
  paste0("- Complete SIM01--SIM30 panel: ", complete_30),
  paste0("- Interim incomplete aggregation explicitly allowed: ", allow_incomplete),
  paste0("- Validation-only mode: ", validate_only),
  paste0("- Equality tolerance: ", format(equality_tolerance, scientific = TRUE)),
  paste0("- Frozen method-parameter SHA-256: `", expected_parameter_hash, "`"),
  "",
  "Every seed must have 18 completed analysis units, 504 method-metric rows, exact",
  "matrix identities/order, a single label-blinded hidden-reference assignment, and",
  "valid artifact checksums. Invalid seeds are listed in `seed_validation.csv` and",
  "are never silently included.",
  if (length(setdiff(seed_ids, valid_seeds))) c(
    "",
    paste0("Invalid or incomplete seeds: ", paste(setdiff(seed_ids, valid_seeds), collapse = ", "))
  ) else character()
)
writeLines(validation_report, file.path(validation_output_dir, "VALIDATION.md"))

if (validate_only) {
  message(
    "Validation-only audit complete: ", length(valid_seeds), "/30 seeds valid; ",
    "no aggregate estimates were produced. Results: ", validation_output_dir
  )
  quit(save = "no", status = 0L)
}
if (!complete_30 && !allow_incomplete) {
  stop(
    "Aggregation refused: only ", length(valid_seeds), "/30 seed runs passed validation. ",
    "Inspect aggregate/seed_validation.csv. Use --allow-incomplete only for an explicitly interim diagnostic.",
    call. = FALSE
  )
}
if (!length(valid_seeds)) stop("No validated seed is available to aggregate.", call. = FALSE)

bind_valid <- function(name) {
  dplyr::bind_rows(lapply(valid_seeds, function(seed_id) seed_tables[[seed_id]][[name]]))
}
run_manifests <- bind_valid("manifest")
statuses <- bind_valid("status")
errors <- bind_valid("errors")
runtimes <- bind_valid("runtime")
method_metrics <- bind_valid("metrics")
invariants <- bind_valid("invariants")
matrix_checks <- bind_valid("matrix_checks")
matrix_manifests <- bind_valid("matrix_manifest")
metric_equivalences <- bind_valid("metric_equivalences")
exposure <- bind_valid("exposure")
selected_parameters <- bind_valid("selected_parameters")

for (table_name in c(
  "statuses", "errors", "runtimes", "invariants", "matrix_checks",
  "matrix_manifests", "metric_equivalences", "exposure", "selected_parameters"
)) {
  value <- get(table_name)
  if (!"seed_id" %in% names(value)) {
    stop("Combined table ", table_name, " lacks seed_id; per-seed provenance is incomplete.", call. = FALSE)
  }
}

write_csv(run_manifests, file.path(aggregate_dir, "all_run_manifests.csv"))
write_csv(statuses, file.path(aggregate_dir, "all_method_status.csv"))
write_csv(errors, file.path(aggregate_dir, "all_method_errors.csv"))
write_csv(runtimes, file.path(aggregate_dir, "all_runtime.csv"))
write_csv(method_metrics, file.path(aggregate_dir, "all_seed_method_metrics.csv"))
write_csv(invariants, file.path(aggregate_dir, "all_exact_invariants.csv"))
write_csv(matrix_checks, file.path(aggregate_dir, "all_method_matrix_checks.csv"))
write_csv(matrix_manifests, file.path(aggregate_dir, "all_matrix_checksums.csv"))
write_csv(metric_equivalences, file.path(aggregate_dir, "all_metric_equivalence_checks.csv"))
write_csv(exposure, file.path(aggregate_dir, "all_hidden_reference_exposure_audits.csv"))
write_csv(selected_parameters, file.path(aggregate_dir, "all_selected_parameters.csv"))

summary_values <- function(values) {
  finite <- values[is.finite(values)]
  if (!length(finite)) {
    return(c(
      n_finite = 0, mean = NA, sd = NA, median = NA, q25 = NA,
      q75 = NA, iqr = NA, empirical_q025 = NA, empirical_q975 = NA
    ))
  }
  quantiles <- stats::quantile(finite, c(0.025, 0.25, 0.5, 0.75, 0.975), names = FALSE)
  c(
    n_finite = length(finite), mean = mean(finite),
    sd = if (length(finite) > 1L) stats::sd(finite) else NA_real_,
    median = quantiles[[3L]], q25 = quantiles[[2L]], q75 = quantiles[[4L]],
    iqr = quantiles[[4L]] - quantiles[[2L]],
    empirical_q025 = quantiles[[1L]], empirical_q975 = quantiles[[5L]]
  )
}

set.seed(bootstrap_seed)
bootstrap_index_cache <- lapply(seq_len(length(seed_ids)), function(n) {
  matrix(sample.int(n, size = bootstrap_replicates * n, replace = TRUE),
         nrow = bootstrap_replicates, ncol = n)
})

bootstrap_mean_ci <- function(values, replicates) {
  finite <- values[is.finite(values)]
  if (!length(finite)) return(c(bootstrap_mean_q025 = NA_real_, bootstrap_mean_q975 = NA_real_))
  if (length(finite) == 1L) return(c(bootstrap_mean_q025 = finite, bootstrap_mean_q975 = finite))
  if (!identical(as.integer(replicates), bootstrap_replicates)) {
    stop("Bootstrap replicate count differs from the frozen study configuration.", call. = FALSE)
  }
  indices <- bootstrap_index_cache[[length(finite)]]
  estimates <- rowMeans(matrix(finite[indices], nrow = replicates, ncol = length(finite)))
  stats::setNames(
    as.numeric(stats::quantile(estimates, c(0.025, 0.975), names = FALSE)),
    c("bootstrap_mean_q025", "bootstrap_mean_q975")
  )
}

method_metrics$method_order <- match(method_metrics$method, method_order)
method_metrics$metric_order <- match(method_metrics$metric, metric_dictionary$metric)
method_metrics <- method_metrics[order(method_metrics$method_order, method_metrics$metric_order, method_metrics$seed_id), ]
method_groups <- split(
  method_metrics,
  interaction(method_metrics$method_order, method_metrics$metric_order, drop = TRUE, lex.order = TRUE)
)
method_summary <- dplyr::bind_rows(lapply(method_groups, function(group) {
  descriptive <- summary_values(group$value)
  bootstrap <- bootstrap_mean_ci(group$value, bootstrap_replicates)
  data.frame(
    method = group$method[[1L]], metric = group$metric[[1L]],
    metric_direction = group$metric_direction[[1L]], units = group$units[[1L]],
    n_seed_rows = nrow(group), n_finite = unname(descriptive[["n_finite"]]),
    n_missing = sum(!is.finite(group$value)),
    mean = unname(descriptive[["mean"]]), sd = unname(descriptive[["sd"]]),
    median = unname(descriptive[["median"]]), q25 = unname(descriptive[["q25"]]),
    q75 = unname(descriptive[["q75"]]), iqr = unname(descriptive[["iqr"]]),
    empirical_q025 = unname(descriptive[["empirical_q025"]]),
    empirical_q975 = unname(descriptive[["empirical_q975"]]),
    bootstrap_mean_q025 = unname(bootstrap[["bootstrap_mean_q025"]]),
    bootstrap_mean_q975 = unname(bootstrap[["bootstrap_mean_q975"]]),
    bootstrap_seed = bootstrap_seed, bootstrap_replicates = bootstrap_replicates,
    bootstrap_unit = "independent simulation seed",
    complete_30_seed_panel = complete_30,
    stringsAsFactors = FALSE
  )
}))
method_summary$method_order <- match(method_summary$method, method_order)
method_summary$metric_order <- match(method_summary$metric, metric_dictionary$metric)
method_summary <- method_summary[order(method_summary$method_order, method_summary$metric_order), ]
method_summary <- method_summary[, setdiff(names(method_summary), c("method_order", "metric_order"))]
write_csv(method_summary, file.path(aggregate_dir, "method_metric_seed_summary.csv"))

paired_rows <- list()
row_index <- 0L
for (contrast_index in seq_len(nrow(contrasts))) {
  contrast <- contrasts[contrast_index, , drop = FALSE]
  for (metric_index in seq_len(nrow(metric_dictionary))) {
    metric_info <- metric_dictionary[metric_index, , drop = FALSE]
    a <- method_metrics[
      method_metrics$method == contrast$method_a & method_metrics$metric == metric_info$metric,
      c("seed_id", "value"), drop = FALSE
    ]
    b <- method_metrics[
      method_metrics$method == contrast$method_b & method_metrics$metric == metric_info$metric,
      c("seed_id", "value"), drop = FALSE
    ]
    names(a)[[2L]] <- "value_a"
    names(b)[[2L]] <- "value_b"
    paired <- merge(data.frame(seed_id = valid_seeds), a, by = "seed_id", all.x = TRUE, sort = FALSE)
    paired <- merge(paired, b, by = "seed_id", all.x = TRUE, sort = FALSE)
    paired <- paired[match(valid_seeds, paired$seed_id), , drop = FALSE]
    paired$delta_a_minus_b <- paired$value_a - paired$value_b
    paired$complete_pair <- is.finite(paired$value_a) & is.finite(paired$value_b)
    paired$expected_direction_supported <- ifelse(
      !paired$complete_pair, NA,
      if (metric_info$direction == "lower") paired$delta_a_minus_b < 0 else paired$delta_a_minus_b > 0
    )
    paired$tied <- paired$complete_pair & paired$delta_a_minus_b == 0
    paired$contrast_id <- contrast$contrast_id
    paired$method_a <- contrast$method_a
    paired$method_b <- contrast$method_b
    paired$metric <- metric_info$metric
    paired$metric_direction <- metric_info$direction
    paired$units <- metric_info$units
    paired$rationale <- contrast$rationale
    row_index <- row_index + 1L
    paired_rows[[row_index]] <- paired[, c(
      "seed_id", "contrast_id", "method_a", "method_b", "metric",
      "metric_direction", "units", "value_a", "value_b", "delta_a_minus_b",
      "complete_pair", "expected_direction_supported", "tied", "rationale"
    )]
  }
}
paired_contrasts <- dplyr::bind_rows(paired_rows)
write_csv(paired_contrasts, file.path(aggregate_dir, "all_seed_paired_contrasts.csv"))

contrast_groups <- split(
  paired_contrasts,
  interaction(
    match(paired_contrasts$contrast_id, contrasts$contrast_id),
    match(paired_contrasts$metric, metric_dictionary$metric),
    drop = TRUE, lex.order = TRUE
  )
)
contrast_summary <- dplyr::bind_rows(lapply(contrast_groups, function(group) {
  values <- group$delta_a_minus_b[group$complete_pair]
  descriptive <- summary_values(values)
  bootstrap <- bootstrap_mean_ci(values, bootstrap_replicates)
  expected <- group$expected_direction_supported[group$complete_pair]
  data.frame(
    contrast_id = group$contrast_id[[1L]],
    method_a = group$method_a[[1L]], method_b = group$method_b[[1L]],
    metric = group$metric[[1L]], metric_direction = group$metric_direction[[1L]],
    units = group$units[[1L]], n_seed_rows = nrow(group),
    n_complete_pairs = unname(descriptive[["n_finite"]]),
    n_missing_pairs = sum(!group$complete_pair),
    mean_delta_a_minus_b = unname(descriptive[["mean"]]),
    sd_delta = unname(descriptive[["sd"]]),
    median_delta = unname(descriptive[["median"]]),
    q25_delta = unname(descriptive[["q25"]]), q75_delta = unname(descriptive[["q75"]]),
    iqr_delta = unname(descriptive[["iqr"]]),
    empirical_q025_delta = unname(descriptive[["empirical_q025"]]),
    empirical_q975_delta = unname(descriptive[["empirical_q975"]]),
    bootstrap_mean_q025 = unname(bootstrap[["bootstrap_mean_q025"]]),
    bootstrap_mean_q975 = unname(bootstrap[["bootstrap_mean_q975"]]),
    n_expected_direction = sum(expected, na.rm = TRUE),
    proportion_expected_direction = if (length(expected)) mean(expected) else NA_real_,
    n_tied = sum(group$tied[group$complete_pair]),
    bootstrap_seed = bootstrap_seed, bootstrap_replicates = bootstrap_replicates,
    bootstrap_unit = "paired independent simulation-seed delta",
    complete_30_seed_panel = complete_30,
    rationale = group$rationale[[1L]], stringsAsFactors = FALSE
  )
}))
contrast_summary$contrast_order <- match(contrast_summary$contrast_id, contrasts$contrast_id)
contrast_summary$metric_order <- match(contrast_summary$metric, metric_dictionary$metric)
contrast_summary <- contrast_summary[order(contrast_summary$contrast_order, contrast_summary$metric_order), ]
contrast_summary <- contrast_summary[, setdiff(names(contrast_summary), c("contrast_order", "metric_order"))]
write_csv(contrast_summary, file.path(aggregate_dir, "paired_contrast_seed_summary.csv"))

truth_metrics <- c(
  "truth_sample_profile_pearson_mean",
  "truth_metabolite_profile_icc_mean",
  "truth_log1p_rmse"
)
selective_contrasts <- c(
  "selective_drift_vs_forced_drift",
  "selective_batch_vs_forced_batch",
  "selective_both_vs_forced_both"
)
clean_truth_selective <- contrast_summary[
  contrast_summary$contrast_id %in% selective_contrasts &
    contrast_summary$metric %in% truth_metrics,
  , drop = FALSE
]
write_csv(
  clean_truth_selective,
  file.path(aggregate_dir, "selective_vs_forced_clean_truth_recovery.csv")
)

headline_metrics <- c(
  "heldout_qc_cv_mean", "residual_gam_deviance_mean",
  "residual_ljung_box_proportion", "batch_weighted_pc_r2",
  "truth_sample_profile_pearson_mean", "truth_metabolite_profile_icc_mean",
  "truth_log1p_rmse"
)
headline <- contrast_summary[contrast_summary$metric %in% headline_metrics, , drop = FALSE]
write_csv(headline, file.path(aggregate_dir, "headline_claim_stability.csv"))

cross_seed_checks <- data.frame(
  check = c(
    "all_30_seed_ids_validated", "display_seed_is_sim01",
    "one_unique_bundle_config_per_seed", "one_unique_raw_bundle_per_seed",
    "one_unique_truth_bundle_per_seed",
    "one_frozen_parameter_hash", "one_code_bundle_hash",
    "one_winn_commit", "one_winn_worktree_state",
    "all_seed_metric_keys_unique", "all_seed_metric_row_counts_504",
    "all_methods_complete", "all_hidden_labels_excluded",
    "all_exact_invariants_passed", "all_metric_alias_invariants_passed",
    "same_bundle_used_for_all_metrics_within_seed"
  ),
  passed = c(
    complete_30,
    sum(run_manifests$display_seed) == 1L && run_manifests$seed_id[run_manifests$display_seed] == "SIM01",
    length(unique(run_manifests$bundle_config_sha256)) == length(valid_seeds),
    length(unique(run_manifests$bundle_raw_sha256)) == length(valid_seeds),
    length(unique(run_manifests$bundle_truth_sha256)) == length(valid_seeds),
    length(unique(run_manifests$frozen_parameters_sha256)) == 1L &&
      unique(run_manifests$frozen_parameters_sha256) == expected_parameter_hash,
    length(unique(run_manifests$code_bundle_sha256)) == 1L,
    length(unique(run_manifests$winn_commit)) == 1L &&
      !is.na(unique(run_manifests$winn_commit)) && nzchar(unique(run_manifests$winn_commit)),
    length(unique(run_manifests$winn_status_sha256)) == 1L &&
      !is.na(unique(run_manifests$winn_status_sha256)) && nzchar(unique(run_manifests$winn_status_sha256)),
    !anyDuplicated(paste(method_metrics$seed_id, method_metrics$method, method_metrics$metric, sep = "::")),
    all(table(method_metrics$seed_id) == 18L * 28L),
    all(grepl("^completed", statuses$status)),
    all(exposure$n_hidden_controls_supplied == 0L &
          !exposure$hidden_labels_visible & !exposure$hidden_used_for_tuning),
    all(invariants$passed),
    all(metric_equivalences$passed),
    all(vapply(valid_seeds, function(seed_id) {
      length(unique(seed_tables[[seed_id]]$metrics$seed_id)) == 1L &&
        unique(seed_tables[[seed_id]]$metrics$seed_id) == seed_id
    }, logical(1)))
  ),
  observed = c(
    paste(valid_seeds, collapse = ";"),
    paste(run_manifests$seed_id[run_manifests$display_seed], collapse = ";"),
    paste(unique(run_manifests$bundle_config_sha256), collapse = ";"),
    paste(unique(run_manifests$bundle_raw_sha256), collapse = ";"),
    paste(unique(run_manifests$bundle_truth_sha256), collapse = ";"),
    paste(unique(run_manifests$frozen_parameters_sha256), collapse = ";"),
    paste(unique(run_manifests$code_bundle_sha256), collapse = ";"),
    paste(unique(run_manifests$winn_commit), collapse = ";"),
    paste(unique(run_manifests$winn_status_sha256), collapse = ";"),
    nrow(method_metrics), paste(table(method_metrics$seed_id), collapse = ";"),
    paste(table(statuses$status), collapse = ";"),
    sum(exposure$n_hidden_controls_supplied),
    paste(table(invariants$passed), collapse = ";"),
    paste(table(metric_equivalences$passed), collapse = ";"),
    paste(unique(method_metrics$seed_id), collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write_csv(cross_seed_checks, file.path(aggregate_dir, "cross_seed_validation_checks.csv"))

critical_cross_seed_checks <- cross_seed_checks$check != "all_30_seed_ids_validated"
if (!all(cross_seed_checks$passed[critical_cross_seed_checks]) ||
    (!allow_incomplete && !all(cross_seed_checks$passed))) {
  stop(
    "Cross-seed validation failed after metric aggregation: ",
    paste(cross_seed_checks$check[!cross_seed_checks$passed], collapse = ", "),
    call. = FALSE
  )
}

bootstrap_manifest <- data.frame(
  seed = bootstrap_seed,
  replicates = bootstrap_replicates,
  resampling_unit = "independent simulation seed",
  paired_contrasts_resample = "paired seed-level deltas",
  index_generation = "one deterministic common bootstrap-index matrix per observed seed count",
  feature_level_resampling_used = FALSE,
  valid_seed_count = length(valid_seeds),
  complete_30_seed_panel = complete_30,
  stringsAsFactors = FALSE
)
write_csv(bootstrap_manifest, file.path(aggregate_dir, "bootstrap_protocol.csv"))

aggregate_manifest <- list(
  schema = "canonical_seed_stability_aggregate_v1",
  status = if (complete_30) "complete" else "interim_incomplete_explicitly_allowed",
  seed_ids_expected = seed_ids,
  seed_ids_included = valid_seeds,
  invalid_seed_ids = setdiff(seed_ids, valid_seeds),
  complete_30_seed_panel = complete_30,
  method_count = length(method_order),
  metric_count = nrow(metric_dictionary),
  primary_contrast_count = nrow(contrasts),
  paired_metric_rows = nrow(paired_contrasts),
  bootstrap = list(
    seed = bootstrap_seed, replicates = bootstrap_replicates,
    unit = "independent simulation seed",
    paired = TRUE
  ),
  frozen_parameters_sha256 = expected_parameter_hash,
  aggregator_sha256 = sha_file(script_path),
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
jsonlite::write_json(
  aggregate_manifest, file.path(aggregate_dir, "aggregate_manifest.json"),
  auto_unbox = TRUE, pretty = TRUE, na = "null"
)

readme <- c(
  "# Canonical repeated-simulation seed-stability outputs",
  "",
  paste0("This aggregation contains ", length(valid_seeds), " validated independent canonical simulation seeds."),
  if (complete_30) {
    "The prespecified SIM01--SIM30 panel is complete."
  } else {
    "This is an explicitly incomplete interim aggregation and must not support manuscript claims."
  },
  "",
  "`method_metric_seed_summary.csv` reports each method and metric across seeds: mean/SD,",
  "median/IQR, 2.5th--97.5th empirical interval, and a bootstrap confidence interval",
  "for the seed-level mean. `paired_contrast_seed_summary.csv` reports the same summaries",
  "for prespecified paired differences (method A minus method B), plus the strict proportion",
  "of seeds supporting the metric's expected direction. Lower-is-better metrics support A",
  "when the delta is below zero; higher-is-better metrics support A when it is above zero.",
  "Ties are recorded separately.",
  "",
  "Bootstrap resampling uses paired independent simulation seeds, never features. Gate metrics",
  "remain NA for methods without package-native gate diagnostics. Invalid or failed seeds are",
  "reported in `seed_validation.csv` and `method_completion_and_failures.csv`; they are never",
  "silently dropped. No heterogeneous empirical datasets are pooled in this workstream."
)
writeLines(readme, file.path(aggregate_dir, "README.md"))

# Hash aggregate products last; exclude this manifest from its own contents.
artifact_manifest_path <- file.path(aggregate_dir, "aggregate_artifact_checksums.csv")
aggregate_files <- list.files(aggregate_dir, full.names = TRUE, recursive = FALSE)
aggregate_files <- setdiff(aggregate_files, artifact_manifest_path)
aggregate_files <- aggregate_files[!file.info(aggregate_files)$isdir]
aggregate_artifacts <- data.frame(
  relative_path = basename(aggregate_files),
  bytes = file.info(aggregate_files)$size,
  sha256 = vapply(aggregate_files, sha_file, character(1)),
  stringsAsFactors = FALSE
)
write_csv(aggregate_artifacts, artifact_manifest_path)

message(
  "Canonical seed aggregation complete: ", length(valid_seeds), "/30 validated seeds, ",
  nrow(method_summary), " method summaries, ", nrow(contrast_summary),
  " paired-contrast summaries."
)
