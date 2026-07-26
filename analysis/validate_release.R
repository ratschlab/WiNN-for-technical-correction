#!/usr/bin/env Rscript

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root) || !dir.exists(release_root)) {
  stop("WINN_RELEASE_ROOT must identify the release directory.", call. = FALSE)
}
if (!requireNamespace("jsonlite", quietly = TRUE) ||
    !requireNamespace("digest", quietly = TRUE)) {
  stop("jsonlite and digest are required.", call. = FALSE)
}

result_root <- file.path(release_root, "results")
validation_root <- file.path(result_root, "final", "validation")
dir.create(validation_root, recursive = TRUE, showWarnings = FALSE)

count_files <- function(path, pattern) {
  if (!dir.exists(path)) return(0L)
  length(list.files(path, pattern = pattern, recursive = TRUE, full.names = TRUE))
}
add_check <- local({
  rows <- list()
  function(check = NULL, passed = NULL, detail = NULL, finish = FALSE) {
    if (isTRUE(finish)) return(do.call(rbind, rows))
    rows[[length(rows) + 1L]] <<- data.frame(
      check = check, passed = isTRUE(passed), detail = as.character(detail),
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }
})

families <- data.frame(
  family = c(
    "primary_methods", "dataset_ablation_bundles", "simulation_seed_methods",
    "simulation_seed_ablation_bundles", "reference_split_methods",
    "partial_confounding", "primary_evaluations", "ablation_evaluations",
    "simulation_seed_evaluations", "reference_split_evaluations"
  ),
  path = c(
    "primary", "ablations", "simulation_seeds", "simulation_seed_ablations",
    "reference_splits", "partial_confounding", "evaluation/primary",
    "evaluation/ablations", "evaluation/simulation_seeds",
    "evaluation/reference_splits"
  ),
  expected = c(54L, 6L, 90L, 10L, 540L, 160L, 6L, 6L, 10L, 60L),
  stringsAsFactors = FALSE
)
families$observed <- vapply(families$path, function(path) {
  count_files(file.path(result_root, path), "^complete[.]json$")
}, integer(1))
families$complete <- families$observed == families$expected
write.csv(families, file.path(validation_root, "task_completeness.csv"), row.names = FALSE)
add_check(
  "all_prespecified_tasks_complete", all(families$complete),
  paste(paste0(families$family, "=", families$observed, "/", families$expected), collapse = "; ")
)

logical_ledger_path <- file.path(validation_root, "logical_task_ledger.csv")
execution_ledger_path <- file.path(validation_root, "execution_ledger.csv")
logical_ledger <- if (file.exists(logical_ledger_path)) {
  read.csv(logical_ledger_path, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame()
execution_ledger <- if (file.exists(execution_ledger_path)) {
  read.csv(execution_ledger_path, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame()
ledger_ok <-
  nrow(logical_ledger) == 942L &&
  !anyDuplicated(logical_ledger$task_id) &&
  all(logical_ledger$complete_marker_present) &&
  !any(logical_ledger$family == "unclassified") &&
  nrow(execution_ledger) > 0L &&
  !any(execution_ledger$family == "unclassified") &&
  !any(execution_ledger$status == "running")
add_check(
  "execution_ledgers_cover_prespecified_tasks_and_attempts",
  ledger_ok,
  paste(nrow(logical_ledger), "logical tasks and", nrow(execution_ledger), "scheduler attempts audited")
)

failure_paths <- list.files(result_root, pattern = "^failed[.]json$", recursive = TRUE, full.names = TRUE)
unresolved_failure <- vapply(failure_paths, function(path) {
  !file.exists(file.path(dirname(path), "complete.json"))
}, logical(1))
failure_table <- data.frame(
  path = substring(failure_paths, nchar(release_root) + 2L),
  resolved_by_complete_output = !unresolved_failure,
  stringsAsFactors = FALSE
)
write.csv(failure_table, file.path(validation_root, "failure_artifact_audit.csv"), row.names = FALSE)
add_check(
  "no_unresolved_failure_artifacts", !any(unresolved_failure),
  paste0(sum(!unresolved_failure), " resolved attempt artifact(s); ", sum(unresolved_failure), " unresolved")
)

expected_sha <- "71a0964cee2778b2e5789d20621147e074c7945e813cf76af2ceeb696104aae1"
manifest_paths <- list.files(
  result_root, pattern = "^run_manifest[.]json$", recursive = TRUE, full.names = TRUE
)
package_rows <- lapply(manifest_paths, function(path) {
  value <- tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
  package <- if (is.list(value)) value$package else NULL
  sha <- if (is.list(package)) as.character(package$source_archive_sha256) else NA_character_
  version <- if (is.list(package)) as.character(package$version) else NA_character_
  data.frame(
    path = substring(path, nchar(release_root) + 2L),
    version = if (length(version)) version else NA_character_,
    source_archive_sha256 = if (length(sha)) sha else NA_character_,
    stringsAsFactors = FALSE
  )
})
package_table <- do.call(rbind, package_rows)
analytical <- grepl(
  "results/(primary|ablations|simulation_seeds|simulation_seed_ablations|reference_splits|partial_confounding)/",
  package_table$path
)
package_table$is_frozen_build <- !analytical |
  (!is.na(package_table$version) & !is.na(package_table$source_archive_sha256) &
     package_table$version == "0.1.4" & package_table$source_archive_sha256 == expected_sha)
write.csv(package_table, file.path(validation_root, "package_provenance_audit.csv"), row.names = FALSE, na = "")
add_check(
  "all_analytical_runs_use_frozen_winn", all(package_table$is_frozen_build),
  paste(sum(analytical), "analytical run manifests audited")
)

matrix_families <- c("primary", "reference_splits", "simulation_seeds")
matrix_paths <- unlist(lapply(matrix_families, function(family) {
  list.files(
    file.path(result_root, family), pattern = "^corrected_matrix[.]rds$",
    recursive = TRUE, full.names = TRUE
  )
}), use.names = FALSE)
matrix_domain_rows <- vector("list", length(matrix_paths))
for (index in seq_along(matrix_paths)) {
  value <- readRDS(matrix_paths[index])
  matrix_domain_rows[[index]] <- data.frame(
    path = substring(matrix_paths[index], nchar(release_root) + 2L),
    minimum = min(value),
    negative_values = sum(value < 0),
    nonfinite_values = sum(!is.finite(value)),
    stringsAsFactors = FALSE
  )
  rm(value)
  if (index %% 25L == 0L) gc(FALSE)
}
matrix_domain_audit <- do.call(rbind, matrix_domain_rows)
write.csv(
  matrix_domain_audit,
  file.path(validation_root, "corrected_matrix_domain_audit.csv"),
  row.names = FALSE
)
matrix_domain_ok <-
  nrow(matrix_domain_audit) == 684L &&
  all(matrix_domain_audit$negative_values == 0L) &&
  all(matrix_domain_audit$nonfinite_values == 0L)
add_check(
  "all_corrected_matrices_are_finite_and_nonnegative",
  matrix_domain_ok,
  paste(nrow(matrix_domain_audit), "primary, split, and repeated-simulation matrices audited")
)

flooring_path <- file.path(
  result_root, "final", "tables", "intensity_domain_safeguard.csv"
)
flooring <- if (file.exists(flooring_path)) {
  read.csv(flooring_path, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame()
floor_applied <- if ("floor_applied" %in% names(flooring)) {
  if (is.logical(flooring$floor_applied)) flooring$floor_applied else
    tolower(as.character(flooring$floor_applied)) == "true"
} else logical()
flooring_ok <-
  nrow(flooring) == 684L && length(floor_applied) == 684L &&
  all(flooring$negative_values_floored[!floor_applied] == 0L) &&
  all(flooring$features_with_flooring[!floor_applied] == 0L) &&
  all(flooring$samples_with_flooring[!floor_applied] == 0L) &&
  all(flooring$negative_values_floored[floor_applied] > 0L) &&
  all(flooring$features_with_flooring[floor_applied] > 0L) &&
  all(flooring$samples_with_flooring[floor_applied] > 0L) &&
  all(flooring$minimum_before_floor[floor_applied] < 0)
add_check(
  "negative_competitor_extrapolations_are_recorded_before_flooring",
  flooring_ok,
  paste(sum(floor_applied), "of", nrow(flooring), "correction runs used the intensity-domain floor")
)

selection <- read.csv(
  file.path(release_root, "analysis", "config", "endpoint_free_selection_manifest.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
unavailable_fields <- c(
  "hidden_qc_values_available", "biological_labels_available",
  "replicate_identities_available", "preservation_endpoints_available",
  "batch_r2_available", "final_rank_available"
)
selection_availability_ok <- all(unavailable_fields %in% names(selection)) &&
  all(vapply(unavailable_fields, function(field) {
    values <- selection[[field]]
    values <- if (is.logical(values)) values else tolower(as.character(values)) == "true"
    !any(values, na.rm = TRUE)
  }, logical(1)))
forbidden_terms <- c(
  "held-out qc", "biological r2", "biological association", "replicate pearson",
  "metabolite icc", "preservation endpoint", "final benchmark rank"
)
criterion_text <- tolower(paste(selection$selection_criterion, collapse = " "))
criterion_ok <- !any(vapply(
  forbidden_terms, function(term) grepl(term, criterion_text, fixed = TRUE), logical(1)
))
add_check(
  "selection_manifest_excludes_evaluation_endpoints",
  selection_availability_ok && criterion_ok,
  paste(length(unavailable_fields), "endpoint-availability fields are uniformly false")
)

selection_provenance_fields <- c(
  "selection_seed", "source_package", "source_version",
  "information_available_during_selection", "departure_from_native_default",
  "departure_justification", "configuration_frozen_before_evaluation"
)
selection_frozen <- if ("configuration_frozen_before_evaluation" %in% names(selection)) {
  value <- selection$configuration_frozen_before_evaluation
  if (is.logical(value)) value else tolower(as.character(value)) == "true"
} else {
  FALSE
}
selection_provenance_ok <-
  nrow(selection) == 54L &&
  !anyDuplicated(paste(selection$dataset_key, selection$method_id)) &&
  all(selection_provenance_fields %in% names(selection)) &&
  all(is.finite(as.numeric(selection$selection_seed))) &&
  all(nzchar(selection$source_package)) &&
  all(nzchar(selection$source_version)) &&
  all(nzchar(selection$information_available_during_selection)) &&
  all(nzchar(selection$departure_justification)) &&
  all(selection_frozen)
add_check(
  "selection_manifest_has_frozen_software_and_seed_provenance",
  selection_provenance_ok,
  paste(nrow(selection), "dataset-method selection records audited")
)

primary_candidates <- read.csv(
  file.path(release_root, "analysis", "config", "primary_competitor_candidate_scores.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
candidate_key <- paste(primary_candidates$dataset_key, primary_candidates$method_id)
candidate_groups <- split(primary_candidates, candidate_key)
candidate_group_ok <- vapply(candidate_groups, function(value) {
  manifest_row <- selection[
    selection$dataset_key == value$dataset_key[[1L]] &
      selection$method_id == value$method_id[[1L]],
    , drop = FALSE
  ]
  nrow(manifest_row) == 1L &&
    nrow(value) == manifest_row$candidate_count &&
    sum(value$selected) == 1L &&
    !any(value$hidden_qc_values_available) &&
    !any(value$biological_labels_available) &&
    !any(value$replicate_identities_available) &&
    !any(value$preservation_endpoints_available)
}, logical(1))
add_check(
  "primary_competitor_rankings_are_complete_and_training_qc_only",
  length(candidate_groups) == 17L && all(candidate_group_ok),
  paste(nrow(primary_candidates), "candidate rows in", length(candidate_groups), "selected configurations")
)

reference_audit_path <- file.path(validation_root, "reference_separation_audit.csv")
if (!file.exists(reference_audit_path)) {
  stop("The release-wide reference-separation audit is missing.", call. = FALSE)
}
reference_audit <- read.csv(reference_audit_path, stringsAsFactors = FALSE, check.names = FALSE)
add_check(
  "training_and_hidden_references_are_separated",
  nrow(reference_audit) > 0L && all(reference_audit$passed),
  paste(nrow(reference_audit), "role, manifest, and recorded-control checks audited")
)

metric_paths <- list.files(
  file.path(result_root, "evaluation"), pattern = "^method_metrics[.]csv$",
  recursive = TRUE, full.names = TRUE
)
metric_tables <- lapply(metric_paths, read.csv, stringsAsFactors = FALSE, check.names = FALSE)
metric_notes <- do.call(rbind, lapply(seq_along(metric_tables), function(index) {
  value <- metric_tables[[index]]
  value <- value[value$metric == "batch_weighted_pc_r2_categorical", , drop = FALSE]
  data.frame(
    path = substring(metric_paths[index], nchar(release_root) + 2L),
    n_rows = nrow(value),
    all_categorical = nrow(value) > 0L && all(value$notes == "Batch is explicitly categorical."),
    stringsAsFactors = FALSE
  )
}))
write.csv(metric_notes, file.path(validation_root, "categorical_batch_r2_audit.csv"), row.names = FALSE)
add_check(
  "batch_r2_is_categorical", all(metric_notes$all_categorical),
  paste(sum(metric_notes$n_rows), "batch R-squared rows audited")
)

metric_crosswalk <- read.csv(
  file.path(release_root, "analysis", "config", "metric_crosswalk.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
partial_metric_paths <- list.files(
  file.path(result_root, "partial_confounding"), pattern = "^method_metrics[.]csv$",
  recursive = TRUE, full.names = TRUE
)
partial_metric_tables <- lapply(
  partial_metric_paths, read.csv, stringsAsFactors = FALSE, check.names = FALSE
)
observed_metric_names <- unique(c(
  unlist(lapply(metric_tables, `[[`, "metric"), use.names = FALSE),
  unlist(lapply(partial_metric_tables, `[[`, "metric"), use.names = FALSE)
))
undocumented_metrics <- setdiff(observed_metric_names, metric_crosswalk$metric)
add_check(
  "metric_crosswalk_covers_saved_metrics", !length(undocumented_metrics),
  if (length(undocumented_metrics)) {
    paste("missing:", paste(undocumented_metrics, collapse = ", "))
  } else {
    paste(length(observed_metric_names), "unique saved metrics documented")
  }
)

all_evaluation_metrics <- do.call(rbind, metric_tables)
icc_rows <- all_evaluation_metrics[
  all_evaluation_metrics$metric == "metabolite_icc_a1_median", , drop = FALSE
]
repeatability_rows <- all_evaluation_metrics[
  all_evaluation_metrics$metric == "feature_repeatability_ratio_median", , drop = FALSE
]
metric_definition_ok <-
  nrow(icc_rows) > 0L && nrow(repeatability_rows) > 0L &&
  all(icc_rows$dataset_key == "mtbls79" & icc_rows$units == "ICC(A,1)") &&
  all(repeatability_rows$dataset_key %in% c("clinical_fiams", "batchcorr_set1") &
        repeatability_rows$units == "repeatability ratio")
add_check(
  "replicate_metrics_have_distinct_definitions", metric_definition_ok,
  paste(nrow(icc_rows), "ICC rows and", nrow(repeatability_rows), "repeatability-ratio rows audited")
)

primary <- read.csv(
  file.path(result_root, "final", "tables", "primary_method_metrics_long.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
winn_parameters <- read.csv(
  file.path(result_root, "final", "tables", "primary_winn_selected_parameters.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
auto_rows <- winn_parameters$method_id %in% c("winn_auto_qc", "winn_auto_batch_qc")
winn_parameter_ok <-
  nrow(winn_parameters) == 18L &&
  !anyDuplicated(paste(winn_parameters$dataset_key, winn_parameters$method_id)) &&
  all(winn_parameters$ljung_box_fitdf == 0L) &&
  all(is.finite(winn_parameters$training_qc_quality_score[auto_rows])) &&
  all(is.na(winn_parameters$training_qc_quality_score[!auto_rows]))
add_check(
  "all_primary_winn_modes_save_selected_parameters",
  winn_parameter_ok,
  paste(nrow(winn_parameters), "dataset-mode parameter records audited")
)
runtime_policy <- read.csv(
  file.path(result_root, "final", "tables", "runtime_measurement_policy.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
reused_fits <- runtime_policy$reuse_decision == "reuse_exact" &
  runtime_policy$method_id != "raw"
runtime_policy_ok <- all(is.finite(runtime_policy$execution_runtime_sec)) &&
  all(is.na(runtime_policy$reported_method_runtime_sec[reused_fits])) &&
  all(grepl("cache retrieval time excluded", runtime_policy$runtime_provenance[reused_fits], fixed = TRUE))
add_check(
  "reused_matrix_loading_time_is_not_reported_as_method_runtime",
  runtime_policy_ok,
  paste(sum(reused_fits), "reused correction fits have reportable runtime set to NA")
)
reuse_audit <- read.csv(
  file.path(result_root, "final", "tables", "competitor_reuse_audit.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
reuse_audit_ok <-
  nrow(reuse_audit) == 36L &&
  all(reuse_audit$final_validation == "passed") &&
  all(reuse_audit$decision_matches_run_manifest) &&
  all(reuse_audit$n_features > 0L & reuse_audit$n_samples > 0L) &&
  all(grepl("^[0-9a-f]{64}$", reuse_audit$input_object_sha256)) &&
  all(grepl("^[0-9a-f]{64}$", reuse_audit$output_object_sha256)) &&
  all(grepl("^[0-9a-f]{64}$", reuse_audit$matrix_file_sha256))
add_check(
  "competitor_reuse_decisions_match_validated_run_artifacts",
  reuse_audit_ok,
  paste(nrow(reuse_audit), "dataset-method reuse records audited")
)
coverage <- primary[
  primary$metric %in% c("retained_features", "retained_samples"),
  c("dataset_key", "method_id", "metric", "value", "denominator"), drop = FALSE
]
coverage$matches_declared_panel <- as.numeric(coverage$value) == as.numeric(coverage$denominator)
write.csv(coverage, file.path(validation_root, "output_coverage_audit.csv"), row.names = FALSE)
add_check(
  "all_primary_outputs_preserve_declared_panels", all(coverage$matches_declared_panel),
  paste(sum(coverage$matches_declared_panel), "of", nrow(coverage), "feature/sample coverage rows")
)

figure_manifest <- read.csv(
  file.path(result_root, "final", "figures", "figure_manifest.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
figure_ok <- vapply(seq_len(nrow(figure_manifest)), function(index) {
  source_ok <- file.exists(file.path(result_root, "final", figure_manifest$source_csv[index]))
  pdf_ok <- !nzchar(figure_manifest$pdf[index]) ||
    file.exists(file.path(result_root, "final", figure_manifest$pdf[index]))
  png_ok <- !nzchar(figure_manifest$png[index]) ||
    file.exists(file.path(result_root, "final", figure_manifest$png[index]))
  source_ok && pdf_ok && png_ok
}, logical(1))
add_check(
  "all_figures_have_source_data", all(figure_ok),
  paste(sum(figure_ok), "of", length(figure_ok), "figure/source entries")
)

checks <- add_check(finish = TRUE)
write.csv(checks, file.path(validation_root, "release_validation.csv"), row.names = FALSE)
utils::capture.output(sessionInfo(), file = file.path(validation_root, "sessionInfo.txt"))
jsonlite::write_json(
  list(
    schema = "endpoint_free_release_validation_v1",
    status = if (all(checks$passed)) "passed" else "failed",
    passed = sum(checks$passed), total = nrow(checks),
    frozen_package_sha256 = expected_sha,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ),
  file.path(validation_root, "validation_manifest.json"), auto_unbox = TRUE, pretty = TRUE
)
if (!all(checks$passed)) {
  stop(
    "Release validation failed: ",
    paste(checks$check[!checks$passed], collapse = ", "), call. = FALSE
  )
}

final_root <- file.path(result_root, "final")
final_files <- list.files(final_root, recursive = TRUE, full.names = TRUE)
manifest_path <- file.path(final_root, "file_manifest.csv")
final_files <- final_files[
  file.info(final_files)$isdir %in% FALSE &
    normalizePath(final_files, mustWork = FALSE) != normalizePath(manifest_path, mustWork = FALSE)
]
final_manifest <- data.frame(
  path = substring(final_files, nchar(final_root) + 2L),
  bytes = as.numeric(file.info(final_files)$size),
  sha256 = vapply(final_files, digest::digest, character(1), file = TRUE, algo = "sha256"),
  stringsAsFactors = FALSE
)
write.csv(final_manifest, manifest_path, row.names = FALSE)
message("Release validation passed (", nrow(checks), " checks).")
