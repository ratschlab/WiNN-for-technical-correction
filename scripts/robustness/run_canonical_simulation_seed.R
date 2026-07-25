#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("scripts/robustness/run_canonical_simulation_seed.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
script_dir <- dirname(script_path)
args <- commandArgs(trailingOnly = TRUE)

value_arg <- function(prefix, default = NULL) {
  found <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (!length(found)) return(default)
  if (length(found) != 1L) stop("Supply exactly one ", prefix, " argument.", call. = FALSE)
  sub(paste0("^", prefix, "="), "", found[[1L]])
}
seed_id <- value_arg("--seed-id")
if (is.null(seed_id) || !grepl("^SIM(0[1-9]|[12][0-9]|30)$", seed_id)) {
  stop("Supply --seed-id=SIM01 through --seed-id=SIM30.", call. = FALSE)
}
force <- "--force" %in% args
dry_run <- "--dry-run" %in% args
winn_only <- "--winn-only" %in% args
if (force && winn_only) {
  stop("Use either --force or --winn-only, not both.", call. = FALSE)
}

project_library <- Sys.getenv("WINN_ROBUSTNESS_R_LIB", unset = "")
if (nzchar(project_library)) {
  if (!dir.exists(project_library)) stop("WINN_ROBUSTNESS_R_LIB does not exist: ", project_library, call. = FALSE)
  .libPaths(c(normalizePath(project_library, mustWork = TRUE), .libPaths()))
}

result_root <- file.path(repo_root, "results", "robustness", "04_simulation_seed_stability")
config_dir <- file.path(result_root, "config")
canonical_run_dir <- file.path(result_root, "runs", seed_id)
preflight_run_dir <- file.path(result_root, "validation", "preflight", seed_id)
# All validation/provenance work is staged outside the canonical run directory.
# This guarantees that a dry run or an exact-reuse check cannot mutate a
# completed, checksummed analysis.
run_dir <- preflight_run_dir
bundle_dir <- file.path(repo_root, "results", "robustness", "03_canonical_simulation", "bundles", seed_id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(run_dir, "logs", "analysis.log")
writeLines(character(), log_path)
log_line <- function(...) {
  value <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " ", paste0(..., collapse = ""))
  cat(value, "\n")
  cat(value, "\n", file = log_path, append = TRUE)
}
run_started <- Sys.time()
log_line("Starting canonical repeated-simulation unit ", seed_id,
         "; force=", force, "; dry_run=", dry_run,
         "; winn_only=", winn_only, ".")

required_packages <- c(
  "devtools", "digest", "jsonlite", "dplyr", "tibble", "limma", "lmtest", "mgcv", "sva",
  "qcrlscR", "statTarget", "TIGERr", "malbacR", "pmartR", "winn"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
devtools::load_all(file.path(repo_root, "winn"), quiet = TRUE)
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source(file.path(repo_root, "scripts", "benchmark_helpers.R"), local = FALSE)
source(file.path(repo_root, "scripts", "public_benchmark_method_helpers.R"), local = FALSE)
source(file.path(repo_root, "scripts", "weighted_pc_r2.R"), local = FALSE)
source(file.path(repo_root, "scripts", "winn_ablation_helpers.R"), local = FALSE)
source(file.path(repo_root, "scripts", "run_order_drift_helpers.R"), local = FALSE)
source(file.path(script_dir, "canonical_cache.R"), local = FALSE)
source(file.path(script_dir, "canonical_seed_method_engine.R"), local = FALSE)
source(file.path(script_dir, "canonical_simulation_metrics.R"), local = FALSE)

sha_file <- canonical_sha256_file
write_csv <- function(value, path) {
  utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
}

parameter_path <- file.path(config_dir, "method_parameters.csv")
parameter_hash_path <- file.path(config_dir, "method_parameters.sha256")
study_config_path <- file.path(config_dir, "study_config.json")
metric_dictionary_path <- file.path(config_dir, "metric_dictionary.csv")
for (path in c(parameter_path, parameter_hash_path, study_config_path, metric_dictionary_path)) {
  if (!file.exists(path)) stop("Missing prepared study configuration: ", path, call. = FALSE)
}
expected_parameter_hash <- trimws(readLines(parameter_hash_path, warn = FALSE)[[1L]])
observed_parameter_hash <- sha_file(parameter_path)
if (!identical(expected_parameter_hash, observed_parameter_hash)) {
  stop("Prepared method parameters failed SHA-256 validation.", call. = FALSE)
}
method_parameters <- read.csv(parameter_path, stringsAsFactors = FALSE, check.names = FALSE)
metric_dictionary <- read.csv(metric_dictionary_path, stringsAsFactors = FALSE, check.names = FALSE)
study_config <- jsonlite::read_json(study_config_path, simplifyVector = TRUE)
if (!seed_id %in% study_config$seed_ids) stop(seed_id, " is not in the prepared study config.", call. = FALSE)

required_bundle_files <- c(
  "raw_intensity.rds", "clean_ground_truth.rds", "sample_metadata.csv",
  "feature_truth.csv", "feature_plate_truth.csv",
  "hidden_reference_assignment.csv", "seed_ledger.csv", "bundle_files.csv",
  "bundle_provenance.json"
)
missing_bundle <- required_bundle_files[!file.exists(file.path(bundle_dir, required_bundle_files))]
if (length(missing_bundle)) stop(seed_id, " bundle is incomplete: ", paste(missing_bundle, collapse = ", "), call. = FALSE)
bundle_files <- read.csv(file.path(bundle_dir, "bundle_files.csv"), stringsAsFactors = FALSE, check.names = FALSE)
bundle_paths <- file.path(bundle_dir, bundle_files$relative_path)
bundle_hash_ok <- all(file.exists(bundle_paths)) &&
  all(as.numeric(file.info(bundle_paths)$size) == bundle_files$bytes) &&
  identical(
    unname(vapply(bundle_paths, sha_file, character(1))),
    unname(bundle_files$sha256)
  )
if (!bundle_hash_ok) stop(seed_id, " bundle file manifest validation failed.", call. = FALSE)
bundle_provenance <- jsonlite::read_json(file.path(bundle_dir, "bundle_provenance.json"), simplifyVector = TRUE)

x <- readRDS(file.path(bundle_dir, "raw_intensity.rds"))
truth <- readRDS(file.path(bundle_dir, "clean_ground_truth.rds"))
meta <- read.csv(file.path(bundle_dir, "sample_metadata.csv"), stringsAsFactors = FALSE, check.names = FALSE)
feature_truth <- read.csv(file.path(bundle_dir, "feature_truth.csv"), stringsAsFactors = FALSE, check.names = FALSE)
hidden_assignment <- read.csv(file.path(bundle_dir, "hidden_reference_assignment.csv"), stringsAsFactors = FALSE, check.names = FALSE)
seed_ledger <- read.csv(file.path(bundle_dir, "seed_ledger.csv"), stringsAsFactors = FALSE, check.names = FALSE)
required_method_rng_components <- unique(method_parameters$rng_component)

input_checks <- data.frame(
  check = c(
    "bundle_manifest_hashes", "raw_dimensions", "truth_dimensions",
    "raw_truth_dimnames", "metadata_sample_order", "numeric_finite_inputs",
    "five_plates", "fifty_total_qcs", "ten_hidden_qcs",
    "forty_training_qcs", "two_hidden_per_plate", "eight_training_per_plate",
    "hidden_assignment_matches_bundle", "unique_ids_and_orders",
    "feature_truth_alignment", "method_rng_seed_ledger",
    "analysis_schema", "parameter_freeze_hash"
  ),
  passed = FALSE,
  detail = "",
  stringsAsFactors = FALSE
)
hidden_ids <- meta$sample_id[meta$is_hidden_reference]
training_ids <- meta$sample_id[meta$is_qc & !meta$is_hidden_reference]
split_counts <- aggregate(
  cbind(hidden = meta$sample_id %in% hidden_ids, training = meta$sample_id %in% training_ids),
  list(batch = meta$batch), sum
)
hidden_positions <- meta$within_batch_order[meta$sample_id %in% hidden_ids]
input_checks$passed <- c(
  bundle_hash_ok,
  identical(dim(x), c(1000L, 500L)),
  identical(dim(truth), c(1000L, 500L)),
  identical(dimnames(x), dimnames(truth)),
  identical(colnames(x), meta$sample_id),
  is.numeric(x) && is.numeric(truth) && all(is.finite(x)) && all(is.finite(truth)),
  length(unique(meta$batch)) == 5L,
  sum(meta$is_qc) == 50L,
  length(hidden_ids) == 10L,
  length(training_ids) == 40L,
  nrow(split_counts) == 5L && all(split_counts$hidden == 2L),
  nrow(split_counts) == 5L && all(split_counts$training == 8L),
  identical(
    meta$is_hidden_reference,
    meta$sample_id %in% hidden_assignment$sample_id[hidden_assignment$is_hidden_reference]
  ) && all(meta$is_qc[meta$is_hidden_reference]),
  !anyDuplicated(meta$sample_id) && identical(meta$run_order, seq_len(ncol(x))) &&
    all(vapply(split(meta$within_batch_order, meta$batch), function(value) {
      identical(value, seq_len(100L))
    }, logical(1))),
  identical(rownames(x), feature_truth$metabolite) && !anyDuplicated(feature_truth$metabolite),
  all(seed_ledger$seed_id == seed_id) && !anyDuplicated(seed_ledger$component) &&
    all(required_method_rng_components %in% seed_ledger$component) &&
    all(is.finite(seed_ledger$effective_seed[match(required_method_rng_components, seed_ledger$component)])),
  nrow(method_parameters) == 18L && !anyDuplicated(method_parameters$method) &&
    nrow(metric_dictionary) == 28L && !anyDuplicated(metric_dictionary$metric),
  identical(observed_parameter_hash, expected_parameter_hash)
)
input_checks$detail <- c(
  paste(nrow(bundle_files), "bundle artifacts"), paste(dim(x), collapse = "x"),
  paste(dim(truth), collapse = "x"), as.character(identical(dimnames(x), dimnames(truth))),
  as.character(identical(colnames(x), meta$sample_id)), "numeric finite matrices",
  paste(unique(meta$batch), collapse = ";"), as.character(sum(meta$is_qc)),
  as.character(length(hidden_ids)), as.character(length(training_ids)),
  paste(split_counts$hidden, collapse = ";"), paste(split_counts$training, collapse = ";"),
  paste0("positions=", paste(hidden_positions, collapse = ";")),
  "unique sample IDs; global order 1:500; within-plate order 1:100",
  paste(nrow(feature_truth), "aligned feature-truth rows"),
  paste(required_method_rng_components, collapse = ";"),
  paste0(nrow(method_parameters), " methods; ", nrow(metric_dictionary), " metrics"),
  observed_parameter_hash
)
write_csv(input_checks, file.path(run_dir, "input_validation.csv"))
if (!all(input_checks$passed)) {
  stop("Input validation failed: ", paste(input_checks$check[!input_checks$passed], collapse = ", "), call. = FALSE)
}

meta_blind <- meta
meta_blind$class <- ifelse(meta_blind$sample_id %in% training_ids, "QC", "Sample")
meta_blind$is_qc <- meta_blind$sample_id %in% training_ids
meta_blind$sample_type <- ifelse(meta_blind$is_qc, "control", "sample")
reference_assignment <- data.frame(
  seed_id = seed_id,
  sample_id = meta$sample_id,
  batch = meta$batch,
  run_order = meta$run_order,
  within_batch_order = meta$within_batch_order,
  evaluation_role = ifelse(
    meta$sample_id %in% hidden_ids, "hidden_test",
    ifelse(meta$sample_id %in% training_ids, "training_control", "study")
  ),
  method_visible_class = meta_blind$class,
  supplied_as_control = meta$sample_id %in% training_ids,
  stringsAsFactors = FALSE
)
write_csv(reference_assignment, file.path(run_dir, "reference_assignment.csv"))
write_csv(split_counts, file.path(run_dir, "reference_split_counts.csv"))

analysis_units <- method_parameters$method
qc_methods <- c("QC-RLSC", "QC-RFSC", "TIGER", "SERRF", "WINN auto (QC)", "WINN auto-batch (QC)")
exposure <- data.frame(
  seed_id = seed_id,
  method = analysis_units,
  n_training_controls_supplied = ifelse(analysis_units %in% qc_methods, length(training_ids), 0L),
  n_hidden_controls_supplied = 0L,
  hidden_labels_visible = FALSE,
  hidden_used_for_tuning = FALSE,
  one_output_for_all_endpoints = TRUE,
  hidden_ids_hash = canonical_sha256_object(hidden_ids),
  stringsAsFactors = FALSE
)
write_csv(exposure, file.path(run_dir, "method_exposure_protocol.csv"))

source_files <- c(
  runner = script_path,
  method_engine = file.path(script_dir, "canonical_seed_method_engine.R"),
  metric_engine = file.path(script_dir, "canonical_simulation_metrics.R"),
  cache = file.path(script_dir, "canonical_cache.R"),
  benchmark_helpers = file.path(repo_root, "scripts", "benchmark_helpers.R"),
  public_helpers = file.path(repo_root, "scripts", "public_benchmark_method_helpers.R"),
  weighted_r2 = file.path(repo_root, "scripts", "weighted_pc_r2.R"),
  ablation_helpers = file.path(repo_root, "scripts", "winn_ablation_helpers.R"),
  drift_helpers = file.path(repo_root, "scripts", "run_order_drift_helpers.R"),
  winn_core = file.path(repo_root, "winn", "R", "winn.R"),
  winn_ablation = file.path(repo_root, "winn", "R", "winn-ablation.R")
)
source_hashes <- vapply(source_files, sha_file, character(1))
code_bundle_sha256 <- canonical_sha256_object(as.list(source_hashes))
winn_git <- function(arguments) {
  suppressWarnings(tryCatch(
    system2(
      "git", c("-C", shQuote(file.path(repo_root, "winn")), arguments),
      stdout = TRUE, stderr = FALSE
    ),
    error = function(e) character()
  ))
}
winn_commit_output <- winn_git(c("rev-parse", "HEAD"))
winn_commit <- if (length(winn_commit_output)) winn_commit_output[[1L]] else NA_character_
winn_status <- winn_git(c("status", "--short"))
if (is.na(winn_commit) || !nzchar(winn_commit)) {
  stop("Unable to record the local WiNN Git commit; refusing an unaudited run.", call. = FALSE)
}
global_context <- list(
  schema = "canonical_seed_stability_run_v2",
  seed_id = seed_id,
  bundle_config_sha256 = bundle_provenance$config_sha256,
  bundle_file_sha256 = stats::setNames(bundle_files$sha256, bundle_files$relative_path),
  frozen_parameters_sha256 = observed_parameter_hash,
  source_sha256 = as.list(source_hashes),
  code_bundle_sha256 = code_bundle_sha256,
  winn_commit = winn_commit,
  winn_status_sha256 = canonical_sha256_object(winn_status),
  training_ids_sha256 = canonical_sha256_object(training_ids),
  hidden_ids_sha256 = canonical_sha256_object(hidden_ids),
  label_blinding = "only training controls labeled QC; hidden references labeled ordinary Sample",
  refresh_scope = if (winn_only) "winn_only_with_frozen_competitors" else "all_methods"
)
analysis_context_sha256 <- canonical_sha256_object(global_context)

input_checksums <- data.frame(
  role = c(names(source_files), "method_parameters", "bundle_provenance", "raw_intensity", "clean_ground_truth", "sample_metadata", "feature_truth", "hidden_assignment", "seed_ledger"),
  path = c(
    unname(source_files), parameter_path, file.path(bundle_dir, "bundle_provenance.json"),
    file.path(bundle_dir, "raw_intensity.rds"), file.path(bundle_dir, "clean_ground_truth.rds"),
    file.path(bundle_dir, "sample_metadata.csv"), file.path(bundle_dir, "feature_truth.csv"),
    file.path(bundle_dir, "hidden_reference_assignment.csv"), file.path(bundle_dir, "seed_ledger.csv")
  ),
  stringsAsFactors = FALSE
)
input_checksums$bytes <- file.info(input_checksums$path)$size
input_checksums$sha256 <- vapply(input_checksums$path, sha_file, character(1))
input_checksums$path <- vapply(input_checksums$path, function(path) {
  if (startsWith(path, paste0(repo_root, .Platform$file.sep))) substring(path, nchar(repo_root) + 2L) else path
}, character(1))
write_csv(input_checksums, file.path(run_dir, "input_checksums.csv"))

completion_path <- file.path(canonical_run_dir, "analysis_complete.json")
artifact_path <- file.path(canonical_run_dir, "artifact_checksums.csv")
validate_completed_run <- function() {
  if (!file.exists(completion_path) || !file.exists(artifact_path)) return(FALSE)
  completion <- tryCatch(jsonlite::read_json(completion_path, simplifyVector = TRUE), error = function(e) NULL)
  completion_valid <- !is.null(completion) &&
    identical(completion$schema, "canonical_seed_complete_v1") &&
    identical(completion$seed_id, seed_id) &&
    identical(completion$analysis_context_sha256, analysis_context_sha256) &&
    isTRUE(completion$invariants_passed) && isTRUE(completion$metrics_complete) &&
    as.integer(completion$method_count) == 18L &&
    as.integer(completion$metric_rows) == 18L * nrow(metric_dictionary)
  if (!completion_valid) return(FALSE)
  artifacts <- tryCatch(read.csv(artifact_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(artifacts) || !all(c("relative_path", "bytes", "sha256") %in% names(artifacts))) return(FALSE)
  required_completed_artifacts <- c(
    "analysis_complete.json", "run_manifest.csv", "method_status.csv",
    "errors.csv", "runtime.csv", "method_metrics.csv", "validation_checks.csv",
    "method_matrix_checks.csv", "matrix_checksums.csv", "metric_equivalence_checks.csv",
    "method_exposure_protocol.csv", "selected_parameters.csv", "sessionInfo.txt"
  )
  if (!all(required_completed_artifacts %in% artifacts$relative_path)) return(FALSE)
  paths <- file.path(canonical_run_dir, artifacts$relative_path)
  artifact_valid <- all(file.exists(paths)) && all(file.info(paths)$size == artifacts$bytes) &&
    identical(unname(vapply(paths, sha_file, character(1))), unname(artifacts$sha256))
  if (!artifact_valid) return(FALSE)
  manifest <- tryCatch(
    read.csv(file.path(canonical_run_dir, "run_manifest.csv"), stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  !is.null(manifest) && nrow(manifest) == 1L &&
    identical(as.character(manifest$seed_id), seed_id) &&
    identical(as.character(manifest$analysis_context_sha256), analysis_context_sha256) &&
    identical(as.character(manifest$status), "completed") &&
    as.integer(manifest$methods_completed) == 18L &&
    as.integer(manifest$metrics_rows) == 18L * nrow(metric_dictionary) &&
    isTRUE(as.logical(manifest$invariants_passed)) &&
    isTRUE(as.logical(manifest$hidden_exclusion_passed))
}

if (dry_run) {
  dry <- data.frame(
    seed_id = seed_id,
    bundle_valid = all(input_checks$passed),
    package_dependencies_available = !length(missing_packages),
    parameter_freeze_valid = identical(observed_parameter_hash, expected_parameter_hash),
    hidden_reference_protocol_valid = length(hidden_ids) == 10L && length(training_ids) == 40L && !any(meta_blind$sample_id[meta_blind$class == "QC"] %in% hidden_ids),
    source_files_hashable = all(file.exists(source_files)),
    existing_completed_run_reusable = validate_completed_run(),
    analysis_context_sha256 = analysis_context_sha256,
    stringsAsFactors = FALSE
  )
  validation_dir <- file.path(result_root, "validation")
  dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(dry, file.path(validation_dir, paste0("dry_run_", seed_id, ".csv")))
  message("Dry-run validation passed for ", seed_id, "; no correction method was executed.")
  quit(save = "no", status = 0L)
}

completed_run_reusable <- validate_completed_run()
if (!force && !winn_only && completed_run_reusable) {
  message(seed_id, " already has a complete content-matched run; reusing it without recomputation.")
  quit(save = "no", status = 0L)
}

canonical_has_files <- dir.exists(canonical_run_dir) &&
  length(list.files(canonical_run_dir, all.files = TRUE, no.. = TRUE)) > 0L
canonical_has_completion_artifact <- file.exists(completion_path) || file.exists(artifact_path)
if (canonical_has_files && !force && !winn_only && canonical_has_completion_artifact) {
  stop(
    seed_id, " has a completed-looking run that does not match the current content context. ",
    "Refusing to overwrite it; use --force to archive that run and recompute.",
    call. = FALSE
  )
}
if (canonical_has_files && force) {
  archive_parent <- file.path(result_root, "archive")
  dir.create(archive_parent, recursive = TRUE, showWarnings = FALSE)
  archive_dir <- file.path(
    archive_parent,
    paste0(seed_id, "_", format(Sys.time(), "%Y%m%dT%H%M%S"), "_", substr(analysis_context_sha256, 1L, 12L))
  )
  if (file.exists(archive_dir) || !file.rename(canonical_run_dir, archive_dir)) {
    stop("Could not archive the existing canonical run before forced recomputation.", call. = FALSE)
  }
}

if (isTRUE(winn_only)) {
  required_frozen_cache <- file.path(
    canonical_run_dir, "cache",
    paste0(c("Raw", "ComBat", "QC_RLSC", "QC_RFSC", "TIGER", "SERRF"), ".rds")
  )
  if (!all(file.exists(required_frozen_cache))) {
    stop(
      "WiNN-only refresh requires all frozen non-WiNN caches in the existing run: ",
      paste(required_frozen_cache[!file.exists(required_frozen_cache)], collapse = ", "),
      call. = FALSE
    )
  }
}

run_dir <- canonical_run_dir
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
resume_stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
old_log <- file.path(run_dir, "logs", "analysis.log")
if (file.exists(old_log)) {
  file.rename(old_log, file.path(run_dir, "logs", paste0("prior_analysis_", resume_stamp, ".log")))
}
old_failure <- file.path(run_dir, "run_failure.json")
if (file.exists(old_failure)) {
  file.rename(old_failure, file.path(run_dir, "logs", paste0("prior_run_failure_", resume_stamp, ".json")))
}
log_path <- file.path(run_dir, "logs", "analysis.log")
writeLines(character(), log_path)
log_line("Beginning correction and metric execution for ", seed_id, ".")
for (audit_file in c(
  "input_validation.csv", "reference_assignment.csv", "reference_split_counts.csv",
  "method_exposure_protocol.csv", "input_checksums.csv"
)) {
  source_path <- file.path(preflight_run_dir, audit_file)
  if (!file.copy(source_path, file.path(run_dir, audit_file), overwrite = TRUE)) {
    stop("Could not stage preflight audit artifact: ", audit_file, call. = FALSE)
  }
}

run_canonical_cache_negative_tests(file.path(run_dir, "cache_key_negative_tests.csv"))
execution <- run_canonical_seed_methods(
  x = x, meta_blind = meta_blind, training_ids = training_ids,
  seed_ledger = seed_ledger, method_parameters = method_parameters,
  run_dir = run_dir, global_context = global_context,
  force = force, winn_only = winn_only, log_line = log_line
)

method_order <- method_parameters$method
all_methods_completed <- nrow(execution$status) == length(method_order) &&
  identical(execution$status$method, method_order) &&
  all(grepl("^completed", execution$status$status))
matrix_checks <- dplyr::bind_rows(lapply(names(execution$matrices), function(method) {
  matrix <- execution$matrices[[method]]
  data.frame(
    seed_id = seed_id, method = method,
    numeric_matrix = is.matrix(matrix) && is.numeric(matrix),
    exact_dimensions = identical(dim(matrix), dim(x)),
    exact_feature_order = identical(rownames(matrix), rownames(x)),
    exact_sample_order = identical(colnames(matrix), colnames(x)),
    all_finite = all(is.finite(matrix)),
    log1p_valid = all(matrix > -1),
    stringsAsFactors = FALSE
  )
}))
write_csv(matrix_checks, file.path(run_dir, "method_matrix_checks.csv"))

max_difference <- function(a, b) max(abs(a - b), na.rm = TRUE)
has <- function(values) all(values %in% names(execution$matrices))
invariants <- data.frame(
  seed_id = seed_id,
  check = c(
    "all_18_analysis_units_attempted", "all_18_analysis_units_completed",
    "raw_equals_c0", "fixed_equals_c4", "c4_equals_gss",
    "c4_equals_independent_current_winn_fixed",
    "all_matrices_exactly_aligned", "same_hidden_references_all_methods",
    "hidden_references_absent_from_controls_and_tuning",
    "single_output_per_method_for_all_metrics", "canonical_bundle_identity_recorded"
  ),
  observed = c(
    nrow(execution$status), sum(grepl("^completed", execution$status$status)),
    if (has(c("Raw", "C0_RAW"))) max_difference(execution$matrices$Raw, execution$matrices$C0_RAW) else Inf,
    if (has(c("WINN default (no QC)", "C4_FULL_FIXED"))) max_difference(execution$matrices[["WINN default (no QC)"]], execution$matrices$C4_FULL_FIXED) else Inf,
    if (has(c("C4_FULL_FIXED", "G_SS"))) max_difference(execution$matrices$C4_FULL_FIXED, execution$matrices$G_SS) else Inf,
    if (has(c("WINN default (no QC)", "C4_FULL_FIXED"))) max_difference(execution$matrices$C4_FULL_FIXED, execution$matrices[["WINN default (no QC)"]]) else Inf,
    nrow(matrix_checks), length(unique(exposure$hidden_ids_hash)),
    sum(exposure$n_hidden_controls_supplied), all(exposure$one_output_for_all_endpoints),
    bundle_provenance$config_sha256
  ),
  tolerance = c(NA, NA, rep(1e-8, 4L), rep(NA, 5L)),
  stringsAsFactors = FALSE
)
invariants$passed <- c(
  nrow(execution$status) == 18L && all(execution$status$attempted),
  all_methods_completed,
  is.finite(as.numeric(invariants$observed[3])) && as.numeric(invariants$observed[3]) <= 1e-8,
  is.finite(as.numeric(invariants$observed[4])) && as.numeric(invariants$observed[4]) <= 1e-8,
  is.finite(as.numeric(invariants$observed[5])) && as.numeric(invariants$observed[5]) <= 1e-8,
  is.finite(as.numeric(invariants$observed[6])) && as.numeric(invariants$observed[6]) <= 1e-8,
  nrow(matrix_checks) == 18L && all(as.matrix(matrix_checks[, -(1:2), drop = FALSE])),
  length(unique(exposure$hidden_ids_hash)) == 1L,
  all(exposure$n_hidden_controls_supplied == 0L & !exposure$hidden_labels_visible & !exposure$hidden_used_for_tuning),
  all(exposure$one_output_for_all_endpoints),
  nzchar(bundle_provenance$config_sha256)
)
write_csv(invariants, file.path(run_dir, "validation_checks.csv"))

if (!all(invariants$passed)) {
  failure <- list(
    schema = "canonical_seed_failure_v1", seed_id = seed_id,
    analysis_context_sha256 = analysis_context_sha256,
    failed_checks = invariants$check[!invariants$passed],
    method_failures = execution$status$method[!grepl("^completed", execution$status$status)],
    stopped_before_metrics = TRUE,
    recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  jsonlite::write_json(failure, file.path(run_dir, "run_failure.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
  stop(seed_id, " failed required invariants; full metric evaluation was not run.", call. = FALSE)
}

metrics <- compute_canonical_seed_metrics(
  seed_id = seed_id, method_matrices = execution$matrices,
  raw = x, truth = truth, meta = meta, hidden_ids = hidden_ids,
  feature_truth = feature_truth, runtime = execution$runtime,
  ablation_diagnostics = execution$ablation_diagnostics,
  ablation_configurations = execution$ablation_configurations,
  metric_dictionary = metric_dictionary, run_dir = run_dir,
  log_line = log_line
)
write_csv(metrics$detailed_row_counts, file.path(run_dir, "metric_detail_row_counts.csv"))

metric_values <- function(method) {
  rows <- metrics$method_metrics[metrics$method_metrics$method == method, , drop = FALSE]
  rows$value[match(metric_dictionary$metric, rows$metric)]
}
compare_metric_alias <- function(check, method_a, method_b) {
  value_a <- metric_values(method_a)
  value_b <- metric_values(method_b)
  evaluated <- metric_dictionary$metric != "runtime_sec"
  same_na <- identical(is.na(value_a[evaluated]), is.na(value_b[evaluated]))
  finite <- evaluated & is.finite(value_a) & is.finite(value_b)
  max_abs_difference <- if (any(finite)) max(abs(value_a[finite] - value_b[finite])) else 0
  data.frame(
    seed_id = seed_id, check = check, method_a = method_a, method_b = method_b,
    metrics_compared = sum(evaluated), same_na_pattern = same_na,
    max_abs_difference = max_abs_difference, tolerance = 1e-12,
    passed = same_na && is.finite(max_abs_difference) && max_abs_difference <= 1e-12,
    excluded_metric = "runtime_sec", stringsAsFactors = FALSE
  )
}
metric_equivalences <- dplyr::bind_rows(
  compare_metric_alias("raw_equals_c0_metrics", "Raw", "C0_RAW"),
  compare_metric_alias("fixed_equals_c4_metrics", "WINN default (no QC)", "C4_FULL_FIXED"),
  compare_metric_alias("c4_equals_gss_metrics", "C4_FULL_FIXED", "G_SS")
)
write_csv(metric_equivalences, file.path(run_dir, "metric_equivalence_checks.csv"))
if (!all(metric_equivalences$passed)) {
  jsonlite::write_json(list(
    schema = "canonical_seed_failure_v1", seed_id = seed_id,
    analysis_context_sha256 = analysis_context_sha256,
    failed_checks = metric_equivalences$check[!metric_equivalences$passed],
    method_failures = character(), stopped_after_metrics = TRUE,
    recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ), file.path(run_dir, "run_failure.json"), auto_unbox = TRUE, pretty = TRUE)
  stop(seed_id, " failed metric-alias equivalence checks; completion was not declared.", call. = FALSE)
}

extract_message <- function(text, pattern, default = NA_character_) {
  match <- regexec(pattern, text, perl = TRUE)
  found <- regmatches(text, match)[[1L]]
  if (length(found) >= 2L) trimws(found[[2L]]) else default
}
selected_parameters <- lapply(seq_len(nrow(method_parameters)), function(index) {
  parameter <- method_parameters[index, , drop = FALSE]
  method <- parameter$method[[1L]]
  text <- execution$warnings_messages$messages[match(method, execution$warnings_messages$method)]
  if (!length(text) || is.na(text)) text <- ""
  selected_json <- parameter$parameters_json[[1L]]
  source <- "frozen before SIM01 pilot and repeated-seed execution"
  if (method %in% c("WINN auto (QC)", "WINN auto-batch (QC)")) {
    selected_json <- canonical_seed_plain_json(list(
      parameters = "auto",
      batch_detection = extract_message(text, "Batch detection: ([^\\n]+)"),
      pelt_penalty = suppressWarnings(as.numeric(extract_message(text, "fkPELT penalty: ([0-9.eE+-]+)"))),
      autocorrelation_test = extract_message(text, "Autocorrelation test: ([^ (\\n]+)"),
      autocorrelation_fdr = suppressWarnings(as.numeric(extract_message(text, "Autocorrelation test: [^\\n]+\\(FDR: ([0-9.]+)\\)"))),
      lag = extract_message(text, "Ljung-Box lag: ([^\\n]+)"),
      spline_method = extract_message(text, "Spline method: ([^\\n]+)"),
      batch_fdr = suppressWarnings(as.numeric(extract_message(text, "Batch correction FDR: ([0-9.]+)"))),
      normalization = extract_message(text, "Normalization: ([^\\n]+)"),
      scale_by_batch = identical(toupper(extract_message(text, "Scale by batch: ([^\\n]+)")), "TRUE"),
      training_qc_cv = suppressWarnings(as.numeric(extract_message(text, "Quality metrics - CV: ([0-9.eE+-]+)"))),
      training_qc_correlation = suppressWarnings(as.numeric(extract_message(text, "Correlation: ([0-9.eE+-]+)")))
    ), auto_unbox = TRUE, null = "null", na = "null")
    source <- "selected automatically from this seed's 40 training controls"
  }
  data.frame(
    seed_id = seed_id, method = method,
    selected_parameters_json = as.character(selected_json),
    selection_source = source,
    hidden_references_used = FALSE,
    stringsAsFactors = FALSE
  )
})
selected_parameters <- dplyr::bind_rows(selected_parameters)
if (nrow(selected_parameters) != 18L || !is.character(selected_parameters$selected_parameters_json)) {
  stop("Selected-parameter audit table failed its 18-row plain-character JSON schema check.", call. = FALSE)
}
write_csv(selected_parameters, file.path(run_dir, "selected_parameters.csv"))

software <- list(
  R = R.version.string,
  platform = R.version$platform,
  library_paths = .libPaths(),
  project_library = if (nzchar(project_library)) normalizePath(project_library, mustWork = TRUE) else NULL,
  packages = stats::setNames(lapply(required_packages, function(package) as.character(utils::packageVersion(package))), required_packages),
  source_sha256 = as.list(source_hashes),
  local_winn_commit = winn_commit,
  local_winn_dirty = length(winn_status) > 0L,
  local_winn_status = winn_status
)
jsonlite::write_json(software, file.path(run_dir, "software_state.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
utils::capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))

run_completed <- Sys.time()
manifest_row <- data.frame(
  seed_id = seed_id,
  display_seed = identical(seed_id, study_config$display_seed),
  status = "completed",
  run_started = format(run_started, "%Y-%m-%dT%H:%M:%S%z"),
  run_completed = format(run_completed, "%Y-%m-%dT%H:%M:%S%z"),
  elapsed_wall_sec = as.numeric(difftime(run_completed, run_started, units = "secs")),
  bundle_config_sha256 = bundle_provenance$config_sha256,
  bundle_raw_sha256 = bundle_files$sha256[match("raw_intensity.rds", bundle_files$relative_path)],
  bundle_truth_sha256 = bundle_files$sha256[match("clean_ground_truth.rds", bundle_files$relative_path)],
  hidden_assignment_sha256 = sha_file(file.path(bundle_dir, "hidden_reference_assignment.csv")),
  frozen_parameters_sha256 = observed_parameter_hash,
  code_bundle_sha256 = code_bundle_sha256,
  winn_commit = winn_commit,
  winn_status_sha256 = canonical_sha256_object(winn_status),
  analysis_context_sha256 = analysis_context_sha256,
  n_features = nrow(x), n_samples = ncol(x), n_batches = length(unique(meta$batch)),
  n_qc_total = sum(meta$is_qc), n_hidden = length(hidden_ids), n_training = length(training_ids),
  methods_attempted = nrow(execution$status), methods_completed = sum(grepl("^completed", execution$status$status)),
  metrics_rows = nrow(metrics$method_metrics), invariants_passed = all(invariants$passed),
  metric_equivalences_passed = all(metric_equivalences$passed),
  hidden_exclusion_passed = all(exposure$n_hidden_controls_supplied == 0L & !exposure$hidden_used_for_tuning),
  refresh_scope = if (winn_only) "winn_only_with_frozen_competitors" else "all_methods",
  stringsAsFactors = FALSE
)
write_csv(manifest_row, file.path(run_dir, "run_manifest.csv"))
jsonlite::write_json(as.list(manifest_row[1, ]), file.path(run_dir, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")

completion <- list(
  schema = "canonical_seed_complete_v1", seed_id = seed_id,
  analysis_context_sha256 = analysis_context_sha256,
  bundle_config_sha256 = bundle_provenance$config_sha256,
  invariants_passed = TRUE, metrics_complete = TRUE,
  method_count = 18L, metric_rows = nrow(metrics$method_metrics),
  refresh_scope = if (winn_only) "winn_only_with_frozen_competitors" else "all_methods",
  completed_at = format(run_completed, "%Y-%m-%dT%H:%M:%S%z")
)
jsonlite::write_json(completion, completion_path, auto_unbox = TRUE, pretty = TRUE)

log_line(
  seed_id, " completed: 18/18 analysis units, ", nrow(metrics$method_metrics),
  " metric rows, all invariants passed; wall_sec=",
  sprintf("%.1f", manifest_row$elapsed_wall_sec), "."
)

# Write this last and do not mutate any listed artifact afterwards.
all_files <- list.files(run_dir, recursive = TRUE, full.names = TRUE)
all_files <- all_files[!grepl("(^|/)cache(/|$)", all_files)]
all_files <- setdiff(all_files, artifact_path)
all_files <- all_files[file.info(all_files)$isdir %in% FALSE]
artifacts <- data.frame(
  relative_path = substring(all_files, nchar(run_dir) + 2L),
  bytes = file.info(all_files)$size,
  sha256 = vapply(all_files, sha_file, character(1)),
  stringsAsFactors = FALSE
)
write_csv(artifacts, artifact_path)
message(seed_id, " canonical repeated-simulation run complete; ", nrow(artifacts), " artifacts hashed.")
