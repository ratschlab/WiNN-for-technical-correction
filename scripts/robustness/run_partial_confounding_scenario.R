#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "scripts/robustness/run_partial_confounding_scenario.R",
    mustWork = TRUE
  )
}
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)

value_arg <- function(prefix, default = NULL) {
  found <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (!length(found)) return(default)
  if (length(found) != 1L) stop("Supply exactly one ", prefix, " argument.",
                                call. = FALSE)
  sub(paste0("^", prefix, "="), "", found[[1L]])
}
scenario_id <- value_arg("--scenario-id")
analysis_mode <- value_arg("--analysis-mode", "full")
if (!analysis_mode %in% c("full", "winn_only")) {
  stop("--analysis-mode must be full or winn_only.", call. = FALSE)
}
winn_only <- identical(analysis_mode, "winn_only")
scenario_pattern <- paste0(
  "^CONF(0[1-9]|1[0-9]|20)_",
  "(none_0p00|(batch|run_order|joint)_(0p25|0p50|0p75|0p90|1p00))$"
)
if (is.null(scenario_id) || !grepl(scenario_pattern, scenario_id)) {
  stop("Supply one valid --scenario-id=CONFxx_axis_strength.", call. = FALSE)
}
force <- "--force" %in% args
dry_run <- "--dry-run" %in% args

project_library <- Sys.getenv("WINN_ROBUSTNESS_R_LIB", unset = "")
if (nzchar(project_library)) {
  if (!dir.exists(project_library)) {
    stop("WINN_ROBUSTNESS_R_LIB does not exist: ", project_library,
         call. = FALSE)
  }
  .libPaths(c(normalizePath(project_library, mustWork = TRUE), .libPaths()))
}

result_root <- file.path(
  repo_root, "results", "robustness", "06_partial_confounding"
)
full_config_root <- file.path(result_root, "full_grid")
full_root <- if (winn_only) {
  file.path(result_root, "winn_only_grid")
} else {
  full_config_root
}
config_dir <- file.path(full_config_root, "config")
canonical_run_dir <- file.path(full_root, "runs", scenario_id)
preflight_run_dir <- file.path(
  full_root, "validation", "preflight", scenario_id
)
preflight_has_files <- dir.exists(preflight_run_dir) &&
  length(list.files(preflight_run_dir, all.files = TRUE, no.. = TRUE)) > 0L
canonical_completion_candidate <-
  file.exists(file.path(canonical_run_dir, "analysis_complete.json")) &&
  file.exists(file.path(canonical_run_dir, "artifact_checksums.csv"))
reuse_validation_scratch <- ""
run_dir <- preflight_run_dir
if (preflight_has_files && !force && canonical_completion_candidate) {
  reuse_validation_scratch <- tempfile(
    pattern = paste0("winn_partial_reuse_", scenario_id, "_"),
    tmpdir = tempdir()
  )
  if (!dir.create(reuse_validation_scratch, recursive = TRUE,
                  showWarnings = FALSE)) {
    stop("Could not create completed-run reuse validation scratch space.",
         call. = FALSE)
  }
  run_dir <- reuse_validation_scratch
} else if (preflight_has_files && !force) {
  stop(
    scenario_id,
    " has prior-version or prior preflight output; use --force to archive it before v3 execution.",
    call. = FALSE
  )
}
if (preflight_has_files && force) {
  preflight_archive_parent <- file.path(full_root, "archive", "preflight")
  dir.create(preflight_archive_parent, recursive = TRUE, showWarnings = FALSE)
  preflight_archive <- file.path(
    preflight_archive_parent,
    paste0(scenario_id, "_", format(Sys.time(), "%Y%m%dT%H%M%S"),
           "_prior_version")
  )
  if (file.exists(preflight_archive) ||
      !file.rename(preflight_run_dir, preflight_archive)) {
    stop("Could not archive prior scenario preflight output.", call. = FALSE)
  }
}
cleanup_reuse_validation_scratch <- function() {
  if (nzchar(reuse_validation_scratch) &&
      dir.exists(reuse_validation_scratch)) {
    unlink(reuse_validation_scratch, recursive = TRUE, force = TRUE)
  }
}
dir.create(file.path(run_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(run_dir, "logs", "analysis.log")
writeLines(character(), log_path)
log_line <- function(...) {
  value <- paste0(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " ",
    paste0(..., collapse = "")
  )
  cat(value, "\n")
  cat(value, "\n", file = log_path, append = TRUE)
}
run_started <- Sys.time()
log_line("Starting partial-confounding scenario ", scenario_id,
         "; analysis_mode=", analysis_mode, "; force=", force,
         "; dry_run=", dry_run, ".")

required_packages <- c(
  "devtools", "digest", "jsonlite", "dplyr", "tibble", "limma",
  "lmtest", "mgcv", "sva", "qcrlscR", "statTarget", "TIGERr",
  "malbacR", "pmartR", "winn"
)
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages,
                                               collapse = ", "),
       call. = FALSE)
}
devtools::load_all(file.path(repo_root, "winn"), quiet = TRUE)
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source(file.path(repo_root, "scripts", "benchmark_helpers.R"), local = FALSE)
source(file.path(repo_root, "scripts", "public_benchmark_method_helpers.R"),
       local = FALSE)
source(file.path(repo_root, "scripts", "weighted_pc_r2.R"), local = FALSE)
source(file.path(repo_root, "scripts", "winn_ablation_helpers.R"),
       local = FALSE)
source(file.path(repo_root, "scripts", "run_order_drift_helpers.R"), local = FALSE)
source(file.path(script_dir, "canonical_cache.R"), local = FALSE)
source(file.path(script_dir, "simulation_core.R"), local = FALSE)
source(file.path(script_dir, "partial_confounding_core.R"), local = FALSE)
source(file.path(script_dir, "partial_confounding_portable_hash.R"),
       local = FALSE)
source(file.path(script_dir, "partial_confounding_pilot_metrics.R"),
       local = FALSE)
source(file.path(script_dir, "partial_confounding_full_metrics.R"),
       local = FALSE)
source(file.path(script_dir, "canonical_seed_method_engine.R"), local = FALSE)
if (winn_only) {
  source(file.path(script_dir, "partial_confounding_winn_only_engine.R"),
         local = FALSE)
}

write_csv <- function(value, path) {
  utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
}
sha_file <- canonical_sha256_file

config_manifest_path <- file.path(config_dir, "config_manifest.csv")
if (!file.exists(config_manifest_path)) {
  stop("Full-grid configuration is absent. Run prepare_partial_confounding_full_grid.R.",
       call. = FALSE)
}
config_manifest <- read.csv(
  config_manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
config_paths <- file.path(config_dir, config_manifest$relative_path)
config_valid <- all(file.exists(config_paths)) &&
  all(as.numeric(file.info(config_paths)$size) == config_manifest$bytes) &&
  identical(unname(vapply(config_paths, sha_file, character(1))),
            unname(config_manifest$sha256))
if (!config_valid) stop("Full-grid configuration manifest validation failed.",
                        call. = FALSE)

scenario_order <- read.csv(
  file.path(config_dir, "scenario_order.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
portable_hashes <- read.csv(
  file.path(config_dir, "portable_matrix_hashes.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
scenario_row <- scenario_order[
  scenario_order$scenario_id == scenario_id, , drop = FALSE
]
if (nrow(scenario_row) != 1L || nrow(scenario_order) != 320L) {
  stop("Scenario is absent or duplicated in the frozen order.",
       call. = FALSE)
}
seed_id <- scenario_row$seed_id[[1L]]
if (winn_only && !seed_id %in% sprintf("CONF%02d", 1:10)) {
  stop("WiNN-only partial-confounding grid is frozen to CONF01-CONF10.",
       call. = FALSE)
}
portable_row <- portable_hashes[
  portable_hashes$scenario_id == scenario_id, , drop = FALSE
]
full_method_parameters <- read.csv(
  file.path(config_dir, "method_parameters.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
metric_dictionary <- read.csv(
  file.path(config_dir, "metric_dictionary.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
study_config <- jsonlite::read_json(
  file.path(config_dir, "study_config.json"), simplifyVector = TRUE
)
portable_policy <- partial_confounding_portable_hash_policy()
portable_config_valid <-
  identical(study_config$portable_reconstruction_hash$schema,
            portable_policy$schema) &&
  as.integer(study_config$portable_reconstruction_hash$digits) ==
    portable_policy$digits &&
  isTRUE(all.equal(
    as.numeric(study_config$portable_reconstruction_hash$quantization_unit),
    portable_policy$quantization_unit,
    tolerance = 0
  ))
hidden_seed_table <- read.csv(
  file.path(config_dir, "hidden_reference_seeds.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
hidden_seed_row <- hidden_seed_table[
  hidden_seed_table$seed_id == seed_id, , drop = FALSE
]
parameter_file_hash <- sha_file(file.path(config_dir, "method_parameters.csv"))
expected_parameter_hash <- trimws(readLines(
  file.path(config_dir, "method_parameters.sha256"), warn = FALSE
)[[1L]])
log_reference_columns <- c(
  "raw_log_with_biology_reference_object_sha256",
  "clean_log_with_biology_reference_object_sha256"
)
valid_sha256 <- function(value) {
  length(value) == 1L && !is.na(value) &&
    grepl("^[0-9a-f]{64}$", value)
}
if (
  nrow(full_method_parameters) != 18L ||
    anyDuplicated(full_method_parameters$method) ||
    nrow(metric_dictionary) != 37L || anyDuplicated(metric_dictionary$metric) ||
    nrow(hidden_seed_row) != 1L || nrow(portable_row) != 1L ||
    !all(log_reference_columns %in% names(portable_row)) ||
    !all(log_reference_columns %in% names(scenario_row)) ||
    !identical(
      portable_row$portable_execution_config_sha256_v3[[1L]],
      scenario_row$portable_execution_config_sha256_v3[[1L]]
    ) ||
    !identical(
      portable_row$raw_log_with_biology_quantized_9_sha256[[1L]],
      scenario_row$raw_log_with_biology_quantized_9_sha256[[1L]]
    ) ||
    !identical(
      portable_row$clean_log_with_biology_quantized_9_sha256[[1L]],
      scenario_row$clean_log_with_biology_quantized_9_sha256[[1L]]
    ) ||
    !identical(
      portable_row$raw_log_with_biology_reference_object_sha256[[1L]],
      scenario_row$raw_log_with_biology_reference_object_sha256[[1L]]
    ) ||
    !identical(
      portable_row$clean_log_with_biology_reference_object_sha256[[1L]],
      scenario_row$clean_log_with_biology_reference_object_sha256[[1L]]
    ) ||
    !valid_sha256(
      portable_row$raw_log_with_biology_reference_object_sha256[[1L]]
    ) ||
    !valid_sha256(
      portable_row$clean_log_with_biology_reference_object_sha256[[1L]]
    ) ||
    !portable_config_valid ||
    !identical(parameter_file_hash, expected_parameter_hash)
) {
  stop("Frozen method, metric, or hidden-reference configuration is invalid.",
       call. = FALSE)
}

target_methods <- if (winn_only) {
  c("Raw", "WINN default (no QC)")
} else {
  full_method_parameters$method
}
method_parameters <- full_method_parameters[
  match(target_methods, full_method_parameters$method), , drop = FALSE
]
if (anyNA(method_parameters$method) ||
    !identical(method_parameters$method, target_methods)) {
  stop("Requested method subset is absent from the frozen method ledger.",
       call. = FALSE)
}
method_count <- nrow(method_parameters)
metric_count <- nrow(metric_dictionary)
parameter_hash <- if (winn_only) {
  canonical_sha256_object(method_parameters)
} else {
  parameter_file_hash
}
completion_schema <- if (winn_only) {
  "partial_confounding_winn_only_scenario_complete_v1"
} else {
  "partial_confounding_scenario_complete_v3"
}

bundle_path <- file.path(result_root, scenario_row$bundle_relative_path[[1L]])
if (!file.exists(bundle_path) ||
    !identical(sha_file(bundle_path), scenario_row$bundle_sha256[[1L]])) {
  stop("Scenario bundle is missing or has the wrong SHA-256.", call. = FALSE)
}
bundle <- readRDS(bundle_path)
base <- generate_canonical_simulation(
  bundle$generation_seeds, include_artifact_matrices = TRUE
)
allocation <- bundle$evaluation_truth$phenotype_allocation
feature_truth <- bundle$evaluation_truth$feature_truth
constructed <- construct_partial_confounding_matrices(
  base, allocation, feature_truth, include_intensities = TRUE
)
scenario_validation <- validate_partial_confounding_scenario(
  base, allocation, feature_truth, constructed
)
x <- constructed$raw_intensity
clean <- constructed$clean_ground_truth
raw_log_identity <- partial_confounding_matrix_identity(
  constructed$raw_log_with_biology,
  scenario_row$raw_log_with_biology_quantized_9_sha256[[1L]],
  portable_row$raw_log_with_biology_reference_object_sha256[[1L]]
)
clean_log_identity <- partial_confounding_matrix_identity(
  constructed$clean_log_with_biology,
  scenario_row$clean_log_with_biology_quantized_9_sha256[[1L]],
  portable_row$clean_log_with_biology_reference_object_sha256[[1L]]
)
raw_log_object_sha256 <- raw_log_identity$observed_object_sha256
clean_log_object_sha256 <- clean_log_identity$observed_object_sha256
raw_log_quantized_sha256 <-
  raw_log_identity$observed_authoritative_quantized_sha256
clean_log_quantized_sha256 <-
  clean_log_identity$observed_authoritative_quantized_sha256
raw_intensity_object_sha256_diagnostic <- canonical_sha256_object(x)
clean_ground_truth_object_sha256_diagnostic <- canonical_sha256_object(clean)
raw_intensity_quantized_9_sha256_diagnostic <-
  partial_confounding_portable_matrix_sha256(x)
clean_ground_truth_quantized_9_sha256_diagnostic <-
  partial_confounding_portable_matrix_sha256(clean)

metadata_evaluation <- constructed$correction_metadata
metadata_evaluation$batch <- metadata_evaluation$plate
metadata_evaluation$within_batch_order <- metadata_evaluation$order_in_plate
metadata_evaluation$is_qc_evaluation <-
  metadata_evaluation$sample_type == "control"
hidden_assignment <- select_canonical_hidden_references(
  base$metadata,
  seed = as.integer(hidden_seed_row$hidden_reference_seed[[1L]]),
  hidden_per_plate = as.integer(hidden_seed_row$hidden_per_plate[[1L]])
)
hidden_ids <- hidden_assignment$sample_id[
  hidden_assignment$is_hidden_reference
]
training_ids <- hidden_assignment$sample_id[
  !hidden_assignment$is_hidden_reference
]
meta_blind <- metadata_evaluation[, c(
  "sample_id", "batch", "run_order", "within_batch_order"
), drop = FALSE]
meta_blind$class <- ifelse(meta_blind$sample_id %in% training_ids,
                           "QC", "Sample")
meta_blind$is_qc <- meta_blind$sample_id %in% training_ids
meta_blind$sample_type <- ifelse(meta_blind$is_qc, "control", "sample")
forbidden_fields <- c(
  "phenotype", "phenotype_numeric", "is_case", "responsive",
  "effect_direction", "phenotype_effect_log", "confounding_axis",
  "nominal_strength"
)

rng_component_map <- c(
  "method_rng::Raw" = "method_Raw",
  "method_rng::ComBat" = "method_ComBat",
  "method_rng::QC-RLSC" = "method_QC_RLSC",
  "method_rng::QC-RFSC" = "method_QC_RFSC",
  "method_rng::TIGER" = "method_TIGER",
  "method_rng::SERRF" = "method_SERRF",
  "method_rng::WINN auto (QC)" = "method_WINN_auto_QC",
  "method_rng::WINN auto-batch (QC)" = "method_WINN_auto_batch_QC",
  "method_rng::WINN fixed (QC identities withheld)" =
    "method_WINN_fixed_no_QC"
)
required_rng <- unique(method_parameters$rng_component)
if (!all(required_rng %in% names(rng_component_map))) {
  stop("No partial-confounding seed mapping for every method RNG component.",
       call. = FALSE)
}
seed_ledger_engine <- data.frame(
  component = required_rng,
  effective_seed = as.integer(vapply(required_rng, function(component) {
    bundle$component_seeds[[rng_component_map[[component]]]]
  }, integer(1))),
  stringsAsFactors = FALSE
)

split_counts <- aggregate(
  cbind(
    hidden = hidden_assignment$is_hidden_reference,
    training = !hidden_assignment$is_hidden_reference
  ),
  list(batch = hidden_assignment$plate), sum
)
reconstruction_hashes <- data.frame(
  scenario_id = scenario_id,
  seed_id = seed_id,
  object = c(
    "raw_log_with_biology", "clean_log_with_biology",
    "raw_intensity", "clean_ground_truth"
  ),
  hash_role = c(
    "authoritative_cross_platform_gate", "authoritative_cross_platform_gate",
    "diagnostic_only", "diagnostic_only"
  ),
  portable_hash_schema = portable_policy$schema,
  quantization_digits = portable_policy$digits,
  equivalence_tolerance = portable_policy$quantization_unit,
  expected_quantized_sha256 = c(
    scenario_row$raw_log_with_biology_quantized_9_sha256[[1L]],
    scenario_row$clean_log_with_biology_quantized_9_sha256[[1L]],
    scenario_row$raw_intensity_quantized_9_sha256_diagnostic[[1L]],
    scenario_row$clean_ground_truth_quantized_9_sha256_diagnostic[[1L]]
  ),
  observed_quantized_sha256 = c(
    raw_log_quantized_sha256, clean_log_quantized_sha256,
    raw_intensity_quantized_9_sha256_diagnostic,
    clean_ground_truth_quantized_9_sha256_diagnostic
  ),
  quantized_hash_match_diagnostic = c(
    raw_log_identity$authoritative_hash_match,
    clean_log_identity$authoritative_hash_match,
    identical(
      raw_intensity_quantized_9_sha256_diagnostic,
      scenario_row$raw_intensity_quantized_9_sha256_diagnostic[[1L]]
    ),
    identical(
      clean_ground_truth_quantized_9_sha256_diagnostic,
      scenario_row$clean_ground_truth_quantized_9_sha256_diagnostic[[1L]]
    )
  ),
  reference_platform_object_sha256 = c(
    portable_row$raw_log_with_biology_reference_object_sha256[[1L]],
    portable_row$clean_log_with_biology_reference_object_sha256[[1L]],
    scenario_row$raw_intensity_sha256[[1L]],
    scenario_row$clean_ground_truth_sha256[[1L]]
  ),
  observed_platform_object_sha256 = c(
    raw_log_object_sha256, clean_log_object_sha256,
    raw_intensity_object_sha256_diagnostic,
    clean_ground_truth_object_sha256_diagnostic
  ),
  exact_object_hash_match_diagnostic = c(
    raw_log_identity$exact_object_hash_match_diagnostic,
    clean_log_identity$exact_object_hash_match_diagnostic,
    identical(raw_intensity_object_sha256_diagnostic,
              scenario_row$raw_intensity_sha256[[1L]]),
    identical(clean_ground_truth_object_sha256_diagnostic,
              scenario_row$clean_ground_truth_sha256[[1L]])
  ),
  is_cross_platform_gate = c(TRUE, TRUE, FALSE, FALSE),
  bundle_file_sha256 = scenario_row$bundle_sha256[[1L]],
  stringsAsFactors = FALSE
)
portable_gate_passed <-
  partial_confounding_v3_reconstruction_gate(reconstruction_hashes)
write_csv(reconstruction_hashes,
          file.path(run_dir, "reconstruction_hashes.csv"))
input_checks <- data.frame(
  check = c(
    "config_manifest", "bundle_file_sha256", "bundle_identity",
    "scenario_reconstruction", "portable_hash_policy",
    "raw_log_with_biology_quantized_9_sha256",
    "clean_log_with_biology_quantized_9_sha256",
    "reference_object_hash_metadata_consistent",
    "platform_object_hashes_recorded",
    "dimensions_and_dimnames", "numeric_finite_inputs", "metadata_order",
    "balanced_design_counts", "five_plates", "hidden_total",
    "training_total", "two_hidden_per_plate", "eight_training_per_plate",
    "phenotype_absent_from_method_metadata", "hidden_labels_absent",
    "method_rng_complete", "effect_design_diagnostics_recorded",
    "parameter_and_metric_freeze"
  ),
  passed = c(
    config_valid,
    identical(sha_file(bundle_path), scenario_row$bundle_sha256[[1L]]),
    identical(bundle$scenario$scenario_id[[1L]], scenario_id) &&
      identical(bundle$scenario_config_sha256,
                scenario_row$scenario_config_sha256[[1L]]),
    all(scenario_validation$passed),
    portable_config_valid,
    raw_log_identity$authoritative_hash_match,
    clean_log_identity$authoritative_hash_match,
    identical(scenario_row$raw_intensity_sha256[[1L]],
              bundle$hashes$raw_intensity_sha256) &&
      identical(scenario_row$clean_ground_truth_sha256[[1L]],
                bundle$hashes$clean_ground_truth_sha256) &&
      identical(
        portable_row$raw_log_with_biology_reference_object_sha256[[1L]],
        scenario_row$raw_log_with_biology_reference_object_sha256[[1L]]
      ) &&
      identical(
        portable_row$clean_log_with_biology_reference_object_sha256[[1L]],
        scenario_row$clean_log_with_biology_reference_object_sha256[[1L]]
      ) &&
      identical(
        portable_row$raw_intensity_reference_object_sha256[[1L]],
        bundle$hashes$raw_intensity_sha256
      ) &&
      identical(
        portable_row$clean_ground_truth_reference_object_sha256[[1L]],
        bundle$hashes$clean_ground_truth_sha256
      ),
    all(grepl("^[0-9a-f]{64}$", c(
      raw_log_object_sha256,
      clean_log_object_sha256,
      portable_row$raw_log_with_biology_reference_object_sha256[[1L]],
      portable_row$clean_log_with_biology_reference_object_sha256[[1L]],
      raw_intensity_object_sha256_diagnostic,
      clean_ground_truth_object_sha256_diagnostic,
      raw_intensity_quantized_9_sha256_diagnostic,
      clean_ground_truth_quantized_9_sha256_diagnostic
    ))),
    identical(dim(x), c(1000L, 500L)) &&
      identical(dim(x), dim(clean)) &&
      identical(dimnames(x), dimnames(clean)),
    is.numeric(x) && is.numeric(clean) && all(is.finite(x)) &&
      all(is.finite(clean)) && all(x > -1) && all(clean > -1),
    identical(colnames(x), meta_blind$sample_id) &&
      identical(metadata_evaluation$run_order, seq_len(500L)),
    sum(allocation$sample_type == "sample") == 450L &&
      sum(allocation$sample_type == "control") == 50L &&
      sum(allocation$phenotype == "case", na.rm = TRUE) == 225L &&
      sum(allocation$phenotype == "control", na.rm = TRUE) == 225L,
    length(unique(meta_blind$batch)) == 5L,
    length(hidden_ids) == 10L,
    length(training_ids) == 40L,
    nrow(split_counts) == 5L && all(split_counts$hidden == 2L),
    nrow(split_counts) == 5L && all(split_counts$training == 8L),
    !any(names(meta_blind) %in% forbidden_fields),
    !any(meta_blind$sample_id[meta_blind$is_qc] %in% hidden_ids),
    !anyNA(seed_ledger_engine$effective_seed) &&
      !anyDuplicated(seed_ledger_engine$component),
    is.finite(scenario_row$design_rank[[1L]]) &&
      is.finite(scenario_row$design_ncol[[1L]]) &&
      scenario_row$design_rank[[1L]] <= scenario_row$design_ncol[[1L]] &&
      is.finite(scenario_row$design_condition_number[[1L]]),
    nrow(method_parameters) == method_count &&
      nrow(metric_dictionary) == metric_count &&
      identical(parameter_file_hash, expected_parameter_hash)
  ),
  detail = "",
  stringsAsFactors = FALSE
)
input_checks$detail <- c(
  paste(nrow(config_manifest), "frozen config files"),
  scenario_row$bundle_sha256[[1L]], scenario_id,
  paste(sum(scenario_validation$passed), "/",
        nrow(scenario_validation), sep = ""),
  paste0(portable_policy$schema, "; digits=", portable_policy$digits,
         "; tolerance=", portable_policy$quantization_unit),
  paste0(raw_log_quantized_sha256, "; exact-log-object-match-diagnostic=",
         reconstruction_hashes$exact_object_hash_match_diagnostic[[1L]]),
  paste0(clean_log_quantized_sha256, "; exact-log-object-match-diagnostic=",
         reconstruction_hashes$exact_object_hash_match_diagnostic[[2L]]),
  paste0("raw-log-reference=",
         portable_row$raw_log_with_biology_reference_object_sha256[[1L]],
         "; clean-log-reference=",
         portable_row$clean_log_with_biology_reference_object_sha256[[1L]],
         "; raw-intensity-reference=", bundle$hashes$raw_intensity_sha256,
         "; clean-intensity-reference=",
         bundle$hashes$clean_ground_truth_sha256),
  paste0("raw-log-observed=", raw_log_object_sha256,
         "; clean-log-observed=", clean_log_object_sha256,
         "; raw-intensity-observed=", raw_intensity_object_sha256_diagnostic,
         "; clean-intensity-observed=",
         clean_ground_truth_object_sha256_diagnostic),
  paste(dim(x), collapse = "x"), "finite and log1p-valid",
  "sample IDs and global order exact", "450 study (225/225) + 50 QC",
  paste(unique(meta_blind$batch), collapse = ";"),
  as.character(length(hidden_ids)), as.character(length(training_ids)),
  paste(split_counts$hidden, collapse = ";"),
  paste(split_counts$training, collapse = ";"),
  paste(names(meta_blind), collapse = ";"),
  paste(hidden_ids, collapse = ";"),
  paste(seed_ledger_engine$component, collapse = ";"),
  paste0("rank=", scenario_row$design_rank[[1L]], "/",
         scenario_row$design_ncol[[1L]], "; kappa=",
         signif(scenario_row$design_condition_number[[1L]], 6L)),
  paste0(method_count, " methods; ", metric_count,
         " metrics; selected-ledger-sha256=", parameter_hash,
         "; source-ledger-sha256=", parameter_file_hash)
)
write_csv(input_checks, file.path(run_dir, "input_validation.csv"))
write_csv(scenario_validation,
          file.path(run_dir, "scenario_reconstruction_validation.csv"))
if (!all(input_checks$passed)) {
  stop("Input validation failed: ",
       paste(input_checks$check[!input_checks$passed], collapse = ", "),
       call. = FALSE)
}

reference_assignment <- data.frame(
  scenario_id = scenario_id, seed_id = seed_id,
  sample_id = meta_blind$sample_id, batch = meta_blind$batch,
  run_order = meta_blind$run_order,
  within_batch_order = meta_blind$within_batch_order,
  evaluation_role = ifelse(
    meta_blind$sample_id %in% hidden_ids, "hidden_test",
    ifelse(meta_blind$sample_id %in% training_ids,
           "training_control", "study")
  ),
  method_visible_class = meta_blind$class,
  supplied_as_control = meta_blind$sample_id %in% training_ids,
  stringsAsFactors = FALSE
)
write_csv(reference_assignment, file.path(run_dir, "reference_assignment.csv"))
write_csv(split_counts, file.path(run_dir, "reference_split_counts.csv"))
write_csv(meta_blind, file.path(run_dir, "method_facing_metadata.csv"))

analysis_units <- method_parameters$method
qc_methods <- c(
  "QC-RLSC", "QC-RFSC", "TIGER", "SERRF", "WINN auto (QC)",
  "WINN auto-batch (QC)"
)
exposure <- data.frame(
  scenario_id = scenario_id, seed_id = seed_id, method = analysis_units,
  n_training_controls_supplied = ifelse(
    analysis_units %in% qc_methods, length(training_ids), 0L
  ),
  n_hidden_controls_supplied = 0L,
  hidden_labels_visible = FALSE,
  hidden_used_for_tuning = FALSE,
  phenotype_fields_supplied = "",
  phenotype_used_for_correction_or_tuning = FALSE,
  one_output_for_all_endpoints = TRUE,
  hidden_ids_sha256 = canonical_sha256_object(hidden_ids),
  method_metadata_sha256 = canonical_sha256_object(meta_blind),
  stringsAsFactors = FALSE
)
write_csv(exposure, file.path(run_dir, "method_exposure_protocol.csv"))

source_files <- c(
  runner = script_path,
  method_engine = file.path(script_dir, "canonical_seed_method_engine.R"),
  metric_engine = file.path(script_dir,
                            "partial_confounding_full_metrics.R"),
  effect_metrics = file.path(script_dir,
                             "partial_confounding_pilot_metrics.R"),
  cache = file.path(script_dir, "canonical_cache.R"),
  simulation_core = file.path(script_dir, "simulation_core.R"),
  confounding_core = file.path(script_dir, "partial_confounding_core.R"),
  portable_hash = file.path(
    script_dir, "partial_confounding_portable_hash.R"
  ),
  benchmark_helpers = file.path(repo_root, "scripts", "benchmark_helpers.R"),
  public_helpers = file.path(
    repo_root, "scripts", "public_benchmark_method_helpers.R"
  ),
  weighted_r2 = file.path(repo_root, "scripts", "weighted_pc_r2.R"),
  ablation_helpers = file.path(repo_root, "scripts",
                               "winn_ablation_helpers.R"),
  drift_helpers = file.path(repo_root, "scripts", "run_order_drift_helpers.R"),
  winn_core = file.path(repo_root, "winn", "R", "winn.R"),
  winn_ablation = file.path(repo_root, "winn", "R", "winn-ablation.R")
)
if (winn_only) {
  source_files <- c(
    source_files,
    winn_only_method_engine = file.path(
      script_dir, "partial_confounding_winn_only_engine.R"
    )
  )
}
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
winn_commit <- if (length(winn_commit_output)) {
  winn_commit_output[[1L]]
} else {
  NA_character_
}
winn_status <- winn_git(c("status", "--short"))
if (is.na(winn_commit) || !nzchar(winn_commit)) {
  stop("Unable to record the local WiNN Git commit.", call. = FALSE)
}
method_context <- list(
  schema = if (winn_only) {
    "partial_confounding_winn_only_phenotype_blind_method_context_v1"
  } else {
    "partial_confounding_phenotype_blind_method_context_v2"
  },
  analysis_mode = analysis_mode,
  seed_id = seed_id,
  raw_intensity_object_sha256_diagnostic =
    raw_intensity_object_sha256_diagnostic,
  raw_intensity_quantized_9_sha256_diagnostic =
    raw_intensity_quantized_9_sha256_diagnostic,
  method_metadata_sha256 = canonical_sha256_object(meta_blind),
  frozen_parameters_sha256 = parameter_hash,
  source_sha256 = as.list(source_hashes),
  code_bundle_sha256 = code_bundle_sha256,
  winn_commit = winn_commit,
  winn_status_sha256 = canonical_sha256_object(winn_status),
  training_ids_sha256 = canonical_sha256_object(training_ids),
  hidden_ids_sha256 = canonical_sha256_object(hidden_ids),
  method_seed_ledger_sha256 = canonical_sha256_object(seed_ledger_engine),
  label_blinding = paste(
    "40 training QCs labeled QC; 10 hidden references labeled ordinary",
    "Sample; phenotype and confounding fields absent"
  )
)
method_context_sha256 <- canonical_sha256_object(method_context)
analysis_provenance <- list(
  schema = if (winn_only) {
    "partial_confounding_winn_only_scenario_run_v1"
  } else {
    "partial_confounding_scenario_run_v3"
  },
  analysis_mode = analysis_mode,
  seed_id = seed_id,
  scenario_id = scenario_id,
  scenario_config_sha256 = scenario_row$scenario_config_sha256[[1L]],
  portable_execution_config_sha256_v3 =
    scenario_row$portable_execution_config_sha256_v3[[1L]],
  bundle_file_sha256 = scenario_row$bundle_sha256[[1L]],
  raw_log_with_biology_object_sha256 = raw_log_object_sha256,
  clean_log_with_biology_object_sha256 = clean_log_object_sha256,
  reference_raw_log_with_biology_object_sha256_diagnostic =
    scenario_row$raw_log_with_biology_reference_object_sha256[[1L]],
  reference_clean_log_with_biology_object_sha256_diagnostic =
    scenario_row$clean_log_with_biology_reference_object_sha256[[1L]],
  raw_log_with_biology_quantized_9_sha256 = raw_log_quantized_sha256,
  clean_log_with_biology_quantized_9_sha256 = clean_log_quantized_sha256,
  raw_intensity_object_sha256_diagnostic =
    raw_intensity_object_sha256_diagnostic,
  clean_ground_truth_object_sha256_diagnostic =
    clean_ground_truth_object_sha256_diagnostic,
  raw_intensity_quantized_9_sha256_diagnostic =
    raw_intensity_quantized_9_sha256_diagnostic,
  clean_ground_truth_quantized_9_sha256_diagnostic =
    clean_ground_truth_quantized_9_sha256_diagnostic,
  portable_hash_schema = portable_policy$schema,
  portable_hash_digits = portable_policy$digits,
  portable_equivalence_tolerance = portable_policy$quantization_unit,
  reference_raw_intensity_object_sha256_diagnostic =
    scenario_row$raw_intensity_sha256[[1L]],
  reference_clean_ground_truth_object_sha256_diagnostic =
    scenario_row$clean_ground_truth_sha256[[1L]],
  phenotype_allocation_sha256 =
    bundle$hashes$phenotype_allocation_sha256,
  responsive_truth_sha256 =
    bundle$hashes$responsive_feature_truth_sha256,
  technical_artifact_sha256 =
    bundle$hashes$technical_artifact_log_sha256,
  frozen_parameters_sha256 = parameter_hash,
  config_sha256 = stats::setNames(
    config_manifest$sha256, config_manifest$relative_path
  ),
  source_sha256 = as.list(source_hashes),
  code_bundle_sha256 = code_bundle_sha256,
  winn_commit = winn_commit,
  winn_status_sha256 = canonical_sha256_object(winn_status),
  method_context_sha256 = method_context_sha256
)
analysis_context_sha256 <- canonical_sha256_object(analysis_provenance)

forbidden_context_fields <- c(
  "scenario_id", "scenario_config_sha256", "bundle_file_sha256",
  "clean_log_with_biology_object_sha256",
  "reference_clean_log_with_biology_object_sha256_diagnostic",
  "clean_log_with_biology_quantized_9_sha256",
  "clean_ground_truth_object_sha256_diagnostic",
  "clean_ground_truth_quantized_9_sha256_diagnostic",
  "reference_clean_ground_truth_object_sha256_diagnostic",
  "phenotype_allocation_sha256",
  "responsive_truth_sha256", "technical_artifact_sha256",
  "confounding_axis", "nominal_strength", "phenotype"
)
recursive_context_names <- function(value) {
  if (!is.list(value)) return(character())
  c(names(value), unlist(lapply(value, recursive_context_names),
                         use.names = FALSE))
}
context_names <- recursive_context_names(method_context)
context_values <- as.character(unlist(
  method_context, recursive = TRUE, use.names = FALSE
))
evaluation_only_hashes <- unique(c(
  scenario_row$scenario_config_sha256[[1L]],
  scenario_row$bundle_sha256[[1L]], clean_log_object_sha256,
  scenario_row$clean_log_with_biology_reference_object_sha256[[1L]],
  scenario_row$clean_ground_truth_sha256[[1L]], clean_log_quantized_sha256,
  clean_ground_truth_object_sha256_diagnostic,
  clean_ground_truth_quantized_9_sha256_diagnostic,
  bundle$hashes$phenotype_allocation_sha256,
  bundle$hashes$responsive_feature_truth_sha256,
  bundle$hashes$technical_artifact_log_sha256
))
evaluation_only_hashes <- evaluation_only_hashes[
  !is.na(evaluation_only_hashes) & nzchar(evaluation_only_hashes)
]
forbidden_name_present <- any(forbidden_context_fields %in% context_names)
evaluation_hash_present <- any(context_values %in% evaluation_only_hashes)
cache_context_check <- data.frame(
  check = c(
    "phenotype_and_truth_fields_absent",
    "scenario_and_confounding_fields_absent",
    "evaluation_only_hash_values_absent",
    "raw_method_metadata_split_parameters_sources_present"
  ),
  passed = c(
    !any(c("clean_log_with_biology_object_sha256",
           "clean_log_with_biology_quantized_9_sha256",
           "reference_clean_log_with_biology_object_sha256_diagnostic",
           "clean_ground_truth_object_sha256_diagnostic",
           "clean_ground_truth_quantized_9_sha256_diagnostic",
           "reference_clean_ground_truth_object_sha256_diagnostic",
           "phenotype_allocation_sha256",
           "responsive_truth_sha256", "phenotype") %in%
           context_names),
    !any(c("scenario_id", "scenario_config_sha256", "confounding_axis",
           "nominal_strength") %in% context_names),
    !evaluation_hash_present,
    all(c("raw_intensity_object_sha256_diagnostic",
          "raw_intensity_quantized_9_sha256_diagnostic",
          "method_metadata_sha256",
          "training_ids_sha256", "hidden_ids_sha256",
          "frozen_parameters_sha256", "source_sha256",
          "method_seed_ledger_sha256") %in% context_names)
  ),
  method_context_sha256 = method_context_sha256,
  stringsAsFactors = FALSE
)
if (!all(cache_context_check$passed) ||
    forbidden_name_present || evaluation_hash_present) {
  stop("Phenotype-blind method cache-context validation failed.",
       call. = FALSE)
}
jsonlite::write_json(
  method_context, file.path(run_dir, "method_cache_context.json"),
  auto_unbox = TRUE, pretty = TRUE, na = "null"
)
write_csv(cache_context_check,
          file.path(run_dir, "method_cache_context_validation.csv"))
exposure$method_cache_context_sha256 <- method_context_sha256
exposure$phenotype_or_truth_hash_in_method_cache_context <-
  forbidden_name_present || evaluation_hash_present
write_csv(exposure, file.path(run_dir, "method_exposure_protocol.csv"))

input_checksum_paths <- c(
  unname(source_files), config_paths, bundle_path,
  file.path(result_root, "scenario_manifest.csv"),
  file.path(result_root, "seed_ledger.csv")
)
input_checksums <- data.frame(
  role = c(
    names(source_files), paste0("config::", config_manifest$relative_path),
    "scenario_bundle", "global_scenario_manifest", "global_seed_ledger"
  ),
  path = input_checksum_paths,
  bytes = as.numeric(file.info(input_checksum_paths)$size),
  sha256 = vapply(input_checksum_paths, sha_file, character(1)),
  stringsAsFactors = FALSE
)
input_checksums$path <- vapply(input_checksums$path, function(path) {
  prefix <- paste0(repo_root, .Platform$file.sep)
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else path
}, character(1))
write_csv(input_checksums, file.path(run_dir, "input_checksums.csv"))

completion_path <- file.path(canonical_run_dir, "analysis_complete.json")
artifact_path <- file.path(canonical_run_dir, "artifact_checksums.csv")
validate_completed_run <- function() {
  if (!file.exists(completion_path) || !file.exists(artifact_path)) {
    return(FALSE)
  }
  completion <- tryCatch(
    jsonlite::read_json(completion_path, simplifyVector = TRUE),
    error = function(e) NULL
  )
  completion_valid <- !is.null(completion) &&
    identical(completion$schema, completion_schema) &&
    identical(completion$scenario_id, scenario_id) &&
    identical(completion$analysis_context_sha256,
              analysis_context_sha256) &&
    isTRUE(completion$invariants_passed) &&
    isTRUE(completion$metrics_complete) &&
    isTRUE(completion$full_gam_complete) &&
    isTRUE(completion$detail_rows_complete) &&
    isTRUE(completion$method_cache_context_blinding_passed) &&
    isTRUE(completion$non_estimability_explicit) &&
    isTRUE(completion$portable_reconstruction_hashes_passed) &&
    identical(completion$portable_execution_config_sha256_v3,
              scenario_row$portable_execution_config_sha256_v3[[1L]]) &&
    identical(completion$raw_log_with_biology_quantized_9_sha256,
              raw_log_quantized_sha256) &&
    identical(completion$clean_log_with_biology_quantized_9_sha256,
              clean_log_quantized_sha256) &&
    identical(completion$raw_log_with_biology_object_sha256,
              raw_log_object_sha256) &&
    identical(completion$clean_log_with_biology_object_sha256,
              clean_log_object_sha256) &&
    identical(
      completion$reference_raw_log_with_biology_object_sha256_diagnostic,
      scenario_row$raw_log_with_biology_reference_object_sha256[[1L]]
    ) &&
    identical(
      completion$reference_clean_log_with_biology_object_sha256_diagnostic,
      scenario_row$clean_log_with_biology_reference_object_sha256[[1L]]
    ) &&
    identical(completion$raw_intensity_object_sha256_diagnostic,
              raw_intensity_object_sha256_diagnostic) &&
    identical(completion$clean_ground_truth_object_sha256_diagnostic,
              clean_ground_truth_object_sha256_diagnostic) &&
    identical(completion$method_context_sha256,
              method_context_sha256) &&
    as.integer(completion$method_count) == method_count &&
    as.integer(completion$metric_rows) == method_count * metric_count
  if (!completion_valid) return(FALSE)
  artifacts <- tryCatch(
    read.csv(artifact_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  required_artifacts <- c(
    "analysis_complete.json", "run_manifest.csv", "method_status.csv",
    "errors.csv", "runtime.csv", "method_metrics.csv",
    "validation_checks.csv", "method_matrix_checks.csv",
    "matrix_checksums.csv", "metric_equivalence_checks.csv",
    "method_exposure_protocol.csv", "method_cache_context.json",
    "method_cache_context_validation.csv", "selected_parameters.csv",
    "metric_detail_row_counts.csv", "post_metric_validation.csv",
    "reconstruction_hashes.csv",
    "metrics/effect_design_diagnostics.csv",
    "metrics/run_order_gam_by_feature_batch.csv.gz", "sessionInfo.txt"
  )
  if (is.null(artifacts) ||
      anyDuplicated(artifacts$relative_path) ||
      any(grepl("(^/|(^|/)\\.\\.(/|$))", artifacts$relative_path)) ||
      !all(required_artifacts %in% artifacts$relative_path)) return(FALSE)
  paths <- file.path(canonical_run_dir, artifacts$relative_path)
  disk_files <- list.files(
    canonical_run_dir, recursive = TRUE, full.names = TRUE
  )
  disk_files <- disk_files[!grepl("(^|/)cache(/|$)", disk_files)]
  disk_files <- setdiff(disk_files, artifact_path)
  disk_files <- disk_files[file.info(disk_files)$isdir %in% FALSE]
  disk_relative <- substring(disk_files, nchar(canonical_run_dir) + 2L)
  setequal(disk_relative, artifacts$relative_path) &&
    all(file.exists(paths)) &&
    all(as.numeric(file.info(paths)$size) == artifacts$bytes) &&
    identical(unname(vapply(paths, sha_file, character(1))),
              unname(artifacts$sha256))
}

if (dry_run) {
  dry <- data.frame(
    scenario_id = scenario_id, seed_id = seed_id,
    bundle_and_reconstruction_valid = all(input_checks$passed),
    package_dependencies_available = !length(missing_packages),
    frozen_config_valid = config_valid,
    hidden_reference_protocol_valid = length(hidden_ids) == 10L &&
      length(training_ids) == 40L &&
      !any(meta_blind$sample_id[meta_blind$class == "QC"] %in% hidden_ids),
    phenotype_blinding_valid = !any(names(meta_blind) %in% forbidden_fields),
    method_cache_context_blinding_valid = all(cache_context_check$passed),
    portable_reconstruction_hashes_valid = portable_gate_passed,
    raw_log_exact_object_hash_match_diagnostic =
      reconstruction_hashes$exact_object_hash_match_diagnostic[[1L]],
    clean_log_exact_object_hash_match_diagnostic =
      reconstruction_hashes$exact_object_hash_match_diagnostic[[2L]],
    raw_intensity_exact_object_hash_match_diagnostic =
      reconstruction_hashes$exact_object_hash_match_diagnostic[[3L]],
    clean_ground_truth_exact_object_hash_match_diagnostic =
      reconstruction_hashes$exact_object_hash_match_diagnostic[[4L]],
    raw_log_with_biology_quantized_9_sha256 =
      raw_log_quantized_sha256,
    clean_log_with_biology_quantized_9_sha256 =
      clean_log_quantized_sha256,
    raw_intensity_quantized_9_sha256_diagnostic =
      raw_intensity_quantized_9_sha256_diagnostic,
    clean_ground_truth_quantized_9_sha256_diagnostic =
      clean_ground_truth_quantized_9_sha256_diagnostic,
    existing_completed_run_reusable = validate_completed_run(),
    analysis_context_sha256 = analysis_context_sha256,
    stringsAsFactors = FALSE
  )
  dry_dir <- file.path(full_root, "validation", "dry_run")
  dir.create(dry_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(dry, file.path(dry_dir, paste0(scenario_id, ".csv")))
  message("Dry-run validation passed for ", scenario_id,
          "; no correction method was executed.")
  cleanup_reuse_validation_scratch()
  quit(save = "no", status = 0L)
}

if (!force && validate_completed_run()) {
  message(scenario_id,
          " already has a complete content-matched run; reusing it.")
  cleanup_reuse_validation_scratch()
  quit(save = "no", status = 0L)
}
canonical_has_files <- dir.exists(canonical_run_dir) &&
  length(list.files(canonical_run_dir, all.files = TRUE, no.. = TRUE)) > 0L
if (canonical_has_files && !force) {
  cleanup_reuse_validation_scratch()
  stop(
    scenario_id,
    " has v1, failed, incomplete, or otherwise non-reusable output; use ",
    "--force to archive the entire run before v3 recomputation.",
    call. = FALSE
  )
}
if (canonical_has_files && force) {
  archive_parent <- file.path(full_root, "archive")
  dir.create(archive_parent, recursive = TRUE, showWarnings = FALSE)
  archive_dir <- file.path(
    archive_parent,
    paste0(scenario_id, "_", format(Sys.time(), "%Y%m%dT%H%M%S"),
           "_", substr(analysis_context_sha256, 1L, 12L))
  )
  if (file.exists(archive_dir) ||
      !file.rename(canonical_run_dir, archive_dir)) {
    stop("Could not archive the prior scenario run.", call. = FALSE)
  }
}

if (nzchar(reuse_validation_scratch)) {
  cleanup_reuse_validation_scratch()
  stop("Internal completed-run reuse validation state reached execution.",
       call. = FALSE)
}

run_dir <- canonical_run_dir
dir.create(file.path(run_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(run_dir, "logs", "analysis.log")
writeLines(character(), log_path)
log_line("Beginning correction and full metrics for ", scenario_id, ".")
preflight_files <- c(
  "input_validation.csv", "scenario_reconstruction_validation.csv",
  "reference_assignment.csv", "reference_split_counts.csv",
  "method_facing_metadata.csv", "method_exposure_protocol.csv",
  "method_cache_context.json", "method_cache_context_validation.csv",
  "input_checksums.csv", "reconstruction_hashes.csv"
)
for (audit_file in preflight_files) {
  if (!file.copy(
    file.path(preflight_run_dir, audit_file), file.path(run_dir, audit_file),
    overwrite = TRUE
  )) {
    stop("Could not stage preflight artifact: ", audit_file,
         call. = FALSE)
  }
}

run_canonical_cache_negative_tests(
  file.path(run_dir, "cache_key_negative_tests.csv")
)
execution <- if (winn_only) {
  run_partial_confounding_winn_only_methods(
    x = x, meta_blind = meta_blind,
    seed_ledger = seed_ledger_engine, method_parameters = method_parameters,
    run_dir = run_dir, global_context = method_context,
    force = force, log_line = log_line
  )
} else {
  run_canonical_seed_methods(
    x = x, meta_blind = meta_blind, training_ids = training_ids,
    seed_ledger = seed_ledger_engine, method_parameters = method_parameters,
    run_dir = run_dir, global_context = method_context,
    force = force, log_line = log_line
  )
}

add_scenario <- function(value) {
  value$scenario_id <- rep(scenario_id, nrow(value))
  value[, c("scenario_id", setdiff(names(value), "scenario_id")),
        drop = FALSE]
}
for (field in c(
  "status", "runtime", "errors", "cache_manifest", "warnings_messages",
  "matrix_manifest"
)) {
  execution[[field]] <- add_scenario(execution[[field]])
}
write_csv(execution$status, file.path(run_dir, "method_status.csv"))
write_csv(execution$runtime, file.path(run_dir, "runtime.csv"))
write_csv(execution$errors, file.path(run_dir, "errors.csv"))
write_csv(execution$cache_manifest, file.path(run_dir, "cache_manifest.csv"))
write_csv(execution$warnings_messages,
          file.path(run_dir, "warnings_messages_errors.csv"))
write_csv(execution$matrix_manifest,
          file.path(run_dir, "matrix_checksums.csv"))

method_order <- method_parameters$method
matrix_checks <- dplyr::bind_rows(lapply(
  names(execution$matrices), function(method) {
    matrix <- execution$matrices[[method]]
    data.frame(
      scenario_id = scenario_id, seed_id = seed_id, method = method,
      numeric_matrix = is.matrix(matrix) && is.numeric(matrix),
      exact_dimensions = identical(dim(matrix), dim(x)),
      exact_feature_order = identical(rownames(matrix), rownames(x)),
      exact_sample_order = identical(colnames(matrix), colnames(x)),
      all_finite = all(is.finite(matrix)),
      log1p_valid = all(matrix > -1),
      stringsAsFactors = FALSE
    )
  }
))
write_csv(matrix_checks, file.path(run_dir, "method_matrix_checks.csv"))

max_difference <- function(a, b) max(abs(a - b), na.rm = TRUE)
has <- function(values) all(values %in% names(execution$matrices))
all_methods_completed <- nrow(execution$status) == method_count &&
  identical(execution$status$method, method_order) &&
  all(grepl("^completed", execution$status$status))
if (winn_only) {
  raw_identity_difference <- if (has("Raw")) {
    max_difference(execution$matrices$Raw, x)
  } else Inf
  invariants <- data.frame(
    scenario_id = scenario_id, seed_id = seed_id,
    check = c(
      "all_2_analysis_units_attempted", "all_2_analysis_units_completed",
      "raw_output_equals_input", "all_matrices_exactly_aligned",
      "same_hidden_references_all_methods",
      "hidden_references_absent_from_controls_and_tuning",
      "phenotype_truth_absent_from_correction_tuning_and_cache",
      "single_output_per_method_for_all_metrics",
      "portable_reconstruction_hashes_match",
      "scenario_bundle_identity_recorded"
    ),
    observed = c(
      nrow(execution$status),
      sum(grepl("^completed", execution$status$status)),
      raw_identity_difference, nrow(matrix_checks),
      length(unique(exposure$hidden_ids_sha256)),
      sum(exposure$n_hidden_controls_supplied),
      sum(exposure$phenotype_used_for_correction_or_tuning),
      all(exposure$one_output_for_all_endpoints),
      sum(reconstruction_hashes$quantized_hash_match_diagnostic[
        reconstruction_hashes$is_cross_platform_gate
      ]),
      scenario_row$portable_execution_config_sha256_v3[[1L]]
    ),
    tolerance = c(NA, NA, 1e-8, rep(NA, 7L)),
    stringsAsFactors = FALSE
  )
  invariants$passed <- c(
    nrow(execution$status) == method_count &&
      all(execution$status$attempted),
    all_methods_completed,
    raw_identity_difference <= 1e-8,
    nrow(matrix_checks) == method_count &&
      all(as.matrix(matrix_checks[, -(1:3), drop = FALSE])),
    length(unique(exposure$hidden_ids_sha256)) == 1L,
    all(exposure$n_hidden_controls_supplied == 0L &
          !exposure$hidden_labels_visible &
          !exposure$hidden_used_for_tuning),
    all(!exposure$phenotype_used_for_correction_or_tuning &
          exposure$phenotype_fields_supplied == "" &
          !exposure$phenotype_or_truth_hash_in_method_cache_context),
    all(exposure$one_output_for_all_endpoints),
    portable_gate_passed,
    nzchar(scenario_row$portable_execution_config_sha256_v3[[1L]])
  )
} else {
  invariants <- data.frame(
    scenario_id = scenario_id, seed_id = seed_id,
    check = c(
      "all_18_analysis_units_attempted", "all_18_analysis_units_completed",
      "raw_equals_c0", "fixed_equals_c4", "c4_equals_gss",
      "fixed_cross_wrapper_equals_independent_ablation",
      "all_matrices_exactly_aligned", "same_hidden_references_all_methods",
      "hidden_references_absent_from_controls_and_tuning",
      "phenotype_truth_absent_from_correction_tuning_and_cache",
      "single_output_per_method_for_all_metrics",
      "portable_reconstruction_hashes_match",
      "scenario_bundle_identity_recorded"
    ),
    observed = c(
      nrow(execution$status),
      sum(grepl("^completed", execution$status$status)),
      if (has(c("Raw", "C0_RAW"))) {
        max_difference(execution$matrices$Raw, execution$matrices$C0_RAW)
      } else Inf,
      if (has(c("WINN default (no QC)", "C4_FULL_FIXED"))) {
        max_difference(execution$matrices[["WINN default (no QC)"]],
                       execution$matrices$C4_FULL_FIXED)
      } else Inf,
      if (has(c("C4_FULL_FIXED", "G_SS"))) {
        max_difference(execution$matrices$C4_FULL_FIXED,
                       execution$matrices$G_SS)
      } else Inf,
      if (has(c("WINN default (no QC)", "C4_FULL_FIXED"))) {
        max_difference(execution$matrices[["WINN default (no QC)"]],
                       execution$matrices$C4_FULL_FIXED)
      } else Inf,
      nrow(matrix_checks), length(unique(exposure$hidden_ids_sha256)),
      sum(exposure$n_hidden_controls_supplied),
      sum(exposure$phenotype_used_for_correction_or_tuning),
      all(exposure$one_output_for_all_endpoints),
      sum(reconstruction_hashes$quantized_hash_match_diagnostic[
        reconstruction_hashes$is_cross_platform_gate
      ]),
      scenario_row$portable_execution_config_sha256_v3[[1L]]
    ),
    tolerance = c(NA, NA, rep(1e-8, 4L), rep(NA, 7L)),
    stringsAsFactors = FALSE
  )
  invariants$passed <- c(
    nrow(execution$status) == method_count && all(execution$status$attempted),
    all_methods_completed,
    as.numeric(invariants$observed[[3L]]) <= 1e-8,
    as.numeric(invariants$observed[[4L]]) <= 1e-8,
    as.numeric(invariants$observed[[5L]]) <= 1e-8,
    as.numeric(invariants$observed[[6L]]) <= 1e-8,
    nrow(matrix_checks) == method_count &&
      all(as.matrix(matrix_checks[, -(1:3), drop = FALSE])),
    length(unique(exposure$hidden_ids_sha256)) == 1L,
    all(exposure$n_hidden_controls_supplied == 0L &
          !exposure$hidden_labels_visible &
          !exposure$hidden_used_for_tuning),
    all(!exposure$phenotype_used_for_correction_or_tuning &
          exposure$phenotype_fields_supplied == "" &
          !exposure$phenotype_or_truth_hash_in_method_cache_context),
    all(exposure$one_output_for_all_endpoints),
    portable_gate_passed,
    nzchar(scenario_row$portable_execution_config_sha256_v3[[1L]])
  )
}
write_csv(invariants, file.path(run_dir, "validation_checks.csv"))
if (!all(invariants$passed)) {
  jsonlite::write_json(list(
    schema = if (winn_only) {
      "partial_confounding_winn_only_scenario_failure_v1"
    } else {
      "partial_confounding_scenario_failure_v3"
    },
    scenario_id = scenario_id, seed_id = seed_id,
    analysis_context_sha256 = analysis_context_sha256,
    failed_checks = invariants$check[!invariants$passed],
    method_failures = execution$status$method[
      !grepl("^completed", execution$status$status)
    ],
    stopped_before_metrics = TRUE,
    recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ), file.path(run_dir, "run_failure.json"),
  auto_unbox = TRUE, pretty = TRUE, na = "null")
  stop(scenario_id, " failed invariants; metrics were not run.",
       call. = FALSE)
}

metrics <- compute_partial_confounding_full_metrics(
  scenario_id = scenario_id, seed_id = seed_id,
  method_matrices = execution$matrices, raw = x, clean = clean,
  allocation = allocation, metadata_evaluation = metadata_evaluation,
  hidden_ids = hidden_ids, feature_truth = feature_truth,
  runtime = execution$runtime,
  ablation_diagnostics = execution$ablation_diagnostics,
  metric_dictionary = metric_dictionary, run_dir = run_dir,
  log_line = log_line
)
write_csv(metrics$detailed_row_counts,
          file.path(run_dir, "metric_detail_row_counts.csv"))
expected_detail_rows <- c(
  phenotype_effect_estimates = method_count * 1000L,
  heldout_qc_cv_by_feature = method_count * 1000L,
  run_order_gam_by_feature_batch = method_count * 1000L * 5L,
  ljung_box_by_feature_batch = method_count * 1000L * 5L,
  clean_recovery_by_feature = method_count * 1000L,
  clean_recovery_by_sample = method_count * 450L
)
observed_detail_rows <- stats::setNames(
  metrics$detailed_row_counts$rows, metrics$detailed_row_counts$table
)
detail_rows_complete <- identical(
  unname(observed_detail_rows[names(expected_detail_rows)]),
  unname(expected_detail_rows)
)

metric_values <- function(method) {
  rows <- metrics$method_metrics[
    metrics$method_metrics$method == method, , drop = FALSE
  ]
  rows$value[match(metric_dictionary$metric, rows$metric)]
}
compare_metric_alias <- function(check, method_a, method_b) {
  value_a <- metric_values(method_a)
  value_b <- metric_values(method_b)
  evaluated <- metric_dictionary$metric != "runtime_sec"
  same_na <- identical(is.na(value_a[evaluated]), is.na(value_b[evaluated]))
  finite <- evaluated & is.finite(value_a) & is.finite(value_b)
  maximum <- if (any(finite)) {
    max(abs(value_a[finite] - value_b[finite]))
  } else {
    0
  }
  data.frame(
    scenario_id = scenario_id, seed_id = seed_id, check = check,
    method_a = method_a, method_b = method_b,
    metrics_compared = sum(evaluated), same_na_pattern = same_na,
    max_absolute_difference = maximum, tolerance = 1e-12,
    passed = same_na && is.finite(maximum) && maximum <= 1e-12,
    excluded_metric = "runtime_sec", stringsAsFactors = FALSE
  )
}
metric_equivalences <- if (winn_only) {
  data.frame(
    scenario_id = scenario_id, seed_id = seed_id,
    check = "ablation_alias_checks_not_applicable_winn_only",
    method_a = "", method_b = "", metrics_compared = 0L,
    same_na_pattern = TRUE, max_absolute_difference = 0,
    tolerance = 0, passed = TRUE, excluded_metric = "all",
    stringsAsFactors = FALSE
  )
} else {
  dplyr::bind_rows(
    compare_metric_alias("raw_equals_c0_metrics", "Raw", "C0_RAW"),
    compare_metric_alias(
      "fixed_equals_c4_metrics", "WINN default (no QC)", "C4_FULL_FIXED"
    ),
    compare_metric_alias("c4_equals_gss_metrics", "C4_FULL_FIXED", "G_SS")
  )
}
write_csv(metric_equivalences,
          file.path(run_dir, "metric_equivalence_checks.csv"))

gam_rows <- metrics$method_metrics[
  metrics$method_metrics$metric ==
    "residual_run_order_gam_mean_deviance", , drop = FALSE
]
effect_designs <- metrics$effect_designs
metrics_complete <- nrow(metrics$method_metrics) == method_count * metric_count &&
  setequal(unique(metrics$method_metrics$metric), metric_dictionary$metric)
full_gam_complete <- nrow(gam_rows) == method_count &&
  all(gam_rows$status == "calculated_full_feature_batch_profiles") &&
  all(is.finite(gam_rows$value))
non_estimability_explicit <- nrow(effect_designs) == method_count &&
  all(c("estimable", "rank", "n_columns", "condition_number") %in%
        names(effect_designs)) &&
  all(metrics$method_metrics$status[
    metrics$method_metrics$metric %in% c(
      "median_attenuation_ratio_clean_responsive",
      "effect_rmse_vs_clean_responsive", "true_positive_rate"
    )
  ] %in% c("calculated", "not_estimable_rank_deficient"))
post_metric_checks <- data.frame(
  scenario_id = scenario_id, seed_id = seed_id,
  check = c(
    "metric_alias_equalities", paste0("metrics_complete_", method_count,
                                       "x", metric_count),
    "full_gam_complete", "full_detail_row_counts",
    "non_estimability_explicit"
  ),
  passed = c(all(metric_equivalences$passed), metrics_complete,
             full_gam_complete, detail_rows_complete,
             non_estimability_explicit),
  stringsAsFactors = FALSE
)
write_csv(post_metric_checks,
          file.path(run_dir, "post_metric_validation.csv"))
if (!all(post_metric_checks$passed)) {
  jsonlite::write_json(list(
    schema = if (winn_only) {
      "partial_confounding_winn_only_scenario_failure_v1"
    } else {
      "partial_confounding_scenario_failure_v3"
    },
    scenario_id = scenario_id, seed_id = seed_id,
    analysis_context_sha256 = analysis_context_sha256,
    failed_checks = post_metric_checks$check[!post_metric_checks$passed],
    stopped_after_metrics = TRUE,
    recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ), file.path(run_dir, "run_failure.json"),
  auto_unbox = TRUE, pretty = TRUE, na = "null")
  stop(scenario_id, " failed post-metric validation.", call. = FALSE)
}

extract_message <- function(text, pattern, default = NA_character_) {
  match <- regexec(pattern, text, perl = TRUE)
  found <- regmatches(text, match)[[1L]]
  if (length(found) >= 2L) trimws(found[[2L]]) else default
}
selected_parameters <- lapply(seq_len(nrow(method_parameters)), function(i) {
  parameter <- method_parameters[i, , drop = FALSE]
  method <- parameter$method[[1L]]
  text <- execution$warnings_messages$messages[
    match(method, execution$warnings_messages$method)
  ]
  if (!length(text) || is.na(text)) text <- ""
  selected_json <- parameter$parameters_json[[1L]]
  source_label <- "frozen before CONF01 pilot and full-grid execution"
  if (method %in% c("WINN auto (QC)", "WINN auto-batch (QC)")) {
    selected_json <- canonical_seed_plain_json(list(
      parameters = "auto",
      batch_detection = extract_message(text, "Batch detection: ([^\\n]+)"),
      pelt_penalty = suppressWarnings(as.numeric(extract_message(
        text, "fkPELT penalty: ([0-9.eE+-]+)"
      ))),
      autocorrelation_test = extract_message(
        text, "Autocorrelation test: ([^ (\\n]+)"
      ),
      autocorrelation_fdr = suppressWarnings(as.numeric(extract_message(
        text, "Autocorrelation test: [^\\n]+\\(FDR: ([0-9.]+)\\)"
      ))),
      lag = extract_message(text, "Ljung-Box lag: ([^\\n]+)"),
      spline_method = extract_message(text, "Spline method: ([^\\n]+)"),
      batch_fdr = suppressWarnings(as.numeric(extract_message(
        text, "Batch correction FDR: ([0-9.]+)"
      ))),
      normalization = extract_message(text, "Normalization: ([^\\n]+)"),
      scale_by_batch = identical(toupper(extract_message(
        text, "Scale by batch: ([^\\n]+)"
      )), "TRUE"),
      training_qc_cv = suppressWarnings(as.numeric(extract_message(
        text, "Quality metrics - CV: ([0-9.eE+-]+)"
      ))),
      training_qc_correlation = suppressWarnings(as.numeric(extract_message(
        text, "Correlation: ([0-9.eE+-]+)"
      )))
    ), auto_unbox = TRUE, null = "null", na = "null")
    source_label <- "selected automatically from this scenario's 40 training QCs"
  }
  data.frame(
    scenario_id = scenario_id, seed_id = seed_id, method = method,
    selected_parameters_json = as.character(selected_json),
    selection_source = source_label,
    hidden_references_used = FALSE,
    phenotype_used = FALSE,
    stringsAsFactors = FALSE
  )
})
selected_parameters <- dplyr::bind_rows(selected_parameters)
if (nrow(selected_parameters) != method_count ||
    !is.character(selected_parameters$selected_parameters_json)) {
  stop(
    "Selected-parameter audit table failed its method-count plain-character JSON schema check.",
    call. = FALSE
  )
}
write_csv(selected_parameters,
          file.path(run_dir, "selected_parameters.csv"))
write_csv(data.frame(
  method = "ComBat phenotype-aware (positive control)",
  status = "not_run_no_validated_standard_wrapper",
  included_in_primary_ranking = FALSE,
  reason = paste(
    "The unchanged run_combat() wrapper has no phenotype-covariate",
    "interface; no bespoke method was introduced."
  ),
  stringsAsFactors = FALSE
), file.path(run_dir, "positive_control_status.csv"))

degradation_inputs <- metrics$method_metrics[
  metrics$method_metrics$metric %in% c(
    "median_attenuation_ratio_clean_responsive",
    "effect_rmse_vs_clean_responsive", "true_positive_rate"
  ), c("scenario_id", "seed_id", "method", "metric", "value", "status"),
  drop = FALSE
]
degradation_inputs$paired_zero_required_for_final_flag <-
  degradation_inputs$metric == "true_positive_rate"
degradation_inputs$threshold <- ifelse(
  degradation_inputs$metric ==
    "median_attenuation_ratio_clean_responsive", 0.80,
  ifelse(degradation_inputs$metric ==
           "effect_rmse_vs_clean_responsive", 0.20, 0.20)
)
write_csv(degradation_inputs,
          file.path(run_dir, "degradation_inputs.csv"))

software <- list(
  R = R.version.string, platform = R.version$platform,
  library_paths = .libPaths(),
  project_library = if (nzchar(project_library)) {
    normalizePath(project_library, mustWork = TRUE)
  } else NULL,
  packages = stats::setNames(lapply(
    required_packages,
    function(package) as.character(utils::packageVersion(package))
  ), required_packages),
  source_sha256 = as.list(source_hashes),
  local_winn_commit = winn_commit,
  local_winn_dirty = length(winn_status) > 0L,
  local_winn_status = winn_status
)
jsonlite::write_json(
  software, file.path(run_dir, "software_state.json"),
  auto_unbox = TRUE, pretty = TRUE, na = "null"
)
utils::capture.output(sessionInfo(), file = file.path(run_dir,
                                                      "sessionInfo.txt"))

run_completed <- Sys.time()
manifest_row <- data.frame(
  scenario_id = scenario_id, seed_id = seed_id,
  scenario_index = scenario_row$scenario_index[[1L]],
  array_index = scenario_row$array_index[[1L]],
  confounding_axis = scenario_row$confounding_axis[[1L]],
  nominal_strength = scenario_row$nominal_strength[[1L]],
  realized_batch_cramers_v =
    scenario_row$cramers_v_phenotype_batch[[1L]],
  realized_order_abs_cor_pooled =
    scenario_row$abs_cor_phenotype_within_batch_order[[1L]],
  realized_order_weighted_mean_abs_cor =
    scenario_row$weighted_mean_abs_within_plate_cor[[1L]],
  realized_order_max_abs_cor =
    scenario_row$max_abs_within_plate_cor[[1L]],
  nominal_one_interpretation = ifelse(
    scenario_row$nominal_strength[[1L]] == 1,
    "maximal_feasible_balanced_five_plate_not_complete_global", "not_applicable"
  ),
  status = "completed",
  run_started = format(run_started, "%Y-%m-%dT%H:%M:%S%z"),
  run_completed = format(run_completed, "%Y-%m-%dT%H:%M:%S%z"),
  elapsed_wall_sec = as.numeric(difftime(
    run_completed, run_started, units = "secs"
  )),
  bundle_sha256 = scenario_row$bundle_sha256[[1L]],
  portable_execution_config_sha256_v3 =
    scenario_row$portable_execution_config_sha256_v3[[1L]],
  raw_log_with_biology_object_sha256 = raw_log_object_sha256,
  clean_log_with_biology_object_sha256 = clean_log_object_sha256,
  reference_raw_log_with_biology_object_sha256_diagnostic =
    scenario_row$raw_log_with_biology_reference_object_sha256[[1L]],
  reference_clean_log_with_biology_object_sha256_diagnostic =
    scenario_row$clean_log_with_biology_reference_object_sha256[[1L]],
  raw_log_with_biology_quantized_9_sha256 = raw_log_quantized_sha256,
  clean_log_with_biology_quantized_9_sha256 = clean_log_quantized_sha256,
  raw_intensity_object_sha256_diagnostic =
    raw_intensity_object_sha256_diagnostic,
  clean_ground_truth_object_sha256_diagnostic =
    clean_ground_truth_object_sha256_diagnostic,
  raw_intensity_quantized_9_sha256_diagnostic =
    raw_intensity_quantized_9_sha256_diagnostic,
  clean_ground_truth_quantized_9_sha256_diagnostic =
    clean_ground_truth_quantized_9_sha256_diagnostic,
  reference_raw_intensity_object_sha256_diagnostic =
    scenario_row$raw_intensity_sha256[[1L]],
  reference_clean_ground_truth_object_sha256_diagnostic =
    scenario_row$clean_ground_truth_sha256[[1L]],
  raw_intensity_exact_object_hash_match_diagnostic =
    reconstruction_hashes$exact_object_hash_match_diagnostic[[3L]],
  clean_ground_truth_exact_object_hash_match_diagnostic =
    reconstruction_hashes$exact_object_hash_match_diagnostic[[4L]],
  raw_intensity_round9_hash_match_diagnostic =
    reconstruction_hashes$quantized_hash_match_diagnostic[[3L]],
  clean_ground_truth_round9_hash_match_diagnostic =
    reconstruction_hashes$quantized_hash_match_diagnostic[[4L]],
  portable_hash_schema = portable_policy$schema,
  portable_hash_digits = portable_policy$digits,
  portable_equivalence_tolerance = portable_policy$quantization_unit,
  portable_reconstruction_hashes_passed = portable_gate_passed,
  technical_artifact_sha256 =
    bundle$hashes$technical_artifact_log_sha256,
  responsive_truth_sha256 =
    bundle$hashes$responsive_feature_truth_sha256,
  phenotype_allocation_sha256 =
    bundle$hashes$phenotype_allocation_sha256,
  method_seed_ledger_sha256 = canonical_sha256_object(seed_ledger_engine),
  method_context_sha256 = method_context_sha256,
  hidden_assignment_sha256 = canonical_sha256_object(hidden_assignment),
  frozen_parameters_sha256 = parameter_hash,
  code_bundle_sha256 = code_bundle_sha256,
  winn_commit = winn_commit,
  winn_status_sha256 = canonical_sha256_object(winn_status),
  analysis_context_sha256 = analysis_context_sha256,
  n_features = nrow(x), n_samples = ncol(x), n_batches = 5L,
  n_qc_total = 50L, n_hidden = 10L, n_training = 40L,
  methods_attempted = nrow(execution$status),
  methods_completed = sum(grepl("^completed", execution$status$status)),
  metrics_rows = nrow(metrics$method_metrics),
  invariants_passed = all(invariants$passed),
  metric_equivalences_passed = all(metric_equivalences$passed),
  full_gam_complete = full_gam_complete,
  detail_rows_complete = detail_rows_complete,
  non_estimability_explicit = non_estimability_explicit,
  hidden_exclusion_passed = all(
    exposure$n_hidden_controls_supplied == 0L &
      !exposure$hidden_used_for_tuning
  ),
  phenotype_blinding_passed = all(
    !exposure$phenotype_used_for_correction_or_tuning &
      !exposure$phenotype_or_truth_hash_in_method_cache_context
  ),
  stringsAsFactors = FALSE
)
write_csv(manifest_row, file.path(run_dir, "run_manifest.csv"))
jsonlite::write_json(
  as.list(manifest_row[1, ]), file.path(run_dir, "run_manifest.json"),
  auto_unbox = TRUE, pretty = TRUE, na = "null"
)
completion <- list(
  schema = completion_schema,
  scenario_id = scenario_id, seed_id = seed_id,
  analysis_context_sha256 = analysis_context_sha256,
  scenario_config_sha256 = scenario_row$scenario_config_sha256[[1L]],
  portable_execution_config_sha256_v3 =
    scenario_row$portable_execution_config_sha256_v3[[1L]],
  technical_artifact_sha256 =
    bundle$hashes$technical_artifact_log_sha256,
  responsive_truth_sha256 =
    bundle$hashes$responsive_feature_truth_sha256,
  raw_log_with_biology_object_sha256 = raw_log_object_sha256,
  clean_log_with_biology_object_sha256 = clean_log_object_sha256,
  reference_raw_log_with_biology_object_sha256_diagnostic =
    scenario_row$raw_log_with_biology_reference_object_sha256[[1L]],
  reference_clean_log_with_biology_object_sha256_diagnostic =
    scenario_row$clean_log_with_biology_reference_object_sha256[[1L]],
  raw_log_with_biology_quantized_9_sha256 = raw_log_quantized_sha256,
  clean_log_with_biology_quantized_9_sha256 = clean_log_quantized_sha256,
  raw_intensity_object_sha256_diagnostic =
    raw_intensity_object_sha256_diagnostic,
  clean_ground_truth_object_sha256_diagnostic =
    clean_ground_truth_object_sha256_diagnostic,
  raw_intensity_quantized_9_sha256_diagnostic =
    raw_intensity_quantized_9_sha256_diagnostic,
  clean_ground_truth_quantized_9_sha256_diagnostic =
    clean_ground_truth_quantized_9_sha256_diagnostic,
  reference_raw_intensity_object_sha256_diagnostic =
    scenario_row$raw_intensity_sha256[[1L]],
  reference_clean_ground_truth_object_sha256_diagnostic =
    scenario_row$clean_ground_truth_sha256[[1L]],
  portable_hash_schema = portable_policy$schema,
  portable_hash_digits = portable_policy$digits,
  portable_equivalence_tolerance = portable_policy$quantization_unit,
  raw_intensity_exact_object_hash_match_diagnostic =
    reconstruction_hashes$exact_object_hash_match_diagnostic[[3L]],
  clean_ground_truth_exact_object_hash_match_diagnostic =
    reconstruction_hashes$exact_object_hash_match_diagnostic[[4L]],
  raw_intensity_round9_hash_match_diagnostic =
    reconstruction_hashes$quantized_hash_match_diagnostic[[3L]],
  clean_ground_truth_round9_hash_match_diagnostic =
    reconstruction_hashes$quantized_hash_match_diagnostic[[4L]],
  portable_reconstruction_hashes_passed = TRUE,
  method_seed_ledger_sha256 = canonical_sha256_object(seed_ledger_engine),
  method_context_sha256 = method_context_sha256,
  invariants_passed = TRUE, metrics_complete = TRUE,
  full_gam_complete = TRUE, detail_rows_complete = TRUE,
  method_cache_context_blinding_passed = TRUE,
  non_estimability_explicit = TRUE,
  method_count = method_count, metric_rows = nrow(metrics$method_metrics),
  completed_at = format(run_completed, "%Y-%m-%dT%H:%M:%S%z")
)
jsonlite::write_json(
  completion, completion_path, auto_unbox = TRUE, pretty = TRUE
)
log_line(
  scenario_id, " completed: ", method_count, "/", method_count,
  " analysis units, ",
  nrow(metrics$method_metrics),
  " metric rows, full GAM and all invariants complete; wall_sec=",
  sprintf("%.1f", manifest_row$elapsed_wall_sec), "."
)

# Write this last and do not mutate any listed artifact afterwards.
all_files <- list.files(run_dir, recursive = TRUE, full.names = TRUE)
all_files <- all_files[!grepl("(^|/)cache(/|$)", all_files)]
all_files <- setdiff(all_files, artifact_path)
all_files <- all_files[file.info(all_files)$isdir %in% FALSE]
artifacts <- data.frame(
  relative_path = substring(all_files, nchar(run_dir) + 2L),
  bytes = as.numeric(file.info(all_files)$size),
  sha256 = vapply(all_files, sha_file, character(1)),
  stringsAsFactors = FALSE
)
write_csv(artifacts, artifact_path)
message(scenario_id, " full scenario runner complete; ", nrow(artifacts),
        " artifacts hashed.")
