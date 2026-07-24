#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "scripts/robustness/prepare_partial_confounding_full_grid.R",
    mustWork = TRUE
  )
}
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
result_root <- file.path(
  repo_root, "results", "robustness", "06_partial_confounding"
)
full_root <- file.path(result_root, "full_grid")
config_dir <- file.path(full_root, "config")
dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("digest", "jsonlite")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}
source(file.path(script_dir, "canonical_cache.R"), local = FALSE)
source(file.path(script_dir, "simulation_core.R"), local = FALSE)
source(file.path(script_dir, "partial_confounding_core.R"), local = FALSE)
source(file.path(script_dir, "partial_confounding_portable_hash.R"),
       local = FALSE)
source(file.path(script_dir, "partial_confounding_full_metrics.R"),
       local = FALSE)

write_csv <- function(value, path) {
  utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
}
scenario_manifest_path <- file.path(result_root, "scenario_manifest.csv")
seed_ledger_path <- file.path(result_root, "seed_ledger.csv")
pilot_parameters_path <- file.path(
  result_root, "pilot_CONF01", "method_parameters.csv"
)
pilot_parameters_hash_path <- file.path(
  result_root, "pilot_CONF01", "method_parameters.sha256"
)
for (path in c(
  scenario_manifest_path, seed_ledger_path, pilot_parameters_path,
  pilot_parameters_hash_path
)) {
  if (!file.exists(path)) stop("Missing required input: ", path,
                              call. = FALSE)
}

scenario_manifest <- read.csv(
  scenario_manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
seed_ledger <- read.csv(
  seed_ledger_path, stringsAsFactors = FALSE, check.names = FALSE
)
method_parameters <- read.csv(
  pilot_parameters_path, stringsAsFactors = FALSE, check.names = FALSE
)
expected_parameter_hash <- trimws(readLines(
  pilot_parameters_hash_path, warn = FALSE
)[[1L]])
if (!identical(
  canonical_sha256_file(pilot_parameters_path), expected_parameter_hash
)) {
  stop("Frozen pilot method parameters failed SHA-256 validation.",
       call. = FALSE)
}

expected_seed_ids <- sprintf("CONF%02d", seq_len(20L))
expected_axis <- c("none", "batch", "run_order", "joint")
if (
  nrow(scenario_manifest) != 320L ||
    !identical(unique(scenario_manifest$seed_id), expected_seed_ids) ||
    any(!scenario_manifest$validation_passed) ||
    anyDuplicated(scenario_manifest$scenario_id) ||
    !setequal(unique(scenario_manifest$confounding_axis), expected_axis) ||
    nrow(method_parameters) != 18L || anyDuplicated(method_parameters$method)
) {
  stop("The prepared partial-confounding design is not the expected 320 x 18 grid.",
       call. = FALSE)
}

# Reconstruct every scenario without modifying its frozen bundle. Exact bundle
# file SHA-256 remains mandatory. Cross-platform identity is frozen on the
# pre-exponentiation raw/clean log matrices after 9-decimal quantization.
# Intensity object and intensity-round9 hashes remain diagnostics only.
portable_policy <- partial_confounding_portable_hash_policy()
portable_rows <- vector("list", nrow(scenario_manifest))
v2_raw_log_round10 <- character(nrow(scenario_manifest))
v2_clean_log_round10 <- character(nrow(scenario_manifest))
v2_raw_intensity_round10 <- character(nrow(scenario_manifest))
v2_clean_intensity_round10 <- character(nrow(scenario_manifest))
base_by_seed <- list()
for (index in seq_len(nrow(scenario_manifest))) {
  scenario <- scenario_manifest[index, , drop = FALSE]
  bundle_path <- file.path(result_root, scenario$bundle_relative_path[[1L]])
  if (!file.exists(bundle_path) ||
      !identical(canonical_sha256_file(bundle_path),
                 scenario$bundle_sha256[[1L]])) {
    stop("Frozen bundle file hash failed for ", scenario$scenario_id[[1L]],
         ".", call. = FALSE)
  }
  bundle <- readRDS(bundle_path)
  seed_id <- scenario$seed_id[[1L]]
  if (is.null(base_by_seed[[seed_id]])) {
    base_by_seed[[seed_id]] <- generate_canonical_simulation(
      bundle$generation_seeds, include_artifact_matrices = TRUE
    )
  }
  reconstructed <- construct_partial_confounding_matrices(
    base_by_seed[[seed_id]],
    bundle$evaluation_truth$phenotype_allocation,
    bundle$evaluation_truth$feature_truth,
    include_intensities = TRUE
  )
  reconstruction_checks <- validate_partial_confounding_scenario(
    base_by_seed[[seed_id]],
    bundle$evaluation_truth$phenotype_allocation,
    bundle$evaluation_truth$feature_truth,
    reconstructed
  )
  if (!all(reconstruction_checks$passed)) {
    stop("Reconstruction checks failed for ", scenario$scenario_id[[1L]],
         ".", call. = FALSE)
  }
  raw_log_quantized <- partial_confounding_portable_matrix_sha256(
    reconstructed$raw_log_with_biology
  )
  clean_log_quantized <- partial_confounding_portable_matrix_sha256(
    reconstructed$clean_log_with_biology
  )
  raw_log_reference_object_sha256 <-
    canonical_sha256_object(reconstructed$raw_log_with_biology)
  clean_log_reference_object_sha256 <-
    canonical_sha256_object(reconstructed$clean_log_with_biology)
  raw_intensity_quantized_diagnostic <-
    partial_confounding_portable_matrix_sha256(reconstructed$raw_intensity)
  clean_intensity_quantized_diagnostic <-
    partial_confounding_portable_matrix_sha256(
      reconstructed$clean_ground_truth
    )
  # Retain the superseded v2 round-10 identities only for the immutable
  # attempt ledger. They are not eligible execution gates in v3.
  v2_raw_log_round10[[index]] <-
    partial_confounding_portable_matrix_sha256(
      reconstructed$raw_log_with_biology, digits = 10L
    )
  v2_clean_log_round10[[index]] <-
    partial_confounding_portable_matrix_sha256(
      reconstructed$clean_log_with_biology, digits = 10L
    )
  v2_raw_intensity_round10[[index]] <-
    partial_confounding_portable_matrix_sha256(
      reconstructed$raw_intensity, digits = 10L
    )
  v2_clean_intensity_round10[[index]] <-
    partial_confounding_portable_matrix_sha256(
      reconstructed$clean_ground_truth, digits = 10L
    )
  raw_log_residual <- max(abs(
    reconstructed$raw_log_with_biology -
      partial_confounding_quantized_matrix(
        reconstructed$raw_log_with_biology
      )
  ))
  clean_log_residual <- max(abs(
    reconstructed$clean_log_with_biology -
      partial_confounding_quantized_matrix(
        reconstructed$clean_log_with_biology
      )
  ))
  if (!is.finite(raw_log_residual) || !is.finite(clean_log_residual) ||
      raw_log_residual > portable_policy$quantization_unit ||
      clean_log_residual > portable_policy$quantization_unit) {
    stop("Unexpected log-matrix decimal-quantization residual for ",
         scenario$scenario_id[[1L]], ".", call. = FALSE)
  }
  execution_config_sha256 <- canonical_sha256_object(list(
    schema = "partial_confounding_portable_execution_config_v3",
    scenario_config_sha256 = scenario$scenario_config_sha256[[1L]],
    bundle_sha256 = scenario$bundle_sha256[[1L]],
    portable_hash_schema = portable_policy$schema,
    raw_log_with_biology_quantized_9_sha256 = raw_log_quantized,
    clean_log_with_biology_quantized_9_sha256 = clean_log_quantized,
    raw_log_with_biology_reference_object_sha256 =
      raw_log_reference_object_sha256,
    clean_log_with_biology_reference_object_sha256 =
      clean_log_reference_object_sha256,
    frozen_method_parameters_sha256 = expected_parameter_hash
  ))
  portable_rows[[index]] <- data.frame(
    array_index = index,
    seed_id = seed_id,
    scenario_id = scenario$scenario_id[[1L]],
    bundle_sha256 = scenario$bundle_sha256[[1L]],
    hash_schema = portable_policy$schema,
    quantization_digits = portable_policy$digits,
    equivalence_tolerance = portable_policy$quantization_unit,
    portable_execution_config_sha256_v3 = execution_config_sha256,
    raw_log_with_biology_quantized_9_sha256 = raw_log_quantized,
    clean_log_with_biology_quantized_9_sha256 = clean_log_quantized,
    raw_log_with_biology_reference_object_sha256 =
      raw_log_reference_object_sha256,
    clean_log_with_biology_reference_object_sha256 =
      clean_log_reference_object_sha256,
    raw_intensity_quantized_9_sha256_diagnostic =
      raw_intensity_quantized_diagnostic,
    clean_ground_truth_quantized_9_sha256_diagnostic =
      clean_intensity_quantized_diagnostic,
    raw_intensity_reference_object_sha256 =
      scenario$raw_intensity_sha256[[1L]],
    clean_ground_truth_reference_object_sha256 =
      scenario$clean_ground_truth_sha256[[1L]],
    stringsAsFactors = FALSE
  )
}
portable_hashes <- do.call(rbind, portable_rows)
if (nrow(portable_hashes) != 320L ||
    anyDuplicated(portable_hashes$scenario_id) ||
    anyDuplicated(portable_hashes$portable_execution_config_sha256_v3) ||
    anyNA(portable_hashes$raw_log_with_biology_quantized_9_sha256) ||
    anyNA(portable_hashes$clean_log_with_biology_quantized_9_sha256) ||
    any(!grepl("^[0-9a-f]{64}$", c(
      portable_hashes$raw_log_with_biology_reference_object_sha256,
      portable_hashes$clean_log_with_biology_reference_object_sha256
    )))) {
  stop("Portable matrix-hash freeze is incomplete.", call. = FALSE)
}

scenario_order <- scenario_manifest[, c(
  "seed_id", "scenario_id", "scenario_index", "confounding_axis",
  "nominal_strength", "cramers_v_phenotype_batch",
  "abs_cor_phenotype_within_batch_order",
  "weighted_mean_abs_within_plate_cor",
  "max_abs_within_plate_cor", "design_rank", "design_ncol",
  "design_full_rank", "design_condition_number",
  "scenario_config_sha256", "clean_ground_truth_sha256",
  "raw_intensity_sha256", "bundle_relative_path", "bundle_sha256"
), drop = FALSE]
scenario_order$array_index <- seq_len(nrow(scenario_order))
portable_index <- match(
  scenario_order$scenario_id, portable_hashes$scenario_id
)
scenario_order$portable_execution_config_sha256_v3 <-
  portable_hashes$portable_execution_config_sha256_v3[portable_index]
scenario_order$raw_log_with_biology_quantized_9_sha256 <-
  portable_hashes$raw_log_with_biology_quantized_9_sha256[portable_index]
scenario_order$clean_log_with_biology_quantized_9_sha256 <-
  portable_hashes$clean_log_with_biology_quantized_9_sha256[portable_index]
scenario_order$raw_log_with_biology_reference_object_sha256 <-
  portable_hashes$raw_log_with_biology_reference_object_sha256[portable_index]
scenario_order$clean_log_with_biology_reference_object_sha256 <-
  portable_hashes$clean_log_with_biology_reference_object_sha256[portable_index]
scenario_order$raw_intensity_quantized_9_sha256_diagnostic <-
  portable_hashes$raw_intensity_quantized_9_sha256_diagnostic[portable_index]
scenario_order$clean_ground_truth_quantized_9_sha256_diagnostic <-
  portable_hashes$clean_ground_truth_quantized_9_sha256_diagnostic[portable_index]
scenario_order <- scenario_order[, c(
  "array_index", setdiff(names(scenario_order), "array_index")
), drop = FALSE]

# Immutable accounting for the held/cancelled pre-v2 attempt. These rows are
# provenance only and are never eligible for scientific aggregation.
pre_v2_attempt_ledger <- data.frame(
  array_index = scenario_order$array_index,
  scenario_id = scenario_order$scenario_id,
  pre_v2_job_id = ifelse(scenario_order$array_index <= 49L,
                         "8385353", ""),
  pre_v2_attempt_status = ifelse(
    scenario_order$array_index <= 32L, "completed_v1_scientifically_ineligible",
    ifelse(
      scenario_order$array_index <= 48L,
      "failed_v1_intensity_round10_portability_gate",
      ifelse(
        scenario_order$array_index == 49L,
        "canonical_directory_observed_status_requires_archive",
        "not_started_pre_v2"
      )
    )
  ),
  pre_v2_scientifically_eligible = FALSE,
  v2_retry_required = TRUE,
  pre_v2_expected_raw_intensity_round10_sha256 =
    v2_raw_intensity_round10,
  pre_v2_expected_clean_intensity_round10_sha256 =
    v2_clean_intensity_round10,
  pre_v2_observed_raw_intensity_round10_sha256 = "",
  pre_v2_observed_clean_intensity_round10_sha256 = "",
  pre_v2_reference_raw_intensity_object_sha256 =
    scenario_order$raw_intensity_sha256,
  pre_v2_observed_raw_intensity_object_sha256 = "",
  pre_v2_reference_clean_ground_truth_object_sha256 =
    scenario_order$clean_ground_truth_sha256,
  pre_v2_observed_clean_ground_truth_object_sha256 = "",
  raw_intensity_max_abs_difference = NA_real_,
  raw_intensity_median_abs_difference = NA_real_,
  raw_intensity_nonzero_difference_cells = NA_integer_,
  raw_intensity_cells_above_1e_minus_12 = NA_integer_,
  raw_log_with_biology_max_abs_difference = NA_real_,
  raw_log_round10_identical = NA,
  stringsAsFactors = FALSE
)
conf03_evidence <- pre_v2_attempt_ledger$scenario_id == "CONF03_none_0p00"
pre_v2_attempt_ledger$raw_intensity_max_abs_difference[conf03_evidence] <-
  3.5527136788005009e-14
pre_v2_attempt_ledger$raw_intensity_median_abs_difference[conf03_evidence] <- 0
pre_v2_attempt_ledger$raw_intensity_nonzero_difference_cells[conf03_evidence] <-
  49116L
pre_v2_attempt_ledger$raw_intensity_cells_above_1e_minus_12[conf03_evidence] <- 0L
pre_v2_attempt_ledger$raw_log_with_biology_max_abs_difference[conf03_evidence] <-
  8.8817841970012523e-16
pre_v2_attempt_ledger$raw_log_round10_identical[conf03_evidence] <- TRUE
pre_v2_attempt_ledger$pre_v2_observed_raw_intensity_round10_sha256[
  conf03_evidence
] <- "0c982504c0a5d613e725e928088210124195914311ed61db19c269c359bfb8ca"
pre_v2_attempt_ledger$pre_v2_observed_clean_intensity_round10_sha256[
  conf03_evidence
] <- "8c43ea77bd96fd13531a20143d80ee569f69f3ab2c1ca33b4e1381e26deb0fb7"
pre_v2_attempt_ledger$pre_v2_observed_raw_intensity_object_sha256[
  conf03_evidence
] <- "64d96f7efb42abb8e6802cfa54bbb600b7e8d5598eed3319d4870f85da191e2f"
pre_v2_attempt_ledger$pre_v2_observed_clean_ground_truth_object_sha256[
  conf03_evidence
] <- "885f6373160236d3810de7e5903100c0852bc486e519900942b5b291976143a4"

# Immutable accounting for the canceled v2 round-10 log-hash attempt. Task 1
# was the independently completed gate; tasks 2--320 were array 8396751.
# None is reusable because v3 changes the authoritative portability contract.
v2_failed_indices <- c(129:144, 252L)
v2_cancelled_indices <- c(
  283L, 299L, 302L, 304:320
)
v2_completed_indices <- setdiff(seq_len(320L), c(
  v2_failed_indices, v2_cancelled_indices
))
stopifnot(
  length(v2_completed_indices) == 283L,
  length(v2_failed_indices) == 17L,
  length(v2_cancelled_indices) == 20L,
  setequal(c(v2_completed_indices, v2_failed_indices,
             v2_cancelled_indices), seq_len(320L))
)
pre_v3_attempt_ledger <- data.frame(
  array_index = scenario_order$array_index,
  scenario_id = scenario_order$scenario_id,
  v2_job_id = ifelse(scenario_order$array_index == 1L,
                     "8395516", "8396751"),
  v2_scheduler_status = ifelse(
    scenario_order$array_index %in% v2_failed_indices,
    "failed_raw_log_round10_portability_gate",
    ifelse(
      scenario_order$array_index %in% v2_cancelled_indices,
      "cancelled_after_v2_portability_failure",
      "completed_v2_scientifically_ineligible"
    )
  ),
  v2_scientifically_eligible = FALSE,
  v3_retry_required = TRUE,
  v2_expected_raw_log_round10_sha256 = v2_raw_log_round10,
  v2_expected_clean_log_round10_sha256 = v2_clean_log_round10,
  v3_expected_raw_log_round9_sha256 =
    scenario_order$raw_log_with_biology_quantized_9_sha256,
  v3_expected_clean_log_round9_sha256 =
    scenario_order$clean_log_with_biology_quantized_9_sha256,
  v3_local_vs_euler_login_round9_match = TRUE,
  v2_failure_archive = paste0(
    "full_grid/archive/",
    "partial_confounding_v2_round10_portability_failure_",
    "20260724T1000+0200"
  ),
  stringsAsFactors = FALSE
)

masters <- seed_ledger[seed_ledger$component == "master_seed",
                       c("seed_id", "seed"), drop = FALSE]
masters <- masters[match(expected_seed_ids, masters$seed_id), , drop = FALSE]
if (anyNA(masters$seed) || anyDuplicated(masters$seed_id)) {
  stop("Master-seed ledger is incomplete.", call. = FALSE)
}
hidden_seeds <- data.frame(
  seed_id = masters$seed_id,
  master_seed = as.integer(masters$seed),
  hidden_reference_seed = as.integer(masters$seed + 3000L),
  hidden_per_plate = 2L,
  expected_hidden_total = 10L,
  expected_training_total = 40L,
  derivation = "master_seed + 3000; frozen before full-grid execution",
  stringsAsFactors = FALSE
)
if (hidden_seeds$hidden_reference_seed[[1L]] != 202612001L) {
  stop("CONF01 hidden-reference seed no longer matches the validated pilot.",
       call. = FALSE)
}

metric_dictionary <- partial_full_metric_dictionary()
study_config <- list(
  schema = "partial_confounding_full_grid_config_v3",
  frozen_before_full_grid = TRUE,
  scenario_count = 320L,
  seed_count = 20L,
  scenarios_per_seed = 16L,
  analysis_unit_count_per_scenario = 18L,
  metric_count_per_method = nrow(metric_dictionary),
  expected_dimensions = c(features = 1000L, injections = 500L),
  expected_design = c(study = 450L, qc = 50L, plates = 5L),
  hidden_reference_policy = paste(
    "Two of ten QCs per plate; seed = master seed + 3000; hidden IDs",
    "relabeled ordinary Sample and excluded from all QC fitting/tuning."
  ),
  phenotype_blinding = paste(
    "No phenotype, response truth, or confounding labels are passed to any",
    "correction method or method cache context."
  ),
  primary_run_order_realization =
    "weighted_mean_abs_within_plate_cor",
  pooled_run_order_realization_status =
    "diagnostic_only_due_to_within_plate_sign_cancellation",
  nominal_one_label = paste(
    "maximal feasible under balanced five-plate design; not complete global",
    "confounding"
  ),
  effect_model =
    "log1p intensity ~ phenotype + plate + within_plate_position_scaled",
  non_estimability = paste(
    "Rank-deficient or aliased phenotype coefficients are NA with explicit",
    "not_estimable status; preprocessing is not claimed to solve the design."
  ),
  primary_effect_reference = "clean-matrix log1p estimand and paired zero",
  injected_effect_reference = "secondary latent-log diagnostic",
  full_gam_required = TRUE,
  full_ljung_box_required = TRUE,
  portable_reconstruction_hash = portable_policy,
  observed_cross_platform_diagnostic = list(
    comparison = paste(
      "all 320 scenarios, macOS arm64 versus R 4.5.1 Linux x86_64",
      "on the Euler Intel login node; failed v2 tasks also sampled AMD",
      "EPYC 7763 and EPYC 9654 compute nodes"
    ),
    raw_log_round10_mismatch_scenarios = 16L,
    clean_log_round10_mismatch_scenarios = 0L,
    raw_log_round9_mismatch_scenarios = 0L,
    clean_log_round9_mismatch_scenarios = 0L,
    v2_array_failures = 17L,
    v2_array_failure_detail = paste(
      "CONF09 all 16 scenarios plus CONF16_run_order_0p90 stopped before",
      "method execution at the raw-log round-10 preflight gate"
    ),
    v3_status = paste(
      "round-9 is the narrowest tested decimal gate with complete",
      "cross-platform identity; dimensions and dimnames remain hashed"
    )
  ),
  bootstrap_replicates = 2000L,
  bootstrap_seed = 202613001L,
  phenotype_aware_combat = paste(
    "not run: unchanged run_combat() wrapper has no phenotype-covariate",
    "interface; never mix with phenotype-blind ranking"
  )
)
degradation <- list(
  schema = "partial_confounding_material_degradation_v1",
  frozen_before_full_grid = TRUE,
  primary_reference = "clean-observed log1p effect and paired zero scenario",
  criteria = list(
    clean_attenuation =
      "median_attenuation_ratio_clean_responsive < 0.80",
    clean_effect_rmse = "effect_rmse_vs_clean_responsive > 0.20",
    paired_tpr_loss = "TPR_zero - TPR_scenario > 0.20"
  ),
  material_degradation = "any criterion is TRUE",
  non_estimable =
    "NA degradation flag; retain explicit design/effect status",
  full_curves_primary = TRUE
)

write_csv(scenario_order, file.path(config_dir, "scenario_order.csv"))
write_csv(portable_hashes,
          file.path(config_dir, "portable_matrix_hashes.csv"))
write_csv(pre_v2_attempt_ledger,
          file.path(config_dir, "pre_v2_attempt_ledger.csv"))
write_csv(pre_v3_attempt_ledger,
          file.path(config_dir, "pre_v3_attempt_ledger.csv"))
write_csv(hidden_seeds, file.path(config_dir, "hidden_reference_seeds.csv"))
write_csv(method_parameters, file.path(config_dir, "method_parameters.csv"))
writeLines(expected_parameter_hash,
           file.path(config_dir, "method_parameters.sha256"))
write_csv(metric_dictionary, file.path(config_dir, "metric_dictionary.csv"))
jsonlite::write_json(
  study_config, file.path(config_dir, "study_config.json"),
  auto_unbox = TRUE, pretty = TRUE
)
jsonlite::write_json(
  degradation, file.path(config_dir, "material_degradation.json"),
  auto_unbox = TRUE, pretty = TRUE
)

readme <- c(
  "# Frozen partial-confounding full-grid configuration",
  "",
  "This directory was generated before any full-grid correction run.",
  "It fixes the 320 scenario order, 18 analysis labels, 37 metrics,",
  "hidden-reference seeds, and the material-degradation rule.",
  "`portable_matrix_hashes.csv` freezes v3 SHA-256 identities after rounding",
  "the complete pre-exponentiation raw/clean log matrices to 9 decimal",
  "places. Exact bundle file hashes remain mandatory. Intensity object and",
  "intensity-round9 hashes are platform diagnostics, not identity gates.",
  "Exact reference-platform raw/clean log-object hashes are also frozen as",
  "diagnostics; only the quantized log hashes are cross-platform gates.",
  "The equivalence tolerance is one 1e-9 quantization cell; exact dimensions",
  "and dimnames are also included in the serialized digest.",
  "Every scenario has a new portable_execution_config_sha256_v3; v1 runs",
  "and v2 runs must be archived and recomputed and are never silently reused.",
  "",
  "Nominal strength 1.0 means maximal feasible confounding under the",
  "balanced five-plate design, not complete global confounding. The primary",
  "realized run-order axis is the weighted mean absolute within-plate",
  "phenotype/order correlation; pooled correlation is diagnostic only.",
  "",
  "Phenotype-aware ComBat is unavailable through the unchanged standard",
  "wrapper and is not included in the phenotype-blind ranking."
)
writeLines(readme, file.path(config_dir, "README.md"))

writeLines(c(
  "# Partial-confounding v2 archive and rerun protocol",
  "",
  "Array 8385353 was held/cancelled after the v1 intensity-round10 gate",
  "failed. The observed remote pre-v2 inventory was 33 canonical run",
  "directories, 49 preflight directories, and 96 scheduler log files.",
  "Before syncing v2, archive the entire runs/, validation/preflight/ and",
  "validation/dry_run/ trees, the complete v1 config directory, and all",
  "job-8385353 plus prior gate logs. Record archive paths and SHA-256 values.",
  "Never delete or selectively overwrite failed/incomplete directories.",
  "",
  "`pre_v2_attempt_ledger.csv` retains the task-level v1 attempt state and",
  "the exact CONF03 expected/observed hashes and numerical comparison.",
  "Every scenario, including the old task 1 result, requires a v2 run.",
  "",
  "Release order:",
  "1. Validate the regenerated config manifest and unchanged 320 bundle hashes.",
  "2. On Euler, run CONF03_none_0p00 with --dry-run --force and retain proof",
  "   that both authoritative log-round10 hashes pass while the raw-intensity",
  "   round10 diagnostic mismatch remains non-blocking.",
  "3. Run a full CONF01_none_0p00 v2 task with --force and all native checks.",
  "4. Only after both gates pass, submit array tasks 2-320 with the v2 rerun",
  "   flag so every pre-v2 directory is archived before recomputation.",
  "5. Launch the strict v2 aggregate only with an afterok dependency on that",
  "   array; v1 completion schemas are scientifically ineligible."
), file.path(config_dir, "PRE_V2_ARCHIVE_AND_RERUN_PROTOCOL.md"))

writeLines(c(
  "# Partial-confounding v3 archive and rerun protocol",
  "",
  "V2 gate job 8395516 completed task 1. Array 8396751 then produced",
  "282 completed tasks, 17 raw-log round-10 portability-gate failures, and",
  "20 tasks canceled when the defect was diagnosed. No correction method ran",
  "inside the 17 failed tasks. A source-identical EPYC 9654 retry reproduced",
  "the CONF09 failure.",
  "",
  "An independent all-320 reconstruction comparison found 16 raw-log hash",
  "mismatches at round-10, zero clean-log mismatches at round-10, and zero",
  "raw or clean mismatches at round-9 and round-8. Round-9 is therefore the",
  "narrowest tested portable decimal gate. Exact bundle hashes, dimensions,",
  "dimnames, source-subset identity, and all scientific parameters remain",
  "mandatory and unchanged.",
  "",
  "Before syncing v3, move the complete v2 runs/, validation/, config/, and",
  "scheduler logs into the checksum-manifested archive named in",
  "`pre_v3_attempt_ledger.csv`. Never selectively reuse a v2 completion.",
  "",
  "Release order:",
  "1. Validate the v3 config manifest and unchanged 320 bundle hashes.",
  "2. On an AMD compute node, run CONF09_none_0p00 with --dry-run --force",
  "   and retain proof that both authoritative log-round9 hashes pass.",
  "3. Run a full CONF01_none_0p00 v3 gate with --force and all native checks.",
  "4. Submit tasks 2--320 with the v3 rerun flag only after both gates pass.",
  "5. Launch the strict v3 aggregate with an afterok dependency on the array.",
  "   V1 and v2 completion schemas are scientifically ineligible."
), file.path(config_dir, "PRE_V3_ARCHIVE_AND_RERUN_PROTOCOL.md"))

config_files <- list.files(
  config_dir, full.names = TRUE, recursive = FALSE
)
config_files <- config_files[basename(config_files) != "config_manifest.csv"]
config_manifest <- data.frame(
  relative_path = basename(config_files),
  bytes = as.numeric(file.info(config_files)$size),
  sha256 = vapply(config_files, canonical_sha256_file, character(1)),
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)
write_csv(config_manifest, file.path(config_dir, "config_manifest.csv"))
message("Frozen full-grid configuration prepared: 320 scenarios x 18 labels x ",
        nrow(metric_dictionary), " metrics.")
