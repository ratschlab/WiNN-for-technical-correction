#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
for (package in c("digest", "jsonlite", "matrixStats", "irlba")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package, call. = FALSE)
}
source(file.path(repo_root, "analysis", "derivative_biology_helpers.R"), local = TRUE)
source(file.path(repo_root, "analysis", "ablation_stage_helpers.R"), local = TRUE)

root <- file.path(repo_root, "results", "final", "derivative")
datasets <- derivative_dataset_registry()
methods <- derivative_method_registry()
config_path <- file.path(repo_root, "analysis", "config", "derivative_biology_config.json")
config_sha <- digest::digest(config_path, file = TRUE, algo = "sha256")

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required derivative output: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

for (dataset_key in datasets$dataset_key) {
  manifest <- read_required(file.path(root, dataset_key, "derivative_run_manifest.csv"))
  stopifnot(
    nrow(manifest) == 1L, !manifest$test_mode,
    manifest$bootstrap_replicates == 2000L,
    manifest$config_sha256 == config_sha,
    !manifest$correction_methods_executed, !manifest$tuning_executed, !manifest$selection_executed
  )
  estimates <- read_required(file.path(root, dataset_key, "bootstrap_estimates.csv"))
  stopifnot(!anyDuplicated(estimates[c("method_id", "metric")]))
  stopifnot(setequal(unique(estimates$method_id), methods$method_id))
  replicates <- read_required(file.path(root, dataset_key, "bootstrap_replicates.csv"))
  stopifnot(length(unique(replicates$bootstrap)) == 2000L)
  if (dataset_key %in% c("mtbls79", "sacurine", "waveica")) {
    permutation <- read_required(file.path(root, dataset_key, "permutation_null.csv"))
    stopifnot(all(table(interaction(permutation$method_id, permutation$endpoint, drop = TRUE)) == 5000L))
    lobo <- read_required(file.path(root, dataset_key, "leave_one_batch_out_bootstrap_summary.csv"))
    stopifnot(nrow(lobo) > 0L, all(lobo$valid_bootstrap > 0L))
  }
}

# Reconstruct and validate every cumulative ablation trajectory.
cumulative_paths <- list.files(
  file.path(repo_root, "results", "final", "tables"),
  pattern = "^ablation_cumulative_values_.*[.]csv$", full.names = TRUE
)
cumulative <- do.call(rbind, lapply(cumulative_paths, read.csv, stringsAsFactors = FALSE, check.names = FALSE))
impacts <- read_required(file.path(repo_root, "results", "final", "tables", "ablation_stage_impacts.csv"))
validate_ablation_stage_impacts(cumulative, impacts)
stopifnot(!anyDuplicated(impacts[c("dataset_key", "table_id", "metric", "stage")]))
stopifnot(all(!is.na(impacts$metric) & nzchar(impacts$metric)))

# Categorical targets reject numeric encodings, and saved empirical batch values
# are reproduced with explicit factors from the same common frozen panel.
numeric_guard <- try(
  weighted_pc_r2_frozen(matrix(rexp(120), nrow = 10), rep(1:3, each = 4), "categorical"),
  silent = TRUE
)
stopifnot(inherits(numeric_guard, "try-error"))
pre <- read_required(file.path(repo_root, "results", "frozen_endpoint_free_inputs", "provenance", "matrix_checksums_pre.csv"))
for (dataset_key in setdiff(datasets$dataset_key, "simulation")) {
  frozen <- derivative_load_matrices(repo_root, dataset_key, pre)
  meta <- derivative_load_metadata(repo_root, dataset_key, frozen$sample_ids)
  study <- which(meta$is_study)
  observed <- weighted_pc_r2_frozen(
    frozen$matrices$raw[, study, drop = FALSE], factor(meta$batch[study]), "categorical"
  )
  point <- read_required(file.path(root, dataset_key, "biological_point_metrics.csv"))
  saved <- point$value[point$method_id == "raw" & point$metric == "batch_weighted_pc_r2_categorical"]
  stopifnot(length(saved) == 1L, isTRUE(all.equal(observed, saved, tolerance = 1e-12)))
}

# Separate Sacurine endpoints and neutral WaveICA labels are mandatory.
sacurine <- read_required(file.path(root, "sacurine", "sacurine_covariate_associations.csv"))
stopifnot(setequal(unique(sacurine$endpoint), c("age", "bmi", "gender")))
waveica <- read_required(file.path(root, "waveica", "association_features.csv"))
stopifnot(nrow(waveica) > 0L)
wave_meta <- read.csv(
  file.path(repo_root, "data", "public", "waveica_adenocarcinoma", "processed", "waveica_metadata.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
stopifnot(setequal(unique(wave_meta$biological_group[!as.logical(wave_meta$is_qc)]), c("group_0", "group_1")))

# Figure inputs exist and clinical public derivatives contain no direct cohort IDs.
figure_manifest <- read_required(file.path(root, "figures", "figure_manifest.csv"))
stopifnot(all(file.exists(file.path(root, figure_manifest$source_csv))))
stopifnot(all(file.exists(file.path(root, figure_manifest$pdf))))
stopifnot(all(file.exists(file.path(root, figure_manifest$png))))
clinical_files <- list.files(file.path(root, "clinical_fiams"), recursive = TRUE, full.names = TRUE)
clinical_text <- paste(unlist(lapply(clinical_files[grepl("[.]csv$", clinical_files)], readLines, warn = FALSE)), collapse = "\n")
stopifnot(!grepl("HG[0-9]{4}[.]", clinical_text))

# Every matrix declared for this derivative run has identical pre/post hashes.
comparison <- read_required(file.path(root, "provenance", "matrix_checksum_comparison.csv"))
stopifnot(nrow(comparison) == nrow(pre), nrow(comparison) >= 54L, all(comparison$unchanged))

non_estimable <- read_required(file.path(root, "tables", "non_estimable_analyses.csv"))
stopifnot(nrow(non_estimable) > 0L, all(nzchar(non_estimable$reason)))
analysis_manifest <- jsonlite::read_json(file.path(root, "derivative_analysis_manifest.json"), simplifyVector = TRUE)
stopifnot(
  isTRUE(analysis_manifest$frozen_matrices_unchanged),
  identical(analysis_manifest$weighted_pc_solver_equivalence_test, "analysis/test_derivative_pca_solver.R"),
  !analysis_manifest$correction_methods_executed,
  !analysis_manifest$tuning_executed,
  !analysis_manifest$selection_executed
)

# Manuscript-facing derivatives must be complete and traceable to saved results.
manuscript_root <- file.path(root, "manuscript")
required_manuscript <- c(
  "ablation_tables.tex", "full_benchmark_tables.csv", "full_benchmark_tables.tex",
  "main_simulation_rows.csv", "main_simulation_rows.tex",
  "supp_categorical_r2_trajectory.csv", "supp_categorical_r2_trajectory.tex",
  "supp_derivative_biology_rows.tex", "supp_reference_split_stability.csv",
  "supp_reference_split_stability.tex", "supp_partial_confounding.csv",
  "supp_partial_confounding.tex", "manuscript_summary_manifest.csv"
)
stopifnot(all(file.exists(file.path(manuscript_root, required_manuscript))))
manuscript_manifest <- read_required(file.path(manuscript_root, "manuscript_summary_manifest.csv"))
stopifnot(
  all(manuscript_manifest$generated_from_frozen_results),
  !any(manuscript_manifest$correction_methods_executed),
  !any(manuscript_manifest$tuning_executed),
  !any(manuscript_manifest$selection_executed),
  all(file.exists(file.path(manuscript_root, manuscript_manifest$artifact)))
)
reference_summary <- read_required(file.path(manuscript_root, "supp_reference_split_stability.csv"))
stopifnot(nrow(reference_summary) == 6L, all(reference_summary$splits == 10L), all(is.finite(reference_summary$mean_cv_change)))
partial_summary <- read_required(file.path(manuscript_root, "supp_partial_confounding.csv"))
maximal_joint <- partial_summary[
  partial_summary$confounding_type == "joint" & partial_summary$level == 1,
]
source_maximal <- read_required(file.path(repo_root, "results", "final", "tables", "partial_confounding_summary.csv"))
source_maximal <- source_maximal[
  source_maximal$method == "WiNN fixed default (no QC)" &
    source_maximal$confounding_type == "joint" & source_maximal$level == 1 &
    source_maximal$metric %in% c("median_attenuation_ratio_injected_responsive", "true_positive_rate"),
]
stopifnot(
  nrow(maximal_joint) == 1L,
  isTRUE(all.equal(
    maximal_joint$mean.median_attenuation_ratio_injected_responsive,
    source_maximal$mean[source_maximal$metric == "median_attenuation_ratio_injected_responsive"],
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    maximal_joint$mean.true_positive_rate,
    source_maximal$mean[source_maximal$metric == "true_positive_rate"],
    tolerance = 1e-12
  ))
)

cat("Derivative biological-preservation, provenance, and ablation tests passed.\n")
