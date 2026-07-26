#!/usr/bin/env Rscript

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root) || !dir.exists(release_root)) {
  stop("WINN_RELEASE_ROOT must identify the release directory.", call. = FALSE)
}

required_packages <- c("dplyr", "tidyr", "ggplot2", "jsonlite", "digest", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

dataset_order <- c(
  "simulation", "mtbls79", "clinical_fiams", "batchcorr_set1",
  "sacurine", "waveica"
)
dataset_labels <- c(
  simulation = "Simulation", mtbls79 = "MTBLS79",
  clinical_fiams = "Clinical FIA-MS", batchcorr_set1 = "BatchCorr Set 1",
  sacurine = "Sacurine", waveica = "WaveICA"
)
method_order <- c(
  "Raw", "ComBat", "QC-RLSC", "QC-RFSC", "TIGER", "SERRF",
  "WiNN auto (QC)", "WiNN auto-batch (QC)",
  "WiNN fixed default (no QC)"
)
method_palette <- c(
  "Raw" = "#5F6368", "ComBat" = "#0072B2", "QC-RLSC" = "#009E73",
  "QC-RFSC" = "#56B4E9", "TIGER" = "#E69F00", "SERRF" = "#D55E00",
  "WiNN auto (QC)" = "#CC79A7", "WiNN auto-batch (QC)" = "#6F4E7C",
  "WiNN fixed default (no QC)" = "#B79F00"
)
variant_order <- c(
  "C0_RAW", "C1_OUTLIER", "C2_SELECTIVE_DRIFT",
  "C3_SELECTIVE_BATCH", "C4_FULL_FIXED"
)
variant_labels <- c(
  C0_RAW = "Raw", C1_OUTLIER = "Outlier shrinkage",
  C2_SELECTIVE_DRIFT = "Selective drift",
  C3_SELECTIVE_BATCH = "Selective batch", C4_FULL_FIXED = "PQN"
)

normalise_method <- function(value) {
  value <- sub("^WINN default \\(no QC\\)$", "WiNN fixed default (no QC)", value)
  value <- sub("^WiNN default \\(no QC\\)$", "WiNN fixed default (no QC)", value)
  value
}

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

json_scalar <- function(value, name, default = NA_character_) {
  field <- value[[name]]
  if (is.null(field) || !length(field)) default else as.character(field[[1L]])
}

complete_count <- function(path) {
  if (!dir.exists(path)) return(0L)
  length(list.files(path, pattern = "^complete[.]json$", recursive = TRUE, full.names = TRUE))
}

expected <- c(
  primary_evaluations = 6L,
  ablation_evaluations = 6L,
  simulation_seed_evaluations = 10L,
  reference_split_evaluations = 60L,
  partial_confounding_tasks = 160L
)
observed <- c(
  primary_evaluations = complete_count(file.path(release_root, "results", "evaluation", "primary")),
  ablation_evaluations = complete_count(file.path(release_root, "results", "evaluation", "ablations")),
  simulation_seed_evaluations = complete_count(file.path(release_root, "results", "evaluation", "simulation_seeds")),
  reference_split_evaluations = complete_count(file.path(release_root, "results", "evaluation", "reference_splits")),
  partial_confounding_tasks = complete_count(file.path(release_root, "results", "partial_confounding"))
)
completion <- data.frame(
  analysis_family = names(expected), expected = as.integer(expected),
  observed = as.integer(observed), complete = observed == expected,
  stringsAsFactors = FALSE
)
if (!all(completion$complete)) {
  detail <- paste(
    completion$analysis_family[!completion$complete],
    paste0(completion$observed[!completion$complete], "/", completion$expected[!completion$complete]),
    collapse = "; "
  )
  stop("Release aggregation refused because task families are incomplete: ", detail, call. = FALSE)
}

final_dir <- file.path(release_root, "results", "final")
table_dir <- file.path(final_dir, "tables")
figure_dir <- file.path(final_dir, "figures")
source_dir <- file.path(figure_dir, "source_data")
validation_dir <- file.path(final_dir, "validation")
for (path in c(final_dir, table_dir, figure_dir, source_dir, validation_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
write.csv(completion, file.path(validation_dir, "analysis_family_completion.csv"), row.names = FALSE)

read_family <- function(family, filename) {
  values <- lapply(dataset_order, function(dataset_key) {
    path <- file.path(release_root, "results", "evaluation", family, dataset_key, filename)
    value <- read_required(path)
    if (!"dataset_key" %in% names(value)) value$dataset_key <- dataset_key
    value
  })
  bind_rows(values)
}

primary <- read_family("primary", "method_metrics.csv")
primary$method <- normalise_method(primary$method)
primary$dataset <- unname(dataset_labels[primary$dataset_key])
primary$dataset <- factor(primary$dataset, levels = unname(dataset_labels))
primary$method <- factor(primary$method, levels = method_order)
primary <- primary[order(primary$dataset, primary$method, primary$metric), ]

primary_tasks <- read_required(file.path(
  release_root, "analysis", "config", "primary_run_matrix.csv"
))
primary_runtime_policy <- bind_rows(lapply(seq_len(nrow(primary_tasks)), function(index) {
  task <- primary_tasks[index, , drop = FALSE]
  manifest <- jsonlite::read_json(
    file.path(release_root, as.character(task$output_dir), "run_manifest.json"),
    simplifyVector = TRUE
  )
  reuse_decision <- json_scalar(manifest, "reuse_decision")
  execution_runtime <- as.numeric(json_scalar(manifest, "runtime_sec"))
  reportable <- as.character(task$method_id) == "raw" || reuse_decision != "reuse_exact"
  data.frame(
    analysis_family = "primary", task_id = as.character(task$task_id),
    dataset_key = as.character(task$dataset_key), seed_id = NA_character_,
    method_id = as.character(task$method_id), method = as.character(task$method),
    reuse_decision = reuse_decision,
    execution_runtime_sec = execution_runtime,
    reported_method_runtime_sec = if (as.character(task$method_id) == "raw") 0 else if (reportable) execution_runtime else NA_real_,
    runtime_provenance = if (as.character(task$method_id) == "raw") {
      "no correction"
    } else if (reportable) {
      "measured method execution in frozen release"
    } else {
      "not remeasured; cache retrieval time excluded"
    },
    stringsAsFactors = FALSE
  )
}))
runtime_rows <- primary$metric == "runtime_sec"
policy_index <- match(
  paste(primary$dataset_key, primary$method_id),
  paste(primary_runtime_policy$dataset_key, primary_runtime_policy$method_id)
)
primary$value[runtime_rows] <- primary_runtime_policy$reported_method_runtime_sec[policy_index[runtime_rows]]
primary$notes[runtime_rows] <- primary_runtime_policy$runtime_provenance[policy_index[runtime_rows]]
write.csv(primary, file.path(table_dir, "primary_method_metrics_long.csv"), row.names = FALSE, na = "")

primary_wide <- primary |>
  select(dataset_key, dataset, method_id, method, metric, value) |>
  pivot_wider(names_from = metric, values_from = value) |>
  arrange(dataset, method)
write.csv(primary_wide, file.path(table_dir, "primary_method_metrics_wide.csv"), row.names = FALSE, na = "")

by_dataset_dir <- file.path(table_dir, "by_dataset")
dir.create(by_dataset_dir, recursive = TRUE, showWarnings = FALSE)
for (dataset_key in dataset_order) {
  write.csv(
    primary_wide[primary_wide$dataset_key == dataset_key, , drop = FALSE],
    file.path(by_dataset_dir, paste0(dataset_key, "_method_metrics.csv")),
    row.names = FALSE, na = ""
  )
}

winn_modes <- c("winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc")
write.csv(
  primary[primary$method_id %in% winn_modes, , drop = FALSE],
  file.path(table_dir, "winn_mode_comparison.csv"), row.names = FALSE, na = ""
)

coverage_and_counts <- primary |>
  filter(metric %in% c(
    "retained_features", "retained_samples", "drift_profiles_corrected",
    "drift_features_corrected", "batch_features_corrected"
  )) |>
  select(dataset_key, dataset, method_id, method, metric, value, denominator, units)
write.csv(
  coverage_and_counts,
  file.path(table_dir, "feature_retention_and_correction_counts.csv"),
  row.names = FALSE, na = ""
)

selection_manifest <- read_required(file.path(
  release_root, "analysis", "config", "endpoint_free_selection_manifest.csv"
))
write.csv(
  selection_manifest,
  file.path(table_dir, "configuration_selection_manifest.csv"),
  row.names = FALSE, na = ""
)
primary_competitor_candidates <- read_required(file.path(
  release_root, "analysis", "config", "primary_competitor_candidate_scores.csv"
))
write.csv(
  primary_competitor_candidates,
  file.path(table_dir, "primary_competitor_candidate_scores.csv"),
  row.names = FALSE, na = ""
)
reuse_audit <- read.delim(
  file.path(release_root, "analysis", "config", "competitor_reuse_audit.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)

condition_counts <- function(path) {
  if (!file.exists(path)) return(c(warnings = NA_integer_, errors = NA_integer_))
  lines <- readLines(path, warn = FALSE)
  warning_start <- match("WARNINGS", lines)
  message_start <- match("MESSAGES", lines)
  error_start <- match("ERROR", lines)
  warning_lines <- if (is.finite(warning_start) && is.finite(message_start) &&
                       message_start > warning_start + 1L) {
    lines[seq.int(warning_start + 1L, message_start - 1L)]
  } else character()
  error_lines <- if (is.finite(error_start) && error_start < length(lines)) {
    lines[seq.int(error_start + 1L, length(lines))]
  } else character()
  c(
    warnings = sum(nzchar(trimws(warning_lines))),
    errors = sum(nzchar(trimws(error_lines)))
  )
}
primary_provenance <- bind_rows(lapply(seq_len(nrow(primary_tasks)), function(index) {
  task <- primary_tasks[index, , drop = FALSE]
  directory <- file.path(release_root, as.character(task$output_dir))
  manifest <- jsonlite::read_json(
    file.path(directory, "run_manifest.json"), simplifyVector = TRUE
  )
  counts <- condition_counts(file.path(directory, "conditions.log"))
  data.frame(
    task_id = as.character(task$task_id),
    dataset_key = as.character(task$dataset_key),
    method_id = as.character(task$method_id),
    method = as.character(task$method),
    reuse_decision = json_scalar(manifest, "reuse_decision"),
    seed = as.integer(json_scalar(manifest, "seed")),
    execution_runtime_sec = as.numeric(json_scalar(manifest, "runtime_sec")),
    reported_method_runtime_sec = primary_runtime_policy$reported_method_runtime_sec[
      match(as.character(task$task_id), primary_runtime_policy$task_id)
    ],
    runtime_provenance = primary_runtime_policy$runtime_provenance[
      match(as.character(task$task_id), primary_runtime_policy$task_id)
    ],
    wall_sec = as.numeric(json_scalar(manifest, "wall_sec")),
    n_features = as.integer(json_scalar(manifest, "n_features")),
    n_samples = as.integer(json_scalar(manifest, "n_samples")),
    warning_lines = as.integer(counts[["warnings"]]),
    error_lines = as.integer(counts[["errors"]]),
    input_object_sha256 = json_scalar(manifest, "input_object_sha256"),
    output_object_sha256 = json_scalar(manifest, "output_object_sha256"),
    matrix_file_sha256 = json_scalar(manifest, "matrix_file_sha256"),
    conditions_log_sha256 = if (file.exists(file.path(directory, "conditions.log"))) {
      digest::digest(file.path(directory, "conditions.log"), file = TRUE, algo = "sha256")
    } else NA_character_,
    package_version = json_scalar(manifest$package, "version"),
    package_source_sha256 = json_scalar(manifest$package, "source_archive_sha256"),
    started_at = json_scalar(manifest, "started_at"),
    completed_at = json_scalar(manifest, "completed_at"),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  primary_provenance, file.path(table_dir, "primary_run_provenance.csv"),
  row.names = FALSE, na = ""
)

primary_winn_parameters <- bind_rows(lapply(dataset_order, function(dataset_key) {
  bind_rows(lapply(winn_modes, function(method_id) {
    details <- readRDS(file.path(
      release_root, "results", "primary", dataset_key, method_id,
      "method_details.rds"
    ))
    selected <- details$selected_parameters
    scalar <- function(name, default = NA) {
      value <- selected[[name]]
      if (is.null(value) || !length(value)) default else value[[1L]]
    }
    data.frame(
      dataset_key = dataset_key,
      dataset = unname(dataset_labels[dataset_key]),
      method_id = method_id,
      method = method_order[match(method_id, winn_modes) + 6L],
      parameters = as.character(scalar("parameters", "fixed")),
      batch_option = as.character(scalar("batch_option", NA_character_)),
      pelt_penalty = as.numeric(scalar("pelt_penalty", NA_real_)),
      test = as.character(scalar("test", NA_character_)),
      ljung_box_fitdf = as.integer(scalar("ljung_box_fitdf", NA_integer_)),
      lag = as.character(scalar("lag", NA_character_)),
      spline_method = as.character(scalar("spline_method", NA_character_)),
      acorr_fdr = as.numeric(scalar("acorr_fdr", NA_real_)),
      anova_fdr = as.numeric(scalar("anova_fdr", NA_real_)),
      normalization = as.character(scalar("normalization", NA_character_)),
      scale_by_batch = as.logical(scalar("scale_by_batch", NA)),
      training_qc_quality_score = as.numeric(scalar("quality_score", NA_real_)),
      training_qc_mean_cv = as.numeric(scalar("mean_cv", NA_real_)),
      training_qc_mean_correlation = as.numeric(scalar("mean_correlation", NA_real_)),
      detected_batch_count = length(unique(details$batch)),
      batch_source = if (is.null(details$batch_source)) NA_character_ else as.character(details$batch_source),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(
  primary_winn_parameters,
  file.path(table_dir, "primary_winn_selected_parameters.csv"),
  row.names = FALSE, na = ""
)

reuse_audit_final <- reuse_audit |>
  left_join(
    primary_provenance |>
      select(
        dataset_key, method_id, recorded_reuse_decision = reuse_decision,
        seed, n_features, n_samples, input_object_sha256,
        output_object_sha256, matrix_file_sha256, package_version,
        package_source_sha256, warning_lines, error_lines
      ),
    by = c("dataset_key", "method_id")
  ) |>
  mutate(
    decision_matches_run_manifest =
      decision == recorded_reuse_decision |
      public_reproduction_decision == recorded_reuse_decision
  )
write.csv(
  reuse_audit_final, file.path(table_dir, "competitor_reuse_audit.csv"),
  row.names = FALSE, na = ""
)

headline_metrics <- c(
  "heldout_qc_cv_mean", "residual_run_order_gam_mean_deviance",
  "residual_ljung_box_proportion_significant",
  "batch_weighted_pc_r2_categorical", "biological_weighted_pc_r2",
  "truth_metabolite_profile_icc_mean", "truth_sample_profile_pearson_mean",
  "truth_log1p_rmse", "metabolite_icc_a1_median",
  "feature_repeatability_ratio_median",
  "sample_replicate_pearson_median", "genuine_replicate_icc_median",
  "cross_batch_effect_pearson_median", "biological_associated_features",
  "retained_features", "retained_samples", "runtime_sec"
)
combined_primary <- primary |>
  filter(metric %in% headline_metrics) |>
  select(dataset_key, dataset, method_id, method, metric, value, direction, denominator, units, notes)
write.csv(combined_primary, file.path(table_dir, "combined_primary_benchmark.csv"), row.names = FALSE, na = "")

runtime <- primary |>
  filter(metric %in% c("runtime_sec", "retained_features", "retained_samples")) |>
  select(dataset_key, dataset, method_id, method, metric, value, denominator, units)
write.csv(runtime, file.path(table_dir, "runtime_and_retention.csv"), row.names = FALSE, na = "")

ablations <- read_family("ablations", "method_metrics.csv")
ablations$dataset <- unname(dataset_labels[ablations$dataset_key])
write.csv(ablations, file.path(table_dir, "ablation_method_metrics_long.csv"), row.names = FALSE, na = "")

impact_definitions <- data.frame(
  table_id = c("heldout_qc_cv", "gam_deviance", "batch_r2", "metabolite_icc_a1", "sample_pearson", "feature_repeatability_ratio"),
  display_name = c("Held-out QC CV", "Residual run-order GAM deviance", "Batch weighted-PC R-squared", "Metabolite ICC(A,1)", "Sample Pearson", "Feature repeatability ratio"),
  simulation_metric = c("heldout_qc_cv_mean", "residual_run_order_gam_mean_deviance", "batch_weighted_pc_r2_categorical", "truth_metabolite_profile_icc_mean", "truth_sample_profile_pearson_mean", NA_character_),
  empirical_metric = c("heldout_qc_cv_mean", "residual_run_order_gam_mean_deviance", "batch_weighted_pc_r2_categorical", "metabolite_icc_a1_median", "sample_replicate_pearson_median", "feature_repeatability_ratio_median"),
  direction = c("lower", "lower", "lower", "higher", "higher", "higher"),
  stringsAsFactors = FALSE
)

stage_source <- list()
for (row_index in seq_len(nrow(impact_definitions))) {
  definition <- impact_definitions[row_index, ]
  rows <- lapply(dataset_order, function(dataset_key) {
    metric_name <- if (dataset_key == "simulation") definition$simulation_metric else definition$empirical_metric
    metric_unavailable <-
      (definition$table_id == "metabolite_icc_a1" &&
         !dataset_key %in% c("simulation", "mtbls79")) ||
      (definition$table_id == "feature_repeatability_ratio" &&
         !dataset_key %in% c("clinical_fiams", "batchcorr_set1")) ||
      (definition$table_id == "sample_pearson" &&
         dataset_key %in% c("sacurine", "waveica"))
    if (metric_unavailable) {
      values <- rep(NA_real_, length(variant_order))
    } else {
      selected <- ablations[
        ablations$dataset_key == dataset_key & ablations$metric == metric_name,
        c("method_id", "value"), drop = FALSE
      ]
      values <- selected$value[match(variant_order, selected$method_id)]
    }
    data.frame(
      dataset_key = dataset_key, dataset = unname(dataset_labels[dataset_key]),
      table_id = definition$table_id, metric = metric_name,
      direction = definition$direction,
      Raw = values[1], Outlier_shrinkage = values[2],
      Selective_drift = values[3], Selective_batch = values[4], PQN = values[5],
      stringsAsFactors = FALSE
    )
  }) |>
    bind_rows()
  value_path <- file.path(table_dir, paste0("ablation_cumulative_values_", definition$table_id, ".csv"))
  write.csv(rows, value_path, row.names = FALSE, na = "")

  sign <- if (definition$direction == "lower") -1 else 1
  impacts <- rows |>
    transmute(
      dataset_key, dataset, table_id, metric, direction,
      Outlier_shrinkage = sign * (Outlier_shrinkage - Raw),
      Selective_drift = sign * (Selective_drift - Outlier_shrinkage),
      Selective_batch = sign * (Selective_batch - Selective_drift),
      PQN = sign * (PQN - Selective_batch),
      impact_definition = "Positive values favor the newly added stage; negative values indicate deterioration."
    )
  write.csv(
    impacts,
    file.path(table_dir, paste0("ablation_step_impacts_", definition$table_id, ".csv")),
    row.names = FALSE, na = ""
  )
  stage_source[[definition$table_id]] <- impacts |>
    pivot_longer(
      cols = c(Outlier_shrinkage, Selective_drift, Selective_batch, PQN),
      names_to = "stage", values_to = "favorable_impact"
    )
}
stage_source <- bind_rows(stage_source)
write.csv(stage_source, file.path(source_dir, "ablation_step_impacts.csv"), row.names = FALSE, na = "")

gate_ids <- c("G_SS", "G_AS", "G_SA", "G_AA")
gate_values <- ablations |>
  filter(method_id %in% gate_ids) |>
  select(dataset_key, dataset, method_id, method, metric, value, direction, denominator, units)
write.csv(gate_values, file.path(table_dir, "selective_vs_forced_metrics_long.csv"), row.names = FALSE, na = "")

gate_wide <- gate_values |>
  select(dataset_key, dataset, metric, direction, method_id, value) |>
  pivot_wider(names_from = method_id, values_from = value)
gate_contrasts <- gate_wide |>
  mutate(
    direction_sign = ifelse(direction == "lower", -1, ifelse(direction == "higher", 1, NA_real_)),
    selective_drift_benefit = direction_sign * (G_SS - G_AS),
    selective_batch_benefit = direction_sign * (G_SS - G_SA),
    fully_selective_benefit = direction_sign * (G_SS - G_AA),
    contrast_definition = "Positive values favor selective correction over the corresponding forced-correction condition."
  ) |>
  select(-direction_sign)
write.csv(gate_contrasts, file.path(table_dir, "selective_vs_forced_contrasts.csv"), row.names = FALSE, na = "")

read_seed_metrics <- function(family) {
  bind_rows(lapply(sprintf("SIM%02d", 1:10), function(seed_id) {
    path <- file.path(
      release_root, "results", "evaluation", "simulation_seeds", seed_id,
      family, "method_metrics.csv"
    )
    value <- read_required(path)
    if (!"seed_id" %in% names(value)) value$seed_id <- seed_id
    value$method <- normalise_method(value$method)
    value
  }))
}
seed_methods <- read_seed_metrics("methods")
seed_ablations <- read_seed_metrics("ablations")

runtime_methods <- read_required(file.path(
  release_root, "analysis", "config", "method_manifest.csv"
))
simulation_runtime_tasks <- read_required(file.path(
  release_root, "analysis", "config", "simulation_task_manifest.csv"
))
seed_runtime_policy <- bind_rows(lapply(sprintf("SIM%02d", 1:10), function(seed_id) {
  bind_rows(lapply(seq_len(nrow(runtime_methods)), function(index) {
    method_id <- as.character(runtime_methods$method_id[index])
    manifest <- jsonlite::read_json(
      file.path(
        release_root, "results", "simulation_seeds", seed_id, method_id,
        "run_manifest.json"
      ),
      simplifyVector = TRUE
    )
    decision <- json_scalar(manifest, "decision")
    execution_runtime <- as.numeric(json_scalar(manifest, "runtime_sec"))
    reportable <- method_id == "raw" || decision != "reuse_exact"
    data.frame(
      analysis_family = "simulation_seed",
      task_id = simulation_runtime_tasks$task_id[
        match(
          paste(seed_id, method_id),
          paste(simulation_runtime_tasks$seed_id, simulation_runtime_tasks$method_id)
        )
      ],
      dataset_key = "simulation",
      seed_id = seed_id, method_id = method_id, method = as.character(runtime_methods$method[index]),
      reuse_decision = decision,
      execution_runtime_sec = execution_runtime,
      reported_method_runtime_sec = if (method_id == "raw") 0 else if (reportable) execution_runtime else NA_real_,
      runtime_provenance = if (method_id == "raw") {
        "no correction"
      } else if (reportable) {
        "measured method execution in frozen release"
      } else {
        "not remeasured; cache retrieval time excluded"
      },
      stringsAsFactors = FALSE
    )
  }))
}))
seed_runtime_rows <- seed_methods$metric == "runtime_sec"
seed_policy_index <- match(
  paste(seed_methods$seed_id, seed_methods$method_id),
  paste(seed_runtime_policy$seed_id, seed_runtime_policy$method_id)
)
seed_methods$value[seed_runtime_rows] <-
  seed_runtime_policy$reported_method_runtime_sec[seed_policy_index[seed_runtime_rows]]
seed_methods$notes[seed_runtime_rows] <-
  seed_runtime_policy$runtime_provenance[seed_policy_index[seed_runtime_rows]]
write.csv(
  bind_rows(primary_runtime_policy, seed_runtime_policy),
  file.path(table_dir, "runtime_measurement_policy.csv"),
  row.names = FALSE, na = ""
)
write.csv(seed_methods, file.path(table_dir, "simulation_seed_method_metrics.csv"), row.names = FALSE, na = "")
write.csv(seed_ablations, file.path(table_dir, "simulation_seed_ablation_metrics.csv"), row.names = FALSE, na = "")

summarise_values <- function(data, groups) {
  data |>
    group_by(across(all_of(groups))) |>
    summarise(
      n = sum(is.finite(value)), mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE), median = median(value, na.rm = TRUE),
      q1 = quantile(value, 0.25, na.rm = TRUE), q3 = quantile(value, 0.75, na.rm = TRUE),
      minimum = min(value, na.rm = TRUE), maximum = max(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(across(c(mean, sd, median, q1, q3, minimum, maximum), ~ifelse(is.nan(.x) | is.infinite(.x), NA_real_, .x)))
}
seed_summary <- summarise_values(
  seed_methods, c("method_id", "method", "metric", "direction", "units")
)
write.csv(seed_summary, file.path(table_dir, "simulation_seed_summary.csv"), row.names = FALSE, na = "")

reference_metrics <- bind_rows(lapply(dataset_order, function(dataset_key) {
  bind_rows(lapply(sprintf("SPLIT%02d", 1:10), function(split_id) {
    read_required(file.path(
      release_root, "results", "evaluation", "reference_splits",
      dataset_key, split_id, "method_metrics.csv"
    ))
  }))
}))
reference_metrics$method <- normalise_method(reference_metrics$method)
reference_metrics$dataset <- unname(dataset_labels[reference_metrics$dataset_key])
write.csv(reference_metrics, file.path(table_dir, "reference_split_method_metrics.csv"), row.names = FALSE, na = "")
reference_summary <- summarise_values(
  reference_metrics,
  c("dataset_key", "dataset", "method_id", "method", "metric", "direction", "units")
)
write.csv(reference_summary, file.path(table_dir, "reference_split_summary.csv"), row.names = FALSE, na = "")

reference_tasks <- read_required(file.path(
  release_root, "analysis", "config", "reference_split_task_manifest.csv"
))
selection_rows <- vector("list", nrow(reference_tasks))
candidate_rows <- list()
candidate_row_index <- 0L
for (index in seq_len(nrow(reference_tasks))) {
  task <- reference_tasks[index, , drop = FALSE]
  directory <- file.path(release_root, as.character(task$output_dir))
  manifest <- jsonlite::read_json(
    file.path(directory, "run_manifest.json"), simplifyVector = TRUE
  )
  details <- readRDS(file.path(directory, "selected_details.rds"))
  candidate_path <- file.path(directory, "candidate_scores.csv")
  selected_candidate_id <- NA_character_
  selection_signature <- as.character(task$selection_rule)
  selected_training_qc_score <- NA_real_
  candidate_count <- 1L

  if (file.exists(candidate_path)) {
    candidates <- read_required(candidate_path)
    selected_from_details <- if (!is.null(details$candidate_id)) {
      as.character(details$candidate_id)
    } else {
      character()
    }
    candidates$selected <- if ("selected" %in% names(candidates) &&
                               is.logical(candidates$selected)) {
      candidates$selected
    } else if ("selected" %in% names(candidates)) {
      tolower(as.character(candidates$selected)) == "true"
    } else if (length(selected_from_details) == 1L) {
      as.character(candidates$candidate_id) == selected_from_details
    } else {
      stop("Candidate table lacks a selection marker and selected_details lacks candidate_id: ", task$task_id)
    }
    selected <- candidates[candidates$selected %in% TRUE, , drop = FALSE]
    if (nrow(selected) != 1L) {
      stop("Reference task does not have exactly one selected candidate: ", task$task_id)
    }
    selected_candidate_id <- as.character(selected$candidate_id)
    selection_signature <- paste0(
      "candidate=", selected_candidate_id, "; ", as.character(selected$parameters)
    )
    selected_training_qc_score <- as.numeric(selected$qc_score)
    candidate_count <- nrow(candidates)
    candidates$dataset_key <- as.character(task$dataset_key)
    candidates$dataset <- unname(dataset_labels[as.character(task$dataset_key)])
    candidates$split_id <- as.character(task$split_id)
    candidates$task_id <- as.character(task$task_id)
    candidates$method_id <- as.character(task$method_id)
    candidates$method <- as.character(task$method)
    candidates$hidden_qc_values_available <- FALSE
    candidates$biological_labels_available <- FALSE
    candidates$replicate_identities_available <- FALSE
    candidate_row_index <- candidate_row_index + 1L
    candidate_rows[[candidate_row_index]] <- candidates
  } else if (!is.null(details$selected_parameters)) {
    parameters <- details$selected_parameters
    keep <- intersect(
      c(
        "parameters", "batch_option", "pelt_penalty", "test",
        "ljung_box_fitdf", "lag", "spline_method", "acorr_fdr",
        "anova_fdr", "normalization", "scale_by_batch"
      ),
      names(parameters)
    )
    selection_signature <- jsonlite::toJSON(
      parameters[keep], auto_unbox = TRUE, null = "null", digits = NA
    )
    selected_training_qc_score <- if (!is.null(parameters$quality_score)) {
      as.numeric(parameters$quality_score)
    } else {
      NA_real_
    }
  }

  selection_rows[[index]] <- data.frame(
    dataset_key = as.character(task$dataset_key),
    dataset = unname(dataset_labels[as.character(task$dataset_key)]),
    split_id = as.character(task$split_id),
    task_id = as.character(task$task_id),
    method_id = as.character(task$method_id),
    method = as.character(task$method),
    selection_rule = as.character(task$selection_rule),
    selected_candidate_id = selected_candidate_id,
    selection_signature = as.character(selection_signature),
    selected_training_qc_score = selected_training_qc_score,
    candidate_count = candidate_count,
    n_training_qc = as.integer(manifest$n_training_qc),
    n_hidden_qc = as.integer(manifest$n_hidden_qc),
    hidden_qc_values_available = FALSE,
    biological_labels_available = FALSE,
    replicate_identities_available = FALSE,
    stringsAsFactors = FALSE
  )
}
reference_selection <- bind_rows(selection_rows)
write.csv(
  reference_selection,
  file.path(table_dir, "reference_split_selection.csv"),
  row.names = FALSE, na = ""
)
reference_selection_stability <- reference_selection |>
  count(
    dataset_key, dataset, method_id, method, selection_signature,
    name = "splits_selected"
  ) |>
  group_by(dataset_key, method_id) |>
  mutate(
    total_splits = sum(splits_selected),
    selection_proportion = splits_selected / total_splits
  ) |>
  ungroup() |>
  arrange(factor(dataset_key, levels = dataset_order), factor(method, levels = method_order), desc(selection_proportion))
write.csv(
  reference_selection_stability,
  file.path(table_dir, "reference_split_selection_stability.csv"),
  row.names = FALSE, na = ""
)
reference_candidates <- bind_rows(candidate_rows)
write.csv(
  reference_candidates,
  file.path(table_dir, "reference_split_candidate_scores.csv"),
  row.names = FALSE, na = ""
)

floor_task_table <- bind_rows(
  primary_tasks |>
    transmute(
      analysis_family = "primary", task_id, dataset_key,
      seed_id = NA_character_, split_id = NA_character_, method_id, method,
      output_dir, details_file = "method_details.rds"
    ),
  simulation_runtime_tasks |>
    transmute(
      analysis_family = "simulation_seed", task_id,
      dataset_key = "simulation", seed_id, split_id = NA_character_,
      method_id, method, output_dir, details_file = "method_details.rds"
    ),
  reference_tasks |>
    transmute(
      analysis_family = "reference_split", task_id, dataset_key,
      seed_id = NA_character_, split_id, method_id, method, output_dir,
      details_file = "selected_details.rds"
    )
)
floor_rows <- lapply(seq_len(nrow(floor_task_table)), function(index) {
  task <- floor_task_table[index, , drop = FALSE]
  details <- readRDS(file.path(
    release_root, as.character(task$output_dir), as.character(task$details_file)
  ))
  record <- details$intensity_floor
  applied <- !is.null(record) && isTRUE(record$applied)
  data.frame(
    task[, setdiff(names(task), "details_file"), drop = FALSE],
    floor_applied = applied,
    negative_values_floored = if (applied) as.integer(record$n_values) else 0L,
    features_with_flooring = if (applied) as.integer(record$n_features) else 0L,
    samples_with_flooring = if (applied) as.integer(record$n_samples) else 0L,
    minimum_before_floor = if (applied) as.numeric(record$minimum_before_floor) else NA_real_,
    rule = if (applied) as.character(record$rule) else "No negative corrected intensities.",
    stringsAsFactors = FALSE
  )
})
intensity_flooring <- bind_rows(floor_rows)
write.csv(
  intensity_flooring,
  file.path(table_dir, "intensity_domain_safeguard.csv"),
  row.names = FALSE, na = ""
)

tiger_task_table <- floor_task_table |>
  filter(method_id == "tiger")
tiger_fallback_rows <- lapply(seq_len(nrow(tiger_task_table)), function(index) {
  task <- tiger_task_table[index, , drop = FALSE]
  details <- readRDS(file.path(
    release_root, as.character(task$output_dir), as.character(task$details_file)
  ))
  record <- details$tiger_nonfinite_fallback
  applied <- !is.null(record) && isTRUE(record$applied)
  data.frame(
    task[, setdiff(names(task), "details_file"), drop = FALSE],
    fallback_applied = applied,
    nonfinite_values_returned = if (applied) as.integer(record$n_nonfinite_values) else 0L,
    features_restored_to_raw = if (applied) as.integer(record$n_features) else 0L,
    samples_with_nonfinite_output = if (applied) as.integer(record$n_samples) else 0L,
    rule = if (applied) as.character(record$rule) else "TIGER output was finite.",
    stringsAsFactors = FALSE
  )
})
tiger_fallbacks <- bind_rows(tiger_fallback_rows)
write.csv(
  tiger_fallbacks,
  file.path(table_dir, "tiger_nonfinite_feature_fallback.csv"),
  row.names = FALSE, na = ""
)

partial_manifest <- read_required(file.path(
  release_root, "analysis", "config", "partial_confounding_task_manifest.csv"
))
partial_metrics <- bind_rows(lapply(partial_manifest$scenario_id, function(scenario_id) {
  read_required(file.path(
    release_root, "results", "partial_confounding", scenario_id, "method_metrics.csv"
  ))
}))
partial_metrics$method <- normalise_method(partial_metrics$method)
partial_metrics <- left_join(
  partial_metrics,
  partial_manifest[, c("scenario_id", "confounding_type", "level")],
  by = "scenario_id"
)
write.csv(partial_metrics, file.path(table_dir, "partial_confounding_method_metrics.csv"), row.names = FALSE, na = "")
partial_summary <- summarise_values(
  partial_metrics,
  c("confounding_type", "level", "method", "metric", "metric_direction", "units")
)
write.csv(partial_summary, file.path(table_dir, "partial_confounding_summary.csv"), row.names = FALSE, na = "")

zero_reference <- partial_metrics |>
  filter(confounding_type == "none", level == 0) |>
  select(seed_id, method, metric, zero_value = value)
partial_deltas <- partial_metrics |>
  left_join(zero_reference, by = c("seed_id", "method", "metric")) |>
  mutate(delta_vs_paired_zero = value - zero_value)
write.csv(partial_deltas, file.path(table_dir, "partial_confounding_paired_zero_deltas.csv"), row.names = FALSE, na = "")

theme_publication <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#E2E5E9", linewidth = 0.3, linetype = "dotted"),
      axis.text = element_text(colour = "#343A40"),
      axis.title = element_text(colour = "#202428", face = "bold"),
      strip.background = element_rect(fill = "#EEF1F4", colour = NA),
      strip.text = element_text(colour = "#202428", face = "bold"),
      plot.title = element_text(colour = "#202428", face = "bold"),
      legend.position = "top", legend.title = element_text(face = "bold")
    )
}

save_figure <- function(plot, stem, width, height, source) {
  write.csv(source, file.path(source_dir, paste0(stem, ".csv")), row.names = FALSE, na = "")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 320)
}

technical_metrics <- c(
  heldout_qc_cv_mean = "Held-out QC CV (%)",
  residual_run_order_gam_mean_deviance = "Residual GAM deviance",
  residual_ljung_box_proportion_significant = "Ljung-Box positive",
  batch_weighted_pc_r2_categorical = "Batch weighted-PC R²"
)
technical_source <- primary |>
  filter(metric %in% names(technical_metrics)) |>
  mutate(
    dataset = factor(as.character(dataset), levels = unname(dataset_labels)),
    method = factor(as.character(method), levels = method_order),
    metric_label = factor(unname(technical_metrics[metric]), levels = unname(technical_metrics))
  )
technical_plot <- ggplot(technical_source, aes(method, value, colour = method)) +
  geom_point(size = 1.6) +
  geom_line(aes(group = 1), colour = "#C8CDD2", linewidth = 0.35) +
  facet_grid(metric_label ~ dataset, scales = "free_y") +
  scale_colour_manual(values = method_palette, drop = FALSE) +
  labs(title = "Technical-effect metrics across the six benchmarks", x = NULL, y = NULL, colour = "Method") +
  theme_publication(8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_figure(technical_plot, "primary_technical_overview", 13, 8, technical_source)

biology_names <- c(
  truth_metabolite_profile_icc_mean = "Truth metabolite ICC",
  truth_sample_profile_pearson_mean = "Truth sample Pearson",
  metabolite_icc_a1_median = "Metabolite ICC(A,1)",
  feature_repeatability_ratio_median = "Feature repeatability ratio",
  sample_replicate_pearson_median = "Sample Pearson",
  genuine_replicate_icc_median = "Replicate-profile ICC(A,1)",
  biological_weighted_pc_r2 = "Biological weighted-PC R²",
  cross_batch_effect_pearson_median = "Cross-batch effect Pearson"
)
biology_source <- primary |>
  filter(metric %in% names(biology_names)) |>
  mutate(
    dataset = factor(as.character(dataset), levels = unname(dataset_labels)),
    method = factor(as.character(method), levels = method_order),
    metric_label = unname(biology_names[metric])
  )
biology_plot <- ggplot(biology_source, aes(method, value, colour = method)) +
  geom_point(size = 1.6) +
  facet_grid(metric_label ~ dataset, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = method_palette, drop = FALSE) +
  labs(title = "Truth recovery and biological preservation", x = NULL, y = NULL, colour = "Method") +
  theme_publication(8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_figure(biology_plot, "primary_preservation_overview", 13, 9, biology_source)

impact_plot_source <- stage_source |>
  filter(table_id %in% c("heldout_qc_cv", "gam_deviance", "batch_r2")) |>
  mutate(
    dataset = factor(dataset, levels = unname(dataset_labels)),
    stage = factor(stage, levels = c("Outlier_shrinkage", "Selective_drift", "Selective_batch", "PQN")),
    metric_label = factor(
      table_id,
      levels = c("heldout_qc_cv", "gam_deviance", "batch_r2"),
      labels = c("Held-out QC CV", "Residual GAM deviance", "Batch weighted-PC R²")
    )
  )
impact_limit <- max(abs(impact_plot_source$favorable_impact), na.rm = TRUE)
impact_plot <- ggplot(impact_plot_source, aes(stage, dataset, fill = favorable_impact)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  facet_wrap(~metric_label, scales = "free", ncol = 1) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0,
    limits = c(-impact_limit, impact_limit), oob = scales::squish
  ) +
  labs(
    title = "Incremental contribution of WiNN stages",
    subtitle = "Positive values favor adding the stage", x = NULL, y = NULL,
    fill = "Favorable\nchange"
  ) +
  theme_publication(9) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_figure(impact_plot, "ablation_stage_impacts", 8, 8, impact_plot_source)

gate_figure_metrics <- c(
  "heldout_qc_cv_mean", "residual_run_order_gam_mean_deviance",
  "batch_weighted_pc_r2_categorical", "metabolite_icc_a1_median",
  "feature_repeatability_ratio_median", "sample_replicate_pearson_median",
  "truth_metabolite_profile_icc_mean", "truth_sample_profile_pearson_mean"
)
gate_metric_labels <- c(
  heldout_qc_cv_mean = "Held-out QC CV (%)",
  residual_run_order_gam_mean_deviance = "Residual GAM deviance",
  batch_weighted_pc_r2_categorical = "Batch weighted-PC R²",
  metabolite_icc_a1_median = "Metabolite ICC(A,1)",
  feature_repeatability_ratio_median = "Feature repeatability ratio",
  sample_replicate_pearson_median = "Sample Pearson",
  truth_metabolite_profile_icc_mean = "Truth metabolite ICC",
  truth_sample_profile_pearson_mean = "Truth sample Pearson"
)
gate_source <- gate_values |>
  filter(metric %in% gate_figure_metrics) |>
  mutate(
    dataset = factor(dataset, levels = unname(dataset_labels)),
    metric_label = factor(
      unname(gate_metric_labels[metric]),
      levels = unname(gate_metric_labels[gate_figure_metrics])
    ),
    gate = factor(
      method_id, levels = gate_ids,
      labels = c("Selective/Selective", "All/Selective", "Selective/All", "All/All")
    )
  )
gate_palette <- c(
  "Selective/Selective" = "#4477AA", "All/Selective" = "#EE6677",
  "Selective/All" = "#228833", "All/All" = "#AA3377"
)
gate_plot <- ggplot(gate_source, aes(gate, value, colour = gate)) +
  geom_point(size = 1.6) +
  facet_grid(metric_label ~ dataset, scales = "free_y") +
  scale_colour_manual(values = gate_palette, drop = FALSE) +
  labs(title = "Selective and forced WiNN correction", x = NULL, y = NULL, colour = "Drift/Batch gate") +
  theme_publication(8) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_figure(gate_plot, "selective_vs_forced", 13, 11, gate_source)

seed_plot_metrics <- c(
  heldout_qc_cv_mean = "Held-out QC CV (%)",
  residual_run_order_gam_mean_deviance = "Residual GAM deviance",
  batch_weighted_pc_r2_categorical = "Batch weighted-PC R²",
  truth_metabolite_profile_icc_mean = "Truth metabolite ICC",
  truth_sample_profile_pearson_mean = "Truth sample Pearson"
)
seed_source <- seed_methods |>
  filter(metric %in% names(seed_plot_metrics)) |>
  mutate(
    method = factor(method, levels = method_order),
    metric_label = factor(unname(seed_plot_metrics[metric]), levels = unname(seed_plot_metrics))
  )
seed_plot <- ggplot(seed_source, aes(method, value, fill = method)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, linewidth = 0.35) +
  geom_point(aes(group = seed_id), position = position_jitter(width = 0.12), size = 0.55, alpha = 0.65) +
  facet_wrap(~metric_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = method_palette, drop = FALSE) +
  labs(title = "Ten-seed simulation stability", x = NULL, y = NULL, fill = "Method") +
  theme_publication(8) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_figure(seed_plot, "simulation_seed_stability", 11, 7, seed_source)

reference_source <- reference_metrics |>
  filter(metric %in% c("heldout_qc_cv_mean", "residual_run_order_gam_mean_deviance")) |>
  mutate(
    dataset = factor(dataset, levels = unname(dataset_labels)),
    method = factor(method, levels = method_order),
    metric_label = ifelse(metric == "heldout_qc_cv_mean", "Held-out QC CV (%)", "Residual GAM deviance")
  )
reference_plot <- ggplot(reference_source, aes(method, value, fill = method)) +
  geom_boxplot(width = 0.68, outlier.shape = NA, linewidth = 0.3) +
  facet_grid(metric_label ~ dataset, scales = "free_y") +
  scale_fill_manual(values = method_palette, drop = FALSE) +
  labs(title = "Ten-split reference stability", x = NULL, y = NULL, fill = "Method") +
  theme_publication(8) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_figure(reference_plot, "reference_split_stability", 13, 5.8, reference_source)

partial_plot_metrics <- c(
  clean_sample_profile_pearson_mean = "Clean sample-profile Pearson",
  median_attenuation_ratio_injected_responsive = "Responsive-effect ratio",
  batch_weighted_pc_r2_categorical = "Batch weighted-PC R²",
  residual_run_order_gam_mean_deviance = "Residual GAM deviance"
)
partial_source <- partial_summary |>
  filter(confounding_type != "none", metric %in% names(partial_plot_metrics)) |>
  mutate(
    method = factor(method, levels = c("Raw", "WiNN fixed default (no QC)")),
    metric_label = factor(unname(partial_plot_metrics[metric]), levels = unname(partial_plot_metrics)),
    confounding_type = factor(confounding_type, levels = c("batch", "run_order", "joint"))
  )
partial_plot <- ggplot(partial_source, aes(level, median, colour = method, group = method)) +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = method), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
  facet_grid(metric_label ~ confounding_type, scales = "free_y") +
  scale_colour_manual(values = method_palette, drop = FALSE) +
  scale_fill_manual(values = method_palette, drop = FALSE) +
  labs(
    title = "WiNN under partial technical-biological confounding",
    subtitle = "Lines show medians; ribbons show interquartile ranges across 10 seeds",
    x = "Confounding strength", y = NULL, colour = "Method", fill = "Method"
  ) + theme_publication(8)
save_figure(partial_plot, "partial_confounding", 10, 8, partial_source)

figure_manifest <- data.frame(
  figure = sub("[.]csv$", "", basename(list.files(source_dir, pattern = "[.]csv$", full.names = TRUE))),
  source_csv = file.path("figures", "source_data", basename(list.files(source_dir, pattern = "[.]csv$", full.names = TRUE))),
  stringsAsFactors = FALSE
)
figure_manifest$pdf <- ifelse(
  file.exists(file.path(figure_dir, paste0(figure_manifest$figure, ".pdf"))),
  file.path("figures", paste0(figure_manifest$figure, ".pdf")), ""
)
figure_manifest$png <- ifelse(
  file.exists(file.path(figure_dir, paste0(figure_manifest$figure, ".png"))),
  file.path("figures", paste0(figure_manifest$figure, ".png")), ""
)
write.csv(figure_manifest, file.path(figure_dir, "figure_manifest.csv"), row.names = FALSE)

utils::capture.output(sessionInfo(), file = file.path(final_dir, "sessionInfo.txt"))
jsonlite::write_json(
  list(
    schema = "endpoint_free_winn_0.1.4_final_aggregate_v1",
    package_version = "0.1.4",
    package_sha256 = "71a0964cee2778b2e5789d20621147e074c7945e813cf76af2ceeb696104aae1",
    datasets = dataset_order, methods = method_order,
    task_completion = completion,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ),
  file.path(final_dir, "aggregate_manifest.json"), auto_unbox = TRUE, pretty = TRUE
)

all_final_files <- list.files(final_dir, recursive = TRUE, full.names = TRUE)
all_final_files <- all_final_files[
  file.info(all_final_files)$isdir %in% FALSE &
    basename(all_final_files) != "file_manifest.csv"
]
checksum_manifest <- data.frame(
  path = substring(all_final_files, nchar(final_dir) + 2L),
  bytes = as.numeric(file.info(all_final_files)$size),
  sha256 = vapply(all_final_files, digest::digest, character(1), file = TRUE, algo = "sha256"),
  stringsAsFactors = FALSE
)
write.csv(checksum_manifest, file.path(final_dir, "file_manifest.csv"), row.names = FALSE)

message("Final endpoint-free release aggregation completed: ", final_dir)
