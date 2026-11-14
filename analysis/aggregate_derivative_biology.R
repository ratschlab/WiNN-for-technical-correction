#!/usr/bin/env Rscript

# Aggregate derivative-only evaluations. No correction or selection code is sourced.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
required_packages <- c("digest", "jsonlite", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
source(file.path(repo_root, "analysis", "derivative_biology_helpers.R"), local = TRUE)

root <- file.path(repo_root, "results", "final", "derivative")
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")
source_dir <- file.path(figure_dir, "source_data")
provenance_dir <- file.path(root, "provenance")
for (path in c(table_dir, figure_dir, source_dir, provenance_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

datasets <- derivative_dataset_registry()
methods <- derivative_method_registry()
config_path <- file.path(repo_root, "analysis", "config", "derivative_biology_config.json")
config <- jsonlite::read_json(config_path, simplifyVector = TRUE)
config_sha <- digest::digest(config_path, file = TRUE, algo = "sha256")

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing derivative artifact: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

bind_rows_base <- function(values) {
  values <- Filter(function(value) !is.null(value) && nrow(value), values)
  if (!length(values)) return(data.frame())
  columns <- unique(unlist(lapply(values, names), use.names = FALSE))
  values <- lapply(values, function(value) {
    for (column in setdiff(columns, names(value))) value[[column]] <- NA
    value[columns]
  })
  result <- do.call(rbind, values); rownames(result) <- NULL; result
}

read_dataset_file <- function(filename, required = TRUE) {
  values <- lapply(seq_len(nrow(datasets)), function(index) {
    path <- file.path(root, datasets$dataset_key[index], filename)
    if (!file.exists(path)) {
      if (required) stop("Missing derivative artifact: ", path, call. = FALSE)
      return(NULL)
    }
    value <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    value$dataset_key <- datasets$dataset_key[index]
    value$dataset <- datasets$dataset[index]
    value
  })
  bind_rows_base(values)
}

run_manifests <- read_dataset_file("derivative_run_manifest.csv")
if (any(run_manifests$test_mode) || any(run_manifests$bootstrap_replicates != 2000L) ||
    any(run_manifests$config_sha256 != config_sha) ||
    any(run_manifests$correction_methods_executed) || any(run_manifests$tuning_executed) ||
    any(run_manifests$selection_executed)) {
  stop("Derivative run manifests do not satisfy the frozen production design.", call. = FALSE)
}
expected_permutations <- run_manifests$dataset_key %in% c("waveica", "mtbls79", "sacurine")
if (any(run_manifests$permutation_replicates[expected_permutations] != 5000L) ||
    any(run_manifests$permutation_replicates[!expected_permutations] != 0L)) {
  stop("Derivative permutation counts do not match the prespecified design.", call. = FALSE)
}

bootstrap_estimates <- read_dataset_file("bootstrap_estimates.csv")
bootstrap_differences <- read_dataset_file("bootstrap_paired_differences.csv")
write.csv(bootstrap_estimates, file.path(table_dir, "bootstrap_estimates.csv"), row.names = FALSE, na = "")
write.csv(bootstrap_differences, file.path(table_dir, "bootstrap_paired_differences.csv"), row.names = FALSE, na = "")

point_metrics <- read_dataset_file("biological_point_metrics.csv", required = FALSE)
simulation_points <- bootstrap_estimates[bootstrap_estimates$dataset_key == "simulation", c(
  "dataset_key", "dataset", "method_id", "method", "metric", "point_estimate"
), drop = FALSE]
names(simulation_points)[names(simulation_points) == "point_estimate"] <- "value"
point_metrics <- bind_rows_base(list(point_metrics, simulation_points))
write.csv(point_metrics, file.path(table_dir, "biological_point_metrics.csv"), row.names = FALSE, na = "")

permutation_null <- bind_rows_base(lapply(c("mtbls79", "sacurine", "waveica"), function(dataset_key) {
  value <- read_required(file.path(root, dataset_key, "permutation_null.csv"))
  value$dataset_key <- dataset_key; value$dataset <- datasets$dataset[match(dataset_key, datasets$dataset_key)]; value
}))
permutation_summary <- bind_rows_base(lapply(c("mtbls79", "sacurine", "waveica"), function(dataset_key) {
  value <- read_required(file.path(root, dataset_key, "permutation_summary.csv"))
  value$dataset_key <- dataset_key; value$dataset <- datasets$dataset[match(dataset_key, datasets$dataset_key)]; value
}))
write.csv(permutation_null, file.path(table_dir, "permutation_null.csv"), row.names = FALSE, na = "")
write.csv(permutation_summary, file.path(table_dir, "permutation_summary.csv"), row.names = FALSE, na = "")

lobo_point <- bind_rows_base(lapply(c("mtbls79", "sacurine", "waveica"), function(dataset_key) {
  value <- read_required(file.path(root, dataset_key, "leave_one_batch_out_effects.csv"))
  value$dataset_key <- dataset_key; value$dataset <- datasets$dataset[match(dataset_key, datasets$dataset_key)]; value
}))
lobo_bootstrap <- bind_rows_base(lapply(c("mtbls79", "sacurine", "waveica"), function(dataset_key) {
  value <- read_required(file.path(root, dataset_key, "leave_one_batch_out_bootstrap_summary.csv"))
  value$dataset_key <- dataset_key; value$dataset <- datasets$dataset[match(dataset_key, datasets$dataset_key)]; value
}))
write.csv(lobo_point, file.path(table_dir, "leave_one_batch_out_effects.csv"), row.names = FALSE, na = "")
write.csv(lobo_bootstrap, file.path(table_dir, "leave_one_batch_out_bootstrap_summary.csv"), row.names = FALSE, na = "")

replicate_metrics <- bind_rows_base(lapply(c("batchcorr_set1", "clinical_fiams"), function(dataset_key) {
  value <- read_required(file.path(root, dataset_key, "replicate_metrics_by_identity.csv"))
  value$dataset_key <- dataset_key; value$dataset <- datasets$dataset[match(dataset_key, datasets$dataset_key)]; value
}))
matching_metrics <- bind_rows_base(lapply(c("batchcorr_set1", "clinical_fiams"), function(dataset_key) {
  value <- read_required(file.path(root, dataset_key, "matching_metrics_by_identity.csv"))
  value$dataset_key <- dataset_key; value$dataset <- datasets$dataset[match(dataset_key, datasets$dataset_key)]; value
}))
write.csv(replicate_metrics, file.path(table_dir, "replicate_metrics_by_identity.csv"), row.names = FALSE, na = "")
write.csv(matching_metrics, file.path(table_dir, "matching_metrics_by_identity.csv"), row.names = FALSE, na = "")

association_tables <- list()
for (dataset_key in c("mtbls79", "waveica")) {
  value <- read_required(file.path(root, dataset_key, "association_features.csv"))
  value$endpoint <- if (dataset_key == "mtbls79") "cow_vs_sheep" else "biological_group"
  value$variance_value <- value$within_group_variance
  value$dataset_key <- dataset_key; value$dataset <- datasets$dataset[match(dataset_key, datasets$dataset_key)]
  association_tables[[dataset_key]] <- value
}
sacurine_associations <- read_required(file.path(root, "sacurine", "sacurine_covariate_associations.csv"))
sacurine_associations$effect <- sacurine_associations$coefficient
sacurine_associations$variance_value <- sacurine_associations$residual_variance
sacurine_associations$dataset_key <- "sacurine"; sacurine_associations$dataset <- "Sacurine"
association_tables$sacurine <- sacurine_associations
associations <- bind_rows_base(association_tables)
write.csv(associations, file.path(table_dir, "biological_association_features.csv"), row.names = FALSE, na = "")

association_keys <- unique(associations[c("dataset_key", "dataset", "method_id", "method", "endpoint")])
association_summary <- do.call(rbind, lapply(seq_len(nrow(association_keys)), function(index) {
  key <- association_keys[index, , drop = FALSE]
  selected <- associations[
    associations$dataset_key == key$dataset_key & associations$method_id == key$method_id &
      associations$endpoint == key$endpoint, , drop = FALSE
  ]
  data.frame(
    key, tested_features = nrow(selected), significant_features = sum(selected$significant, na.rm = TRUE),
    median_absolute_effect = stats::median(abs(selected$effect), na.rm = TRUE),
    median_absolute_standardized_effect = stats::median(abs(selected$standardized_effect), na.rm = TRUE),
    median_within_group_or_residual_variance = stats::median(selected$variance_value, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(association_summary, file.path(table_dir, "biological_association_summary.csv"), row.names = FALSE, na = "")

metric_registry <- read.csv(
  file.path(repo_root, "analysis", "config", "derivative_metric_direction_registry.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
write.csv(metric_registry, file.path(table_dir, "metric_direction_registry.csv"), row.names = FALSE, na = "")

non_estimable <- data.frame(
  dataset_key = c(
    "simulation", "simulation", "mtbls79", "batchcorr_set1", "batchcorr_set1",
    "clinical_fiams", "clinical_fiams", "sacurine", "waveica"
  ),
  analysis = c(
    "phenotype-label permutation", "leave-one-batch-out biological effects",
    "continuous phenotype association", "phenotype association", "leave-one-batch-out phenotype effects",
    "phenotype association", "leave-one-batch-out phenotype effects",
    "replicate/sample-identity agreement", "replicate/sample-identity agreement"
  ),
  status = "not applicable by design",
  reason = c(
    "Clean simulation truth is observed directly; there is no empirical phenotype-label null to permute.",
    "Simulation uncertainty is paired by completed seed and evaluated against clean truth rather than empirical batch-held-out effects.",
    "The deposited biological endpoint used here is cow versus sheep; no additional continuous phenotype was prespecified.",
    "This dataset supports repeated accession identity, not a deposited phenotype-group endpoint.",
    "Cross-batch identity matching and replicate agreement are the design-supported validation endpoints.",
    "No clinical phenotype is included in the privacy-audited benchmark output.",
    "Cross-plate identity matching and replicate agreement are the permitted design-supported endpoints.",
    "There are no repeated participants in the deposited study design.",
    "There are no repeated biological identities in the deposited study design."
  ),
  stringsAsFactors = FALSE
)
non_estimable$dataset <- datasets$dataset[match(non_estimable$dataset_key, datasets$dataset_key)]
non_estimable <- non_estimable[c("dataset_key", "dataset", "analysis", "status", "reason")]
write.csv(non_estimable, file.path(table_dir, "non_estimable_analyses.csv"), row.names = FALSE, na = "")

palette <- c(
  "Raw" = "#5E6472", "ComBat" = "#287271", "QC-RLSC" = "#2A9D8F",
  "QC-RFSC" = "#577590", "TIGER" = "#E9C46A", "SERRF" = "#E76F51",
  "WiNN auto (QC)" = "#9B5DE5", "WiNN auto-batch (QC)" = "#6D597A",
  "WiNN fixed default (no QC)" = "#B08900"
)
theme_publication <- function(base_size = 9) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E5E7EB", linewidth = 0.3, linetype = "dotted"),
      axis.text = ggplot2::element_text(colour = "#374151"),
      axis.title = ggplot2::element_text(colour = "#111827", face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#F1F3F5", colour = NA),
      strip.text = ggplot2::element_text(face = "bold", colour = "#111827"),
      plot.title = ggplot2::element_text(face = "bold", colour = "#111827"),
      legend.position = "top", legend.title = ggplot2::element_text(face = "bold")
    )
}

save_figure <- function(plot, stem, width, height, source) {
  write.csv(source, file.path(source_dir, paste0(stem, ".csv")), row.names = FALSE, na = "")
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in")
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 320)
}

bootstrap_key <- paste(bootstrap_estimates$dataset_key, bootstrap_estimates$method_id, bootstrap_estimates$metric)
get_boot <- function(dataset_key, metric) {
  selected <- bootstrap_estimates[bootstrap_estimates$dataset_key == dataset_key & bootstrap_estimates$metric == metric, ]
  selected
}
tradeoff_specs <- data.frame(
  dataset_key = c(
    "simulation", "simulation", "mtbls79", "waveica", "sacurine", "sacurine", "sacurine",
    "batchcorr_set1", "batchcorr_set1", "clinical_fiams", "clinical_fiams"
  ),
  x_metric = c(
    "truth_log1p_rmse", "truth_log1p_rmse", rep("batch_weighted_pc_r2_categorical", 9)
  ),
  y_metric = c(
    "truth_metabolite_profile_icc_mean", "truth_sample_profile_pearson_mean",
    "biological_weighted_pc_r2", "biological_weighted_pc_r2",
    "age_weighted_pc_r2", "bmi_weighted_pc_r2", "gender_weighted_pc_r2",
    "sample_replicate_pearson_median", "cross_batch_top1_accuracy",
    "sample_replicate_pearson_median", "cross_batch_top1_accuracy"
  ),
  panel = c(
    "Simulation: truth ICC", "Simulation: truth profile correlation", "MTBLS79: sample type",
    "WaveICA: deposited group", "Sacurine: age", "Sacurine: BMI", "Sacurine: gender",
    "BatchCorr: replicate correlation", "BatchCorr: cross-batch matching",
    "Clinical: replicate correlation", "Clinical: cross-plate matching"
  ), stringsAsFactors = FALSE
)
tradeoff_source <- bind_rows_base(lapply(seq_len(nrow(tradeoff_specs)), function(index) {
  specification <- tradeoff_specs[index, ]
  x <- get_boot(specification$dataset_key, specification$x_metric)
  y <- get_boot(specification$dataset_key, specification$y_metric)
  y <- y[match(x$method_id, y$method_id), ]
  data.frame(
    dataset_key = specification$dataset_key, panel = specification$panel,
    method_id = x$method_id, method = x$method,
    x_metric = specification$x_metric, x = x$point_estimate, x_ci_lower = x$ci_lower, x_ci_upper = x$ci_upper,
    y_metric = specification$y_metric, y = y$point_estimate, y_ci_lower = y$ci_lower, y_ci_upper = y$ci_upper,
    stringsAsFactors = FALSE
  )
}))
tradeoff_source$method <- factor(tradeoff_source$method, levels = methods$method)
tradeoff_plot <- ggplot2::ggplot(tradeoff_source, ggplot2::aes(x, y, colour = method)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = y_ci_lower, ymax = y_ci_upper), width = 0, linewidth = 0.35, alpha = 0.75) +
  ggplot2::geom_segment(ggplot2::aes(x = x_ci_lower, xend = x_ci_upper, y = y, yend = y), linewidth = 0.35, alpha = 0.75) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_wrap(~panel, scales = "free", ncol = 3) +
  ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
  ggplot2::labs(
    title = "Technical removal and biological preservation from frozen outputs",
    subtitle = "Points are prespecified estimates; bars are paired 95% bootstrap intervals (2,000 biological-unit resamples)",
    x = "Residual technical structure or truth error (lower is better)",
    y = "Biological preservation or truth recovery (higher is better)", colour = "Method"
  ) + theme_publication(8)
save_figure(tradeoff_plot, "technical_biological_tradeoffs", 12, 10, tradeoff_source)

wave_source <- tradeoff_source[tradeoff_source$dataset_key == "waveica", , drop = FALSE]
wave_plot <- ggplot2::ggplot(wave_source, ggplot2::aes(x, y, colour = method)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = y_ci_lower, ymax = y_ci_upper), width = 0, linewidth = 0.5) +
  ggplot2::geom_segment(ggplot2::aes(x = x_ci_lower, xend = x_ci_upper, y = y, yend = y), linewidth = 0.5) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
  ggplot2::labs(
    title = "WaveICA technical-biological trade-off",
    subtitle = "Neutral deposited labels group_0/group_1; 2,000 batch-by-group-stratified bootstrap resamples",
    x = "Categorical batch: variance-weighted mean R² across principal components 1-5 (lower)",
    y = "Deposited group: variance-weighted mean R² across principal components 1-5 (higher)", colour = "Method"
  ) + theme_publication(9)
save_figure(wave_plot, "waveica_technical_biological_tradeoff", 8.5, 6, wave_source)

# Updated MTBLS79 overview in the structure of the submitted public-data
# figure. PCA scores and every bar are derived from the same frozen common
# matrix panel used above; no correction code is called.
pre_for_figures <- read_required(file.path(
  repo_root, "results", "frozen_endpoint_free_inputs", "provenance", "matrix_checksums_pre.csv"
))
mtbls_frozen <- derivative_load_matrices(repo_root, "mtbls79", pre_for_figures)
mtbls_meta <- derivative_load_metadata(repo_root, "mtbls79", mtbls_frozen$sample_ids)
mtbls_study <- which(mtbls_meta$is_study)
make_pca_source <- function(method_id, panel) {
  matrix_value <- log1p(mtbls_frozen$matrices[[method_id]][, mtbls_study, drop = FALSE])
  fit <- stats::prcomp(t(matrix_value), center = TRUE, scale. = FALSE, rank. = 2)
  explained <- 100 * fit$sdev[1:2]^2 / sum(fit$sdev^2)
  data.frame(
    row_type = "pca", panel = sprintf("%s: PC1 %.1f%%, PC2 %.1f%%", panel, explained[1], explained[2]), method_id = method_id,
    method = methods$method[match(method_id, methods$method_id)],
    sample_id = mtbls_meta$sample_id[mtbls_study], class = mtbls_meta$class[mtbls_study],
    batch = as.character(mtbls_meta$batch[mtbls_study]),
    x = fit$x[, 1], y = fit$x[, 2],
    x_label = sprintf("PC1 (%.1f%%)", explained[1]), y_label = sprintf("PC2 (%.1f%%)", explained[2]),
    metric = NA_character_, value = NA_real_, stringsAsFactors = FALSE
  )
}
mtbls_bars <- point_metrics[
  point_metrics$dataset_key == "mtbls79" &
    point_metrics$metric %in% c(
      "batch_weighted_pc_r2_categorical", "biological_weighted_pc_r2", "biological_associated_features"
    ),
]
mtbls_bars$panel <- c(
  batch_weighted_pc_r2_categorical = "C  Categorical batch R²",
  biological_weighted_pc_r2 = "D  Cow-sheep R²",
  biological_associated_features = "E  Cow-sheep associations"
)[mtbls_bars$metric]
raw_associations <- mtbls_bars$value[
  mtbls_bars$method_id == "raw" & mtbls_bars$metric == "biological_associated_features"
]
mtbls_bars$value[mtbls_bars$metric == "biological_associated_features"] <-
  mtbls_bars$value[mtbls_bars$metric == "biological_associated_features"] - raw_associations
mtbls_bars$row_type <- "bar"
mtbls_bars$sample_id <- NA_character_
mtbls_bars$class <- NA_character_
mtbls_bars$batch <- NA_character_
mtbls_bars$x <- match(mtbls_bars$method_id, methods$method_id)
mtbls_bars$y <- mtbls_bars$value
mtbls_bars$x_label <- "Method"
mtbls_bars$y_label <- "Saved metric value"
mtbls_figure_source <- bind_rows_base(list(
  make_pca_source("raw", "A  Raw profiles"),
  make_pca_source("winn_fixed_no_qc", "B  Fixed WiNN profiles"),
  mtbls_bars
))
mtbls_figure_source$colour_group <- ifelse(
  mtbls_figure_source$row_type == "pca", mtbls_figure_source$class, mtbls_figure_source$method
)
mtbls_pca_plot <- ggplot2::ggplot(
  mtbls_figure_source[mtbls_figure_source$row_type == "pca", ],
  ggplot2::aes(x, y, colour = class, shape = batch)
) +
  ggplot2::geom_point(size = 1.25, alpha = 0.8) +
  ggplot2::facet_wrap(~panel, scales = "free", ncol = 2) +
  ggplot2::scale_colour_manual(values = c(C = "#C05A6E", S = "#4F6D8A"), drop = TRUE) +
  ggplot2::scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 15, "4" = 18, "5" = 3, "6" = 4, "7" = 8, "8" = 1)) +
  ggplot2::labs(x = "PC1 score", y = "PC2 score", colour = "Class", shape = "Batch") +
  theme_publication(7)
mtbls_bar_source <- mtbls_figure_source[mtbls_figure_source$row_type == "bar", ]
mtbls_bar_source$method <- factor(mtbls_bar_source$method, levels = methods$method)
mtbls_bar_plot <- ggplot2::ggplot(mtbls_bar_source, ggplot2::aes(method, y, fill = method)) +
  ggplot2::geom_col(
    width = 0.72
  ) +
  ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 3) +
  ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
  ggplot2::labs(x = NULL, y = "Saved metric value", fill = "Method") +
  theme_publication(7) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1))
mtbls_plot <- (mtbls_pca_plot / mtbls_bar_plot) +
  patchwork::plot_layout(heights = c(1.05, 1)) +
  patchwork::plot_annotation(
    title = "MTBLS79: technical removal with preserved sample-type structure",
    subtitle = "PCA and method metrics use the same 2,488-feature frozen panel"
  )
save_figure(mtbls_plot, "mtbls79_fixed_benchmark", 10, 8, mtbls_figure_source)

clinical_source <- bootstrap_estimates[
  bootstrap_estimates$dataset_key == "clinical_fiams" &
    bootstrap_estimates$metric %in% c(
      "batch_weighted_pc_r2_categorical", "genuine_replicate_icc_median", "cross_batch_top1_accuracy"
    ),
]
clinical_source$method <- factor(clinical_source$method, levels = methods$method)
clinical_source$metric_label <- factor(
  clinical_source$metric,
  levels = c("batch_weighted_pc_r2_categorical", "genuine_replicate_icc_median", "cross_batch_top1_accuracy"),
  labels = c("Categorical plate R² (lower)", "Genuine replicate ICC(A,1) (higher)", "Cross-plate top-1 matching (higher)")
)
clinical_plot <- ggplot2::ggplot(clinical_source, ggplot2::aes(method, point_estimate, colour = method)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0.18, linewidth = 0.4) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~metric_label, scales = "free_y", ncol = 1) +
  ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
  ggplot2::labs(
    title = "Pair-aware clinical FIA-MS benchmark",
    subtitle = "Point estimates and 95% biological-sample cluster-bootstrap intervals",
    x = NULL, y = NULL, colour = "Method"
  ) + theme_publication(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
save_figure(clinical_plot, "clinical_pair_aware_derivative", 7.5, 8, clinical_source)

ablation_impacts_for_plot <- read_required(file.path(
  repo_root, "results", "final", "tables", "ablation_stage_impacts.csv"
))
batch_cumulative <- ablation_impacts_for_plot[
  ablation_impacts_for_plot$table_id == "batch_r2" & ablation_impacts_for_plot$stage_index == 1L,
  c("dataset_key", "dataset", "cum_raw", "cum_outlier", "cum_drift", "cum_batch", "cum_pqn")
]
batch_trajectory_source <- reshape(
  batch_cumulative,
  varying = c("cum_raw", "cum_outlier", "cum_drift", "cum_batch", "cum_pqn"),
  v.names = "value", timevar = "stage",
  times = c("Raw", "Outlier shrinkage", "Selective drift", "Selective batch", "PQN"),
  direction = "long"
)
rownames(batch_trajectory_source) <- NULL
batch_trajectory_source$stage <- factor(
  batch_trajectory_source$stage,
  levels = c("Raw", "Outlier shrinkage", "Selective drift", "Selective batch", "PQN")
)
dataset_palette <- c(
  "Simulation" = "#495057", "MTBLS79" = "#287271", "Clinical FIA-MS" = "#577590",
  "BatchCorr Set 1" = "#B08900", "Sacurine" = "#9B5DE5", "WaveICA" = "#E76F51"
)
batch_trajectory_plot <- ggplot2::ggplot(
  batch_trajectory_source,
  ggplot2::aes(stage, value, group = dataset, colour = dataset)
) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::scale_colour_manual(values = dataset_palette, drop = FALSE) +
  ggplot2::labs(
    title = "Categorical batch structure through the fixed workflow",
    subtitle = "Every point uses the same frozen matrix panel and factor-coded batch target",
    x = NULL, y = "Categorical batch: variance-weighted mean R² across principal components 1-5",
    colour = "Dataset"
  ) + theme_publication(9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
save_figure(batch_trajectory_plot, "categorical_batch_r2_trajectory", 8.5, 5.5, batch_trajectory_source)

observed_plot <- permutation_summary[c("dataset_key", "dataset", "method_id", "method", "endpoint", "observed_count")]
observed_plot$row_type <- "observed"
observed_plot$permutation <- NA_integer_
observed_plot$significant_features <- observed_plot$observed_count
null_plot <- permutation_null
null_plot$row_type <- "null"
null_plot$observed_count <- NA_real_
permutation_plot_source <- bind_rows_base(list(null_plot, observed_plot))
permutation_plot_source$method <- factor(permutation_plot_source$method, levels = methods$method)
null_plot$method <- factor(null_plot$method, levels = methods$method)
observed_plot$method <- factor(observed_plot$method, levels = methods$method)
endpoint_labels <- c(
  age = "Age", bmi = "BMI", gender = "Gender",
  biological_group = "Deposited group", cow_vs_sheep = "Cow vs sheep"
)
permutation_plot_source$panel <- paste(
  permutation_plot_source$dataset,
  unname(endpoint_labels[permutation_plot_source$endpoint]), sep = ": "
)
null_plot$panel <- paste(null_plot$dataset, unname(endpoint_labels[null_plot$endpoint]), sep = ": ")
observed_plot$panel <- paste(observed_plot$dataset, unname(endpoint_labels[observed_plot$endpoint]), sep = ": ")
permutation_plot <- ggplot2::ggplot(null_plot, ggplot2::aes(method, significant_features, fill = method)) +
  ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA, linewidth = 0.3) +
  ggplot2::geom_point(data = observed_plot, ggplot2::aes(method, observed_count), inherit.aes = FALSE, shape = 21, fill = "white", colour = "#111827", size = 2) +
  ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 3) +
  ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
  ggplot2::labs(
    title = "Permutation calibration of biological association counts",
    subtitle = "Null distributions use 5,000 prespecified constrained permutations; white points are observed counts",
    x = NULL, y = "BH-FDR-significant features", fill = "Method"
  ) + theme_publication(7) +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
save_figure(permutation_plot, "association_permutation_calibration", 12, 8, permutation_plot_source)

lobo_source <- lobo_bootstrap[lobo_bootstrap$metric %in% c("pearson", "spearman", "sign_concordance"), ]
lobo_source$method <- factor(lobo_source$method, levels = methods$method)
lobo_source$metric_label <- factor(
  lobo_source$metric,
  levels = c("pearson", "spearman", "sign_concordance"),
  labels = c("Pearson effect correlation", "Spearman effect correlation", "Effect sign concordance")
)
lobo_source$endpoint_label <- unname(endpoint_labels[lobo_source$endpoint])
lobo_source$panel <- paste(lobo_source$dataset, lobo_source$endpoint_label, sep = ": ")
lobo_plot <- ggplot2::ggplot(lobo_source, ggplot2::aes(method, bootstrap_mean, colour = method)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, linewidth = 0.35) +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::facet_grid(metric_label ~ panel, scales = "free_y") +
  ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
  ggplot2::labs(
    title = "Held-out-batch biological-effect replication",
    subtitle = "Fold-specific estimates with 95% biological-unit bootstrap intervals",
    x = NULL, y = NULL, colour = "Method"
  ) + theme_publication(7) +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
save_figure(lobo_plot, "leave_one_batch_out_effect_replication", 13, 8, lobo_source)

replicate_source <- bootstrap_estimates[
  bootstrap_estimates$dataset_key %in% c("batchcorr_set1", "clinical_fiams") &
    bootstrap_estimates$metric %in% c("sample_replicate_pearson_median", "genuine_replicate_icc_median", "cross_batch_top1_accuracy"),
  , drop = FALSE
]
replicate_source$method <- factor(replicate_source$method, levels = methods$method)
replicate_source$metric_label <- factor(
  replicate_source$metric,
  levels = c("cross_batch_top1_accuracy", "genuine_replicate_icc_median", "sample_replicate_pearson_median"),
  labels = c("Cross-batch/plate top-1 accuracy", "Genuine replicate ICC(A,1)", "Replicate-profile Pearson correlation")
)
replicate_plot <- ggplot2::ggplot(replicate_source, ggplot2::aes(method, point_estimate, colour = method)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, linewidth = 0.4) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_grid(metric_label ~ dataset, scales = "free_y") +
  ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
  ggplot2::labs(
    title = "Repeated-identity preservation",
    subtitle = "Point estimates and 95% accession/patient-cluster bootstrap intervals",
    x = NULL, y = NULL, colour = "Method"
  ) + theme_publication(8) +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
save_figure(replicate_plot, "replicate_preservation", 10, 8, replicate_source)

pre_path <- file.path(repo_root, "results", "frozen_endpoint_free_inputs", "provenance", "matrix_checksums_pre.csv")
post_path <- file.path(repo_root, "results", "frozen_endpoint_free_inputs", "provenance", "matrix_checksums_post.csv")
pre <- read_required(pre_path); post <- read_required(post_path)
comparison <- merge(pre, post, by = c("category", "path"), suffixes = c("_pre", "_post"), all = TRUE)
comparison$unchanged <- with(comparison,
  !is.na(sha256_pre) & !is.na(sha256_post) & bytes_pre == bytes_post & sha256_pre == sha256_post
)
if (nrow(comparison) != nrow(pre) || nrow(comparison) != nrow(post) || !all(comparison$unchanged)) {
  stop("Frozen matrix pre/post checksum comparison failed.", call. = FALSE)
}
pre$derivative_use <- ifelse(
  pre$category == "primary", "direct_matrix_evaluation",
  ifelse(pre$category == "simulation_seed", "saved_seed_metric_evaluation",
    ifelse(pre$category == "ablation", "cumulative_metric_provenance", "checksum_audit_only")
  )
)
write.csv(pre, file.path(provenance_dir, "frozen_matrix_manifest.csv"), row.names = FALSE, na = "")
write.csv(comparison, file.path(provenance_dir, "matrix_checksum_comparison.csv"), row.names = FALSE, na = "")
release_provenance <- data.frame(
  artifact = c(
    "WiNN package version", "WiNN package commit", "WiNN source tarball SHA-256",
    "Frozen analysis commit", "Frozen pre-derivative result file manifest SHA-256",
    "Frozen aggregate manifest SHA-256", "Frozen dataset manifest SHA-256",
    "Derivative configuration SHA-256", "Weighted-PC computation",
    "PCA solver equivalence-test SHA-256"
  ),
  value = c(
    config$provenance$package_version, config$provenance$package_commit,
    config$provenance$package_tarball_sha256,
    config$provenance$analysis_commit,
    digest::digest(file.path(repo_root, "results", "frozen_endpoint_free_inputs", "provenance",
                             "pre_derivative_file_manifest.csv"), file = TRUE, algo = "sha256"),
    digest::digest(file.path(repo_root, "results", "final", "aggregate_manifest.json"),
                   file = TRUE, algo = "sha256"),
    digest::digest(file.path(repo_root, "analysis", "config", "dataset_manifest.csv"),
                   file = TRUE, algo = "sha256"),
    config_sha,
    "Exact top-five eigendecomposition when the smaller matrix dimension is <=300; otherwise deterministic irlba truncated SVD with an exact fallback.",
    digest::digest(file.path(repo_root, "analysis", "test_derivative_pca_solver.R"), file = TRUE, algo = "sha256")
  ), stringsAsFactors = FALSE
)
write.csv(release_provenance, file.path(provenance_dir, "release_provenance.csv"), row.names = FALSE, na = "")

figure_manifest <- data.frame(
  figure = sub("[.]csv$", "", basename(list.files(source_dir, pattern = "[.]csv$", full.names = TRUE))),
  source_csv = file.path("figures", "source_data", basename(list.files(source_dir, pattern = "[.]csv$", full.names = TRUE))),
  stringsAsFactors = FALSE
)
figure_manifest$pdf <- file.path("figures", paste0(figure_manifest$figure, ".pdf"))
figure_manifest$png <- file.path("figures", paste0(figure_manifest$figure, ".png"))
write.csv(figure_manifest, file.path(figure_dir, "figure_manifest.csv"), row.names = FALSE, na = "")

utils::capture.output(sessionInfo(), file = file.path(root, "sessionInfo.txt"))
jsonlite::write_json(list(
  schema = paste0("endpoint_free_winn_", config$provenance$package_version,
                  "_derivative_biology_manifest_v1"),
  package_version = config$provenance$package_version,
  package_commit = config$provenance$package_commit,
  package_tarball_sha256 = config$provenance$package_tarball_sha256,
  analysis_commit = config$provenance$analysis_commit,
  configuration_sha256 = config_sha,
  bootstrap_replicates = 2000,
  permutation_replicates = 5000,
  weighted_pc_solver = "Exact top-five eigendecomposition for smaller matrices; deterministic irlba truncated SVD above dimension 300, with an exact fallback.",
  weighted_pc_solver_equivalence_test = "analysis/test_derivative_pca_solver.R",
  frozen_matrix_count = nrow(pre),
  frozen_matrices_unchanged = all(comparison$unchanged),
  correction_methods_executed = FALSE,
  tuning_executed = FALSE,
  selection_executed = FALSE,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
), file.path(root, "derivative_analysis_manifest.json"), auto_unbox = TRUE, pretty = TRUE)

all_files <- list.files(root, recursive = TRUE, full.names = TRUE)
all_files <- all_files[file.info(all_files)$isdir %in% FALSE & basename(all_files) != "derivative_file_manifest.csv"]
file_manifest <- data.frame(
  path = substring(all_files, nchar(root) + 2L), bytes = as.numeric(file.info(all_files)$size),
  sha256 = vapply(all_files, digest::digest, character(1), file = TRUE, algo = "sha256"),
  stringsAsFactors = FALSE
)
write.csv(file_manifest, file.path(root, "derivative_file_manifest.csv"), row.names = FALSE)

# Refresh the release-level checksum inventory after all derivative artifacts
# have been written. The pre-derivative manifest hash remains recorded above as
# frozen provenance; this operation only inventories the current final folder.
final_root <- file.path(repo_root, "results", "final")
all_final_files <- list.files(final_root, recursive = TRUE, full.names = TRUE)
all_final_files <- all_final_files[
  file.info(all_final_files)$isdir %in% FALSE & basename(all_final_files) != "file_manifest.csv"
]
final_file_manifest <- data.frame(
  path = substring(all_final_files, nchar(final_root) + 2L),
  bytes = as.numeric(file.info(all_final_files)$size),
  sha256 = vapply(all_final_files, digest::digest, character(1), file = TRUE, algo = "sha256"),
  stringsAsFactors = FALSE
)
write.csv(final_file_manifest, file.path(final_root, "file_manifest.csv"), row.names = FALSE)
cat("Aggregated derivative-only biological-preservation analyses.\n")
