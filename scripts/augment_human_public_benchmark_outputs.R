#!/usr/bin/env Rscript

augment_human_public_benchmark_outputs <- function(repo_root, dataset) {
  stopifnot(dataset %in% c("sacurine", "waveica_adenocarcinoma"))
  required <- c("dplyr", "ggplot2", "limma", "patchwork", "tidyr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing output-augmentation package(s): ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(tidyr)
  })
  source(file.path(repo_root, "scripts", "plot_theme.R"), local = TRUE)

  processed_dir <- file.path(repo_root, "data", "public", dataset, "processed")
  result_dir <- file.path(repo_root, "results", dataset)
  figure_dir <- file.path(result_dir, "figures")
  source_dir <- file.path(figure_dir, "source_data")
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  matrix_name <- if (dataset == "sacurine") "sacurine_intensity_matrix.rds" else "waveica_intensity_matrix.rds"
  metadata_name <- if (dataset == "sacurine") "sacurine_metadata.csv" else "waveica_metadata.csv"
  seed <- if (dataset == "sacurine") 20260810L else 20260811L
  x <- readRDS(file.path(processed_dir, matrix_name))
  meta <- read.csv(file.path(processed_dir, metadata_name), check.names = FALSE, stringsAsFactors = FALSE)
  meta$is_qc <- as.logical(meta$is_qc)
  meta <- meta[match(colnames(x), meta$sample_id), , drop = FALSE]
  holdout <- read.csv(file.path(result_dir, "qc_holdout_ids.csv"), stringsAsFactors = FALSE)$sample_id
  training <- read.csv(file.path(result_dir, "qc_training_ids.csv"), stringsAsFactors = FALSE)$sample_id
  meta$qc_role <- ifelse(meta$sample_id %in% holdout, "Hidden QC",
                         ifelse(meta$sample_id %in% training, "Training QC", "Study"))
  meta$display_group <- ifelse(
    meta$is_qc, meta$qc_role,
    if (dataset == "sacurine") paste0("Study: ", meta$gender) else paste0("Study: ", meta$biological_group)
  )
  write.csv(meta, file.path(source_dir, "acquisition_design.csv"), row.names = FALSE, quote = TRUE, na = "")
  theme_benchmark <- function() {
    theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"),
            legend.position = "bottom", plot.title = element_text(face = "bold"))
  }
  acquisition_palette <- c(
    `Hidden QC` = "#D55E00", `Training QC` = "#0072B2",
    `Study: Female` = "#E78AC3", `Study: Male` = "#1B9E77",
    `Study: group_0` = "#E69F00", `Study: group_1` = "#009E73"
  )
  acquisition_plot <- ggplot(meta, aes(global_run_order, as.numeric(factor(batch)), color = display_group)) +
    geom_vline(xintercept = head(tapply(meta$global_run_order, meta$batch, max), -1L) + 0.5, color = "grey55", linetype = 2) +
    geom_point(size = 1.35, alpha = 0.8) +
    scale_color_manual(values = acquisition_palette[unique(meta$display_group)]) +
    scale_y_continuous(breaks = seq_along(unique(meta$batch)), labels = unique(meta$batch)) +
    labs(title = paste(if (dataset == "sacurine") "Sacurine negative-ion" else "WaveICA adenocarcinoma", "acquisition design"),
         x = "Global run order", y = "Supplied batch", color = NULL) + theme_benchmark()
  ggsave(file.path(figure_dir, "acquisition_design.pdf"), acquisition_plot, width = 11, height = 4.8)

  # Randomization and QC-spacing source data and publication figure.
  study <- meta[!meta$is_qc, , drop = FALSE]
  qc <- meta[meta$is_qc, , drop = FALSE] |>
    arrange(batch, within_batch_order) |>
    group_by(batch) |>
    mutate(previous_qc_gap = within_batch_order - lag(within_batch_order)) |>
    ungroup()
  write.csv(qc[, c("sample_id", "batch", "within_batch_order", "global_run_order", "qc_role", "previous_qc_gap")],
            file.path(source_dir, "qc_spacing.csv"), row.names = FALSE, quote = TRUE, na = "")

  if (dataset == "sacurine") {
    composition <- study |> count(batch, gender, name = "n")
    sequence <- bind_rows(
      data.frame(sample_id = study$sample_id, batch = study$batch, within_batch_order = study$within_batch_order, variable = "age", value = study$age),
      data.frame(sample_id = study$sample_id, batch = study$batch, within_batch_order = study$within_batch_order, variable = "bmi", value = study$bmi)
    )
    write.csv(composition, file.path(source_dir, "randomization_composition.csv"), row.names = FALSE, quote = TRUE)
    write.csv(sequence, file.path(source_dir, "randomization_sequence.csv"), row.names = FALSE, quote = TRUE)
    p_comp <- ggplot(composition, aes(batch, n, fill = gender)) + geom_col(position = "fill") +
      scale_y_continuous(labels = scales::label_percent()) + labs(title = "Gender composition", x = NULL, y = "Proportion", fill = "Gender") + theme_benchmark()
    p_seq <- ggplot(sequence, aes(within_batch_order, value, color = batch)) + geom_point(alpha = 0.7, size = 1) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) + facet_grid(variable ~ batch, scales = "free_y") +
      labs(title = "Age and BMI along within-batch order", x = "Within-batch order", y = NULL, color = "Batch") + theme_benchmark()
  } else {
    composition <- study |> count(batch, biological_group, name = "n")
    sequence <- study[, c("sample_id", "batch", "within_batch_order", "biological_group")]
    sequence$group_numeric <- as.integer(factor(sequence$biological_group)) - 1L
    write.csv(composition, file.path(source_dir, "randomization_composition.csv"), row.names = FALSE, quote = TRUE)
    write.csv(sequence, file.path(source_dir, "randomization_sequence.csv"), row.names = FALSE, quote = TRUE)
    p_comp <- ggplot(composition, aes(batch, n, fill = biological_group)) + geom_col(position = "fill") +
      scale_y_continuous(labels = scales::label_percent()) + labs(title = "Biological-group composition", x = NULL, y = "Proportion", fill = "Group") + theme_benchmark()
    p_seq <- ggplot(sequence, aes(within_batch_order, group_numeric, color = biological_group)) + geom_point(alpha = 0.65, size = 0.9) +
      facet_wrap(~batch, ncol = 1) + scale_y_continuous(breaks = c(0, 1), labels = c("group_0", "group_1")) +
      labs(title = "Biological group along within-batch order", x = "Within-batch order", y = NULL, color = "Group") + theme_benchmark()
  }
  p_qc <- ggplot(qc, aes(within_batch_order, as.numeric(factor(batch)), color = qc_role)) +
    geom_point(size = 1.8) + scale_color_manual(values = c(`Training QC` = "#0072B2", `Hidden QC` = "#D55E00")) +
    labs(title = "QC positions and holdout status", x = "Within-batch order", y = "Batch", color = NULL) + theme_benchmark()
  randomization_plot <- (p_comp | p_qc) / p_seq + plot_annotation(title = paste(dataset, "randomization diagnostics"))
  ggsave(file.path(figure_dir, "randomization_diagnostics.pdf"), randomization_plot, width = 12, height = 9)

  method_order <- winn_method_order()
  method_palette <- winn_method_palette()
  cache_path <- function(label) file.path(
    result_dir, "method_cache",
    paste0(gsub("[^A-Za-z0-9]+", "_", tolower(label)), "_", dataset, "_hiddenqc_", seed, "_v1.rds")
  )
  method_order <- method_order[vapply(method_order, function(label) file.exists(cache_path(label)), logical(1))]
  if (!"Raw" %in% method_order) stop("Output augmentation requires a valid Raw method cache.")
  matrices <- setNames(lapply(method_order, function(label) readRDS(cache_path(label))), method_order)
  common_effect_features <- Reduce(intersect, lapply(matrices, rownames))

  # Redraw the core comparison panels with human-readable facet labels.
  metrics_for_plot <- read.csv(file.path(result_dir, "method_metrics.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  metrics_for_plot$method <- factor(metrics_for_plot$method, levels = method_order)
  technical_source <- metrics_for_plot |>
    select(method, heldout_qc_cv_mean, residual_gam_deviance_mean, residual_ljung_box_proportion, batch_weighted_pc_r2) |>
    pivot_longer(-method, names_to = "metric", values_to = "value")
  technical_labels <- c(heldout_qc_cv_mean = "Held-out QC CV (%)", residual_gam_deviance_mean = "Residual run-order GAM deviance",
                        residual_ljung_box_proportion = "Residual Ljung-Box proportion", batch_weighted_pc_r2 = "Batch weighted-PC R²")
  write.csv(technical_source, file.path(source_dir, "technical_comparison.csv"), row.names = FALSE, quote = TRUE)
  p_technical <- ggplot(technical_source, aes(method, value, fill = method)) + geom_col() +
    facet_wrap(~metric, scales = "free_y", labeller = as_labeller(technical_labels)) +
    scale_fill_manual(values = method_palette, guide = "none") +
    labs(title = paste(if (dataset == "sacurine") "Sacurine negative-ion" else "WaveICA adenocarcinoma", "technical comparison"), x = NULL, y = NULL) +
    theme_benchmark() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(file.path(figure_dir, "technical_comparison.pdf"), p_technical, width = 13, height = 7)
  biological_source <- metrics_for_plot |>
    select(method, biological_weighted_pc_r2_mean, biological_associated_features, cross_batch_effect_pearson_median) |>
    pivot_longer(-method, names_to = "metric", values_to = "value")
  biological_labels <- c(biological_weighted_pc_r2_mean = "Biological weighted-PC R²", biological_associated_features = "FDR-associated features",
                         cross_batch_effect_pearson_median = "Median cross-batch effect correlation")
  write.csv(biological_source, file.path(source_dir, "biological_comparison.csv"), row.names = FALSE, quote = TRUE)
  p_biological <- ggplot(biological_source, aes(method, value, fill = method)) + geom_col() +
    facet_wrap(~metric, scales = "free_y", labeller = as_labeller(biological_labels)) +
    scale_fill_manual(values = method_palette, guide = "none") +
    labs(title = paste(if (dataset == "sacurine") "Sacurine negative-ion" else "WaveICA adenocarcinoma", "biological preservation"), x = NULL, y = NULL) +
    theme_benchmark() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(file.path(figure_dir, "biological_comparison.pdf"), p_biological, width = 12, height = 6)

  fit_effects <- function(mat, md, dataset, method, batch_id) {
    z <- log1p(mat[, md$sample_id, drop = FALSE])
    if (dataset == "sacurine") {
      z <- log1p(sweep(mat[, md$sample_id, drop = FALSE], 2, md$osmolality, "/"))
      md$age_z <- (md$age - mean(study$age, na.rm = TRUE)) / sd(study$age, na.rm = TRUE)
      md$bmi_z <- (md$bmi - mean(study$bmi, na.rm = TRUE)) / sd(study$bmi, na.rm = TRUE)
      md$gender <- factor(md$gender); md$sampling_factor <- factor(md$sampling)
      design <- model.matrix(~ age_z + bmi_z + gender + sampling_factor, data = md)
      names_map <- c(age = "age_z", bmi = "bmi_z", gender = grep("^gender", colnames(design), value = TRUE)[1])
    } else {
      md$biological_group <- factor(md$biological_group, levels = c("group_0", "group_1"))
      design <- model.matrix(~ biological_group, data = md)
      names_map <- c(biological_group = grep("^biological_group", colnames(design), value = TRUE)[1])
    }
    fit <- limma::eBayes(limma::lmFit(z, design), robust = TRUE)
    bind_rows(lapply(names(names_map), function(variable) {
      i <- match(names_map[[variable]], colnames(design))
      effect <- fit$coefficients[, i]
      outcome_sd <- apply(z, 1, sd)
      data.frame(method = method, batch = as.character(batch_id), variable = variable,
                 feature_id = rownames(z), coefficient = effect,
                 standardized_effect = effect / outcome_sd, p_value = fit$p.value[, i],
                 stringsAsFactors = FALSE)
    }))
  }

  batch_effects <- bind_rows(lapply(method_order, function(method) bind_rows(lapply(unique(study$batch), function(batch_id) {
    md <- study[study$batch == batch_id, , drop = FALSE]
    fit_effects(matrices[[method]][common_effect_features, , drop = FALSE], md, dataset, method, batch_id)
  }))))
  write.csv(batch_effects, file.path(source_dir, paste0(dataset, "_batch_effect_estimates.csv")), row.names = FALSE, quote = TRUE, na = "")

  if (dataset == "sacurine") {
    wide <- batch_effects |> select(method, variable, feature_id, batch, standardized_effect) |>
      pivot_wider(names_from = batch, values_from = standardized_effect)
    batch_names <- unique(study$batch)
    p_effect <- ggplot(wide, aes(x = .data[[batch_names[1]]], y = .data[[batch_names[2]]])) +
      geom_hline(yintercept = 0, color = "grey85") + geom_vline(xintercept = 0, color = "grey85") +
      geom_point(alpha = 0.5, size = 0.7) + geom_abline(slope = 1, intercept = 0, color = "#D55E00", linewidth = 0.45) +
      facet_grid(variable ~ method, scales = "free") +
      labs(title = "Sacurine age, BMI, and gender effect concordance", x = paste(batch_names[1], "standardized effect"), y = paste(batch_names[2], "standardized effect")) + theme_benchmark() +
      theme(axis.text = element_text(size = 6), strip.text.x = element_text(size = 7))
    ggsave(file.path(figure_dir, "sacurine_effect_concordance.pdf"), p_effect, width = 16, height = 8)
  } else {
    primary <- read.csv(file.path(result_dir, "waveica_group_associations_primary.csv"), check.names = FALSE)
    direction_consistency <- bind_rows(lapply(method_order, function(method) {
      significant_ids <- primary$feature_id[primary$method == method & primary$significant]
      effects <- batch_effects[batch_effects$method == method & batch_effects$feature_id %in% significant_ids, , drop = FALSE]
      by_feature <- split(effects$standardized_effect, effects$feature_id)
      consistent <- vapply(by_feature, function(value) length(value) == length(unique(study$batch)) &&
                             all(is.finite(value)) && length(unique(sign(value))) == 1L, logical(1))
      data.frame(method = method, pooled_significant_features = length(significant_ids),
                 evaluable_in_every_batch = length(consistent),
                 consistent_direction_every_batch = sum(consistent),
                 proportion_consistent_direction_every_batch = if (length(consistent)) mean(consistent) else NA_real_,
                 stringsAsFactors = FALSE)
    }))
    write.csv(direction_consistency, file.path(result_dir, "waveica_cross_batch_direction_consistency.csv"), row.names = FALSE, quote = TRUE, na = "")
    metrics_path <- file.path(result_dir, "method_metrics.csv")
    metrics <- read.csv(metrics_path, check.names = FALSE, stringsAsFactors = FALSE)
    metrics$cross_batch_pooled_significant_direction_consistency <- direction_consistency$proportion_consistent_direction_every_batch[match(metrics$method, direction_consistency$method)]
    write.csv(metrics, metrics_path, row.names = FALSE, quote = TRUE, na = "")
    metrics_long <- tidyr::pivot_longer(metrics, -c(method, repeat_metric_reason), names_to = "metric", values_to = "value")
    write.csv(metrics_long, file.path(result_dir, "method_metrics_long.csv"), row.names = FALSE, quote = TRUE, na = "")
    raw_primary <- primary[primary$method == "Raw", , drop = FALSE]
    panel_features <- head(raw_primary$feature_id[order(-abs(raw_primary$standardized_effect), raw_primary$feature_id)], 12L)
    panel <- batch_effects[batch_effects$feature_id %in% panel_features, , drop = FALSE]
    panel$feature_id <- factor(panel$feature_id, levels = rev(panel_features))
    write.csv(panel, file.path(source_dir, "waveica_batch_effect_panel.csv"), row.names = FALSE, quote = TRUE, na = "")
    p_effect <- ggplot(panel, aes(standardized_effect, feature_id, color = batch)) +
      geom_vline(xintercept = 0, color = "grey80") + geom_point(position = position_dodge(width = 0.55), size = 1.2) +
      facet_wrap(~method, ncol = 3) + labs(title = "WaveICA group effects across supplied batches",
                                           subtitle = "Panel fixed by the 12 largest absolute Raw primary standardized effects",
                                           x = "Within-batch standardized group effect", y = NULL, color = "Batch") + theme_benchmark()
    ggsave(file.path(figure_dir, "waveica_batch_effects.pdf"), p_effect, width = 12, height = 11)
  }

  # Label-invariant alignment between each WiNN segmentation and supplied batches.
  assignments_path <- file.path(result_dir, "winn_auto_batch_assignments.csv")
  selectivity_path <- file.path(result_dir, "winn_selectivity_summary.csv")
  if (file.exists(assignments_path) && file.exists(selectivity_path)) {
    assignments <- read.csv(assignments_path, check.names = FALSE, stringsAsFactors = FALSE)
    adjusted_rand <- function(a, b) {
      tab <- table(a, b); n <- sum(tab)
      choose2 <- function(v) v * (v - 1) / 2
      index <- sum(choose2(tab)); row_sum <- sum(choose2(rowSums(tab))); col_sum <- sum(choose2(colSums(tab)))
      expected <- row_sum * col_sum / choose2(n)
      denominator <- 0.5 * (row_sum + col_sum) - expected
      if (denominator == 0) as.numeric(identical(as.character(a), as.character(b))) else (index - expected) / denominator
    }
    alignment <- bind_rows(lapply(split(assignments, assignments$method), function(d) data.frame(
      method = d$method[1], supplied_segments = length(unique(d$supplied_batch)), inferred_segments = length(unique(d$detected_batch)),
      adjusted_rand_index = adjusted_rand(d$supplied_batch, d$detected_batch),
      exact_boundary_recall = {
        supplied_cp <- which(d$supplied_batch[-1] != d$supplied_batch[-nrow(d)])
        detected_cp <- which(d$detected_batch[-1] != d$detected_batch[-nrow(d)])
        if (!length(supplied_cp)) NA_real_ else mean(supplied_cp %in% detected_cp)
      }, stringsAsFactors = FALSE
    )))
    write.csv(alignment, file.path(result_dir, "winn_batch_alignment.csv"), row.names = FALSE, quote = TRUE, na = "")
    selectivity <- read.csv(selectivity_path, check.names = FALSE, stringsAsFactors = FALSE)
    for (column in setdiff(names(alignment), "method")) selectivity[[column]] <- alignment[[column]][match(selectivity$method, alignment$method)]
    write.csv(selectivity, selectivity_path, row.names = FALSE, quote = TRUE, na = "")
    selectivity_long <- bind_rows(
      data.frame(method = selectivity$method, gate = "Drift corrected", proportion = selectivity$proportion_unique_features_detrended),
      data.frame(method = selectivity$method, gate = "Batch corrected", proportion = selectivity$proportion_features_selected_for_batch),
      data.frame(method = selectivity$method, gate = "Both gates", proportion = selectivity$drift_batch_overlap_features / selectivity$features_tested_for_batch)
    )
    write.csv(selectivity_long, file.path(source_dir, "winn_selectivity.csv"), row.names = FALSE, quote = TRUE)
    p_selectivity <- ggplot(selectivity_long, aes(method, proportion, fill = gate)) +
      geom_col(position = "dodge") + scale_y_continuous(labels = scales::label_percent()) +
      labs(title = paste(if (dataset == "sacurine") "Sacurine negative-ion" else "WaveICA adenocarcinoma", "WiNN selectivity"),
           x = NULL, y = "Proportion of input features", fill = NULL) + theme_benchmark() +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))
    ggsave(file.path(figure_dir, "winn_selectivity.pdf"), p_selectivity, width = 9, height = 5)
  }
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
  root <- dirname(dirname(script_path))
  dataset_arg <- grep("^--dataset=", commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(dataset_arg) != 1L) stop("Supply --dataset=sacurine or --dataset=waveica_adenocarcinoma.")
  augment_human_public_benchmark_outputs(root, sub("^--dataset=", "", dataset_arg))
}
