#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript or R --file.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly = TRUE)
root_arg <- grep("^--root=", args, value = TRUE)

required <- c("dplyr", "ggplot2", "tidyr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing package(s): ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})
source(file.path(repo_root, "scripts", "plot_theme.R"))

datasets <- c("simulation", "mtbls79", "batchcorr_set1", "sacurine", "waveica")
dataset_labels <- c(
  simulation = "Simulation", mtbls79 = "MTBLS79",
  batchcorr_set1 = "BatchCorr Set 1", sacurine = "Sacurine", waveica = "WaveICA"
)
root <- if (length(root_arg) == 1L) {
  normalizePath(sub("^--root=", "", root_arg), mustWork = TRUE)
} else {
  file.path(repo_root, "results", "winn_ablations")
}
combined_dir <- file.path(root, "combined")
figure_dir <- file.path(combined_dir, "figures")
source_dir <- file.path(combined_dir, "figure_source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

read_dataset_table <- function(dataset, relative) {
  path <- file.path(root, dataset, relative)
  if (!file.exists(path)) stop("Required completed output is missing: ", path)
  value <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"dataset" %in% names(value)) value$dataset <- dataset
  value
}

bind_dataset_table <- function(relative) {
  bind_rows(lapply(datasets, read_dataset_table, relative = relative))
}

write_table <- function(value, name, directory = combined_dir) {
  path <- file.path(directory, name)
  write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
  invisible(path)
}

write_svg_from_pdf <- function(pdf_path, svg_path) {
  converter <- Sys.which("pdftocairo")
  if (!nzchar(converter)) {
    warning("pdftocairo is unavailable; SVG was not generated for ", basename(pdf_path))
    return(invisible(FALSE))
  }
  status <- system2(converter, c("-svg", shQuote(pdf_path), shQuote(svg_path)))
  if (!identical(status, 0L) || !file.exists(svg_path) || file.info(svg_path)$size <= 0) {
    warning("Failed to convert ", basename(pdf_path), " to SVG.")
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

metrics <- bind_dataset_table(file.path("metrics", "primary_metrics.csv")) |>
  mutate(
    dataset_label = unname(dataset_labels[dataset]),
    value = as.numeric(value),
    denominator = as.character(denominator)
  )
selection <- bind_dataset_table(file.path("metrics", "selectivity_counts.csv"))
magnitude <- bind_dataset_table(file.path("metrics", "correction_magnitude.csv"))
runtime <- bind_dataset_table(file.path("metrics", "runtime_summary.csv"))
validation <- bind_dataset_table(file.path("config", "validation_checks.csv"))
equality <- bind_dataset_table(file.path("config", "fixed_result_equivalence.csv"))

cumulative <- metrics |> filter(variant_id %in% paste0("C", 0:4, c("_RAW", "_OUTLIER", "_SELECTIVE_DRIFT", "_SELECTIVE_BATCH", "_FULL_FIXED")))
gate <- metrics |> filter(variant_id %in% c("G_SS", "G_AS", "G_SA", "G_AA"))
write_table(cumulative, "cumulative_primary_metrics.csv")
write_table(gate, "gate_primary_metrics.csv")
write_table(selection, "selectivity_counts.csv")
write_table(magnitude, "correction_magnitude.csv")
write_table(runtime, "runtime_summary.csv")
write_table(validation, "validation_checks.csv")
write_table(equality, "fixed_result_equivalence.csv")

cumulative_order <- c("C0_RAW", "C1_OUTLIER", "C2_SELECTIVE_DRIFT", "C3_SELECTIVE_BATCH", "C4_FULL_FIXED")
increments <- cumulative |>
  mutate(variant_id = factor(variant_id, levels = cumulative_order)) |>
  arrange(dataset, metric, variant_id) |>
  group_by(dataset, evaluation_panel, metric, metric_direction, units) |>
  mutate(
    previous_variant = lag(as.character(variant_id)),
    previous_value = lag(value),
    next_variant = as.character(variant_id),
    next_value = value
  ) |>
  filter(!is.na(previous_variant)) |>
  transmute(
    dataset, evaluation_panel, metric, metric_direction,
    previous_variant, next_variant, previous_value, next_value,
    absolute_delta = next_value - previous_value,
    relative_delta = ifelse(is.finite(previous_value) & abs(previous_value) > .Machine$double.eps,
      100 * (next_value - previous_value) / abs(previous_value), NA_real_),
    denominator_notes = ifelse(is.finite(previous_value) & abs(previous_value) > .Machine$double.eps,
      "Relative delta is percent of absolute previous value.",
      "Relative delta omitted because the previous value is zero or nonfinite."),
    units
  ) |>
  ungroup()
write_table(increments, "cumulative_incremental_effects.csv")

gate_wide <- gate |>
  select(dataset, metric, metric_direction, units, variant_id, value) |>
  distinct() |>
  pivot_wider(names_from = variant_id, values_from = value)
required_gate <- c("G_SS", "G_AS", "G_SA", "G_AA")
if (!all(required_gate %in% names(gate_wide))) stop("Gate table is incomplete.")
gate_contrasts <- gate_wide |>
  rowwise() |>
  do({
    row <- .
    values <- c(
      drift_gate_forced_all = mean(c(row$G_AS, row$G_AA)) - mean(c(row$G_SS, row$G_SA)),
      batch_gate_forced_all = mean(c(row$G_SA, row$G_AA)) - mean(c(row$G_SS, row$G_AS)),
      drift_batch_interaction = row$G_AA - row$G_AS - row$G_SA + row$G_SS,
      G_AS_minus_G_SS = row$G_AS - row$G_SS,
      G_SA_minus_G_SS = row$G_SA - row$G_SS,
      G_AA_minus_G_SS = row$G_AA - row$G_SS
    )
    data.frame(
      dataset = row$dataset, metric = row$metric, metric_direction = row$metric_direction,
      units = row$units, contrast = names(values), value = as.numeric(values),
      interpretation_note = "Descriptive contrast only; no feature-level pseudo-replication or hypothesis test.",
      stringsAsFactors = FALSE
    )
  }) |>
  ungroup()
write_table(gate_contrasts, "gate_factorial_contrasts.csv")

preservation_metric <- c(
  simulation = "truth_sample_profile_pearson_mean",
  mtbls79 = "sample_replicate_pearson_median",
  batchcorr_set1 = "genuine_replicate_pearson_median",
  sacurine = "cross_batch_effect_pearson_median",
  waveica = "group_weighted_pc_r2"
)
main_metrics <- c(
  "heldout_qc_cv_mean", "hidden_qc_profile_pearson", "residual_gam_deviance_mean",
  "residual_ljung_box_proportion", "batch_weighted_pc_r2", unique(preservation_metric)
)
publication_main <- metrics |>
  filter(metric %in% main_metrics) |>
  filter(metric %in% c(
    "heldout_qc_cv_mean", "hidden_qc_profile_pearson", "residual_gam_deviance_mean",
    "residual_ljung_box_proportion", "batch_weighted_pc_r2"
  ) | metric == unname(preservation_metric[dataset])) |>
  select(dataset, dataset_label, variant_id, variant_label, analysis_type, metric, value,
    metric_direction, denominator, units, notes)
write_table(publication_main, "publication_main_table.csv")
write_table(metrics, "publication_supplementary_table.csv")

stage_labels <- c(
  C0_RAW = "Raw", C1_OUTLIER = "Outlier", C2_SELECTIVE_DRIFT = "+ selective drift",
  C3_SELECTIVE_BATCH = "+ selective batch", C4_FULL_FIXED = "+ PQN"
)
metric_labels <- c(
  residual_gam_deviance_mean = "Residual GAM deviance",
  batch_weighted_pc_r2 = "Batch/plate weighted-PC R²",
  heldout_qc_cv_mean = "Held-out QC/reference CV (%)"
)

figure_a <- cumulative |>
  filter(metric %in% names(metric_labels) | metric == unname(preservation_metric[dataset])) |>
  mutate(
    stage = factor(variant_id, levels = cumulative_order, labels = unname(stage_labels[cumulative_order])),
    panel = case_when(
      metric %in% names(metric_labels) ~ unname(metric_labels[metric]),
      TRUE ~ "Design-specific preservation"
    ),
    preservation_definition = ifelse(panel == "Design-specific preservation", metric, NA_character_)
  )
write_table(figure_a, "figure_A_cumulative_stage_effects.csv", source_dir)
p_a <- ggplot(figure_a, aes(stage, value, group = dataset, colour = dataset_label)) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.6) +
  facet_wrap(~panel, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Metric value", colour = "Dataset",
    title = "Cumulative WiNN stage effects",
    subtitle = "Preservation endpoints are design-specific and are identified in the source table") +
  scale_colour_manual(values = winn_dataset_palette()) +
  winn_theme_publication(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
ggsave(file.path(figure_dir, "Figure_A_cumulative_stage_effects.pdf"), p_a, width = 10, height = 7)
ggsave(file.path(figure_dir, "Figure_A_cumulative_stage_effects.png"), p_a, width = 10, height = 7, dpi = 300)
write_svg_from_pdf(
  file.path(figure_dir, "Figure_A_cumulative_stage_effects.pdf"),
  file.path(figure_dir, "Figure_A_cumulative_stage_effects.svg")
)

gate_labels <- c(G_SS = "Selective/Selective", G_AS = "All/Selective", G_SA = "Selective/All", G_AA = "All/All")
figure_b_metrics <- gate |>
  filter(metric %in% c(
    "residual_gam_deviance_mean", "batch_weighted_pc_r2", "hidden_qc_profile_pearson",
    character()
  ) | metric == unname(preservation_metric[dataset]))
figure_b_magnitude <- magnitude |>
  filter(variant_id %in% names(gate_labels)) |>
  transmute(
    dataset, dataset_label = unname(dataset_labels[dataset]), variant_id,
    metric = "median_abs_log1p_change_raw", value = median_abs_log1p_change_raw,
    metric_direction = "descriptive", units = "log1p intensity"
  )
figure_b <- bind_rows(figure_b_metrics, figure_b_magnitude) |>
  mutate(
    gate_combination = factor(variant_id, levels = names(gate_labels), labels = unname(gate_labels)),
    panel = case_when(
      metric == "residual_gam_deviance_mean" ~ "Residual GAM deviance",
      metric == "batch_weighted_pc_r2" ~ "Batch/plate weighted-PC R²",
      metric == "hidden_qc_profile_pearson" ~ "Held-out QC profile correlation",
      metric == "median_abs_log1p_change_raw" ~ "Median absolute correction from Raw",
      TRUE ~ "Design-specific preservation"
    )
  )
write_table(figure_b, "figure_B_gate_factorial.csv", source_dir)
p_b <- ggplot(figure_b, aes(gate_combination, value, group = dataset, colour = dataset_label)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.6) +
  facet_wrap(~panel, scales = "free_y", ncol = 2) +
  labs(x = "Drift gate / batch gate", y = "Metric value", colour = "Dataset",
    title = "Selective and forced-all WiNN gates") +
  scale_colour_manual(values = winn_dataset_palette()) +
  winn_theme_publication(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
ggsave(file.path(figure_dir, "Figure_B_gate_factorial.pdf"), p_b, width = 10, height = 7)
ggsave(file.path(figure_dir, "Figure_B_gate_factorial.png"), p_b, width = 10, height = 7, dpi = 300)
write_svg_from_pdf(
  file.path(figure_dir, "Figure_B_gate_factorial.pdf"),
  file.path(figure_dir, "Figure_B_gate_factorial.svg")
)

figure_c <- selection |>
  filter(variant_id == "G_SS") |>
  select(dataset, drift_corrected_profile_proportion, drift_unique_feature_proportion,
    batch_corrected_feature_proportion, pqn_sample_proportion) |>
  pivot_longer(-dataset, names_to = "selection_measure", values_to = "proportion") |>
  mutate(
    dataset_label = unname(dataset_labels[dataset]),
    selection_label = factor(
      selection_measure,
      levels = c(
        "drift_corrected_profile_proportion",
        "drift_unique_feature_proportion",
        "batch_corrected_feature_proportion",
        "pqn_sample_proportion"
      ),
      labels = c(
        "Drift-corrected feature-batch profiles",
        "Features corrected for drift",
        "Features corrected for residual batch",
        "Samples altered by PQN shrinkage"
      )
    )
  )
write_table(figure_c, "figure_C_selectivity_intervention.csv", source_dir)
p_c <- ggplot(figure_c, aes(dataset_label, proportion, fill = selection_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  scale_fill_manual(values = winn_selection_palette()) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Proportion", fill = "Denominator-specific measure",
    title = "Selective WiNN intervention size",
    subtitle = "Drift profiles, unique features, batch-eligible features, and samples use distinct denominators") +
  winn_theme_publication(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
ggsave(file.path(figure_dir, "Figure_C_selectivity_intervention.pdf"), p_c, width = 10, height = 5.8)
ggsave(file.path(figure_dir, "Figure_C_selectivity_intervention.png"), p_c, width = 10, height = 5.8, dpi = 300)
write_svg_from_pdf(
  file.path(figure_dir, "Figure_C_selectivity_intervention.pdf"),
  file.path(figure_dir, "Figure_C_selectivity_intervention.svg")
)

figure_d <- gate |>
  filter(metric == "batch_weighted_pc_r2" | metric == unname(preservation_metric[dataset])) |>
  select(dataset, dataset_label, variant_id, metric, value) |>
  mutate(kind = ifelse(metric == "batch_weighted_pc_r2", "technical", "preservation")) |>
  select(-metric) |>
  pivot_wider(names_from = kind, values_from = value) |>
  mutate(gate_combination = unname(gate_labels[variant_id]), preservation_metric = preservation_metric[dataset])
write_table(figure_d, "figure_D_correction_preservation_tradeoff.csv", source_dir)
p_d <- ggplot(figure_d, aes(technical, preservation, colour = gate_combination, label = variant_id)) +
  geom_point(size = 2) + geom_text(nudge_y = 0.01, size = 2.3, check_overlap = TRUE) +
  scale_colour_manual(values = winn_ablation_gate_palette()) +
  facet_wrap(~dataset_label, scales = "free", ncol = 3) +
  labs(x = "Residual batch/plate weighted-PC R² (lower is better)",
    y = "Design-specific preservation metric", colour = "Drift / batch gate",
    title = "Technical removal and preservation are separate endpoints") +
  winn_theme_publication(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(figure_dir, "Figure_D_correction_preservation_tradeoff.pdf"), p_d, width = 11, height = 7)
ggsave(file.path(figure_dir, "Figure_D_correction_preservation_tradeoff.png"), p_d, width = 11, height = 7, dpi = 300)
write_svg_from_pdf(
  file.path(figure_dir, "Figure_D_correction_preservation_tradeoff.pdf"),
  file.path(figure_dir, "Figure_D_correction_preservation_tradeoff.svg")
)

figure_source <- bind_rows(
  figure_a |> transmute(figure = "A", dataset, variant_id, metric, value, panel),
  figure_b |> transmute(figure = "B", dataset, variant_id, metric, value, panel),
  figure_c |> transmute(figure = "C", dataset, variant_id = "G_SS", metric = selection_measure, value = proportion, panel = "selectivity"),
  figure_d |> transmute(figure = "D", dataset, variant_id, metric = preservation_metric, value = preservation, panel = "preservation")
)
write_table(figure_source, "figure_source_data.csv")

message("Combined WiNN ablation tables and figures written to: ", combined_dir)
