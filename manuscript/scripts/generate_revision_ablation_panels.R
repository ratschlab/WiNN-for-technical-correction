#!/usr/bin/env Rscript
# Supplementary Figures S1 (ordered stage impacts), S2 (selective vs forced gates)
# and S9 (categorical batch R2 trajectory). Ported unchanged in plotting logic from
# WiNN_evaluation_package_2026-09-14/generate_manuscript_panels.R so these panels have
# a committed generator. No correction, tuning or model selection is run here.
#
# Usage:
#   Rscript scripts/generate_revision_ablation_panels.R \
#     --final-root=<release results/final> \
#     --source-root=<dir with SF_ablation_stage_impacts.csv, SF_categorical_batch_r2_trajectory.csv>
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(readr); library(scales)
})
args <- commandArgs(TRUE)
arg <- function(name) {
  hit <- sub(paste0("^--", name, "="), "", grep(paste0("^--", name, "="), args, value = TRUE))
  if (length(hit) != 1L) stop("Supply --", name, "=", call. = FALSE)
  normalizePath(hit, mustWork = TRUE)
}
final_root <- arg("final-root")
source_root <- arg("source-root")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
manuscript_root <- dirname(dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg), fixed = TRUE))))
supp_dir <- file.path(manuscript_root, "figures", "supplementary_panels")
panel_source_dir <- file.path(manuscript_root, "figures", "source_data")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_source_dir, recursive = TRUE, showWarnings = FALSE)

pt_to_mm <- 25.4 / 72.27
line_width <- 0.45 * pt_to_mm
point_stroke <- 0.35 * pt_to_mm
gate_colours <- c(
  "Selective / selective" = "#0072B2", "Forced drift" = "#D55E00",
  "Forced batch" = "#009E73", "Forced drift + batch" = "#7A5195"
)
dataset_labels <- c(
  simulation = "Simulation", mtbls79 = "MTBLS79",
  clinical_fiams = "Clinical FIA-MS", batchcorr_set1 = "BatchCorr Set 1",
  sacurine = "Sacurine", waveica = "WaveICA"
)
theme_panel <- function(base_size = 9.5) {
  base_size <- max(base_size, 9.5)
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(colour = "#596273", linewidth = line_width),
      axis.ticks = element_line(colour = "#596273", linewidth = line_width),
      panel.grid.major = element_line(colour = "#E2E7ED", linewidth = line_width),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(face = "bold", colour = "#172033", size = base_size),
      plot.subtitle = element_text(colour = "#465269", size = base_size - 0.4),
      axis.title = element_text(colour = "#263247", size = base_size - 0.2),
      axis.text = element_text(colour = "#263247", size = base_size - 0.4),
      strip.background = element_rect(fill = "#F0F2F5", colour = NA),
      strip.text = element_text(face = "bold", colour = "#263247", size = base_size - 0.2),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5),
      plot.margin = margin(3, 4, 3, 4)
    )
}
save_panel <- function(plot, stem, width_mm, height_mm, source_data) {
  ggsave(file.path(supp_dir, paste0(stem, ".pdf")), plot, width = width_mm, height = height_mm,
         units = "mm", device = grDevices::pdf, useDingbats = FALSE, bg = "white")
  write_csv(as.data.frame(source_data), file.path(panel_source_dir, paste0(stem, ".csv")), na = "")
  cat("wrote", stem, "\n")
}

# S1: ordered stage impacts ------------------------------------------------------
stage_data <- read_csv(file.path(source_root, "SF_ablation_stage_impacts.csv"), show_col_types = FALSE) |>
  mutate(
    stage_short = factor(recode(stage, Outlier_shrinkage = "MAD", Selective_drift = "Drift",
                                Selective_batch = "Batch", PQN = "PQN"),
                         levels = c("MAD", "Drift", "Batch", "PQN")),
    dataset = factor(dataset, levels = unname(dataset_labels)),
    benefit_native = direction_sign * raw_change,
    benefit_display = if_else(
      metric %in% c("residual_run_order_gam_mean_deviance", "batch_weighted_pc_r2_categorical"),
      100 * benefit_native, benefit_native)
  )
for (spec in list(
  list("A", "heldout_qc_cv_mean", "QC CV improvement", "Percentage points"),
  list("B", "residual_run_order_gam_mean_deviance", "Run-order improvement", "Percentage points"),
  list("C", "batch_weighted_pc_r2_categorical", "Batch-structure improvement", "Percentage points"))) {
  panel_data <- stage_data |> filter(metric == spec[[2]])
  limit <- max(abs(panel_data$benefit_display), na.rm = TRUE)
  plot <- ggplot(panel_data, aes(stage_short, dataset, fill = benefit_display)) +
    geom_tile(colour = "white", linewidth = line_width) +
    scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                         midpoint = 0, limits = c(-limit, limit), oob = squish) +
    labs(title = paste(spec[[1]], spec[[3]]), x = NULL, y = NULL, fill = spec[[4]]) +
    theme_panel(9.2) +
    theme(axis.text.x = element_text(angle = 28, hjust = 1), axis.line = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom")
  save_panel(plot, paste0("S1", spec[[1]], "_stage_", spec[[2]]), 85, 55, panel_data)
}

# S2: selective and forced gates ------------------------------------------------------
gate_data <- read_csv(file.path(final_root, "tables", "selective_vs_forced_metrics_long.csv"),
                      show_col_types = FALSE) |>
  mutate(
    gate_label = recode(method_id, G_SS = "Selective / selective", G_AS = "Forced drift",
                        G_SA = "Forced batch", G_AA = "Forced drift + batch"),
    gate_label = factor(gate_label, levels = names(gate_colours)),
    dataset = factor(dataset, levels = unname(dataset_labels)),
    dataset_short = recode(as.character(dataset), `Clinical FIA-MS` = "Clinical",
                           `BatchCorr Set 1` = "BatchCorr"),
    dataset_short = factor(dataset_short, levels = c("Simulation", "MTBLS79", "Clinical",
                                                     "BatchCorr", "Sacurine", "WaveICA"))
  )
for (spec in list(
  list("A", "heldout_qc_cv_mean", "Held-out QC CV", "Held-out QC CV (%)"),
  list("B", "residual_run_order_gam_mean_deviance", "Residual run-order deviance", "GAM deviance explained"),
  list("C", "batch_weighted_pc_r2_categorical", "Residual batch structure", "Weighted-PC R2"))) {
  panel_data <- gate_data |> filter(metric == spec[[2]])
  plot <- ggplot(panel_data, aes(gate_label, value, colour = gate_label)) +
    geom_point(size = 1.8, shape = 16, stroke = 0) +
    facet_wrap(vars(dataset_short), nrow = 1, scales = "free_y") +
    scale_colour_manual(values = gate_colours, drop = FALSE) +
    labs(title = paste(spec[[1]], spec[[3]]), x = NULL, y = spec[[4]]) +
    theme_panel(9.2) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
  save_panel(plot, paste0("S2", spec[[1]], "_", spec[[2]]), 174, 55, panel_data)
}
seed_gate <- read_csv(file.path(final_root, "tables", "simulation_seed_ablation_metrics.csv"),
                      show_col_types = FALSE) |>
  filter(method_id %in% c("G_SS", "G_AS", "G_SA", "G_AA")) |>
  mutate(
    gate_label = recode(method_id, G_SS = "Selective / selective", G_AS = "Forced drift",
                        G_SA = "Forced batch", G_AA = "Forced drift + batch"),
    gate_label = factor(gate_label, levels = rev(names(gate_colours)))
  )
for (spec in list(
  list("D", "truth_log1p_rmse", "Simulation clean-truth RMSE", "Clean-truth RMSE"),
  list("E", "truth_metabolite_profile_icc_mean", "Simulation feature ICC(A,1)", "ICC(A,1) to clean truth"),
  list("F", "truth_sample_profile_pearson_mean", "Simulation sample correlation", "Pearson r to clean truth"))) {
  panel_data <- seed_gate |> filter(metric == spec[[2]]) |> group_by(method_id, gate_label) |>
    mutate(method_mean = mean(value), method_se = sd(value) / sqrt(n())) |> ungroup()
  summary_data <- panel_data |> distinct(method_id, gate_label, method_mean, method_se)
  plot <- ggplot(panel_data, aes(value, gate_label, colour = gate_label)) +
    geom_segment(data = summary_data,
                 aes(x = method_mean - method_se, xend = method_mean + method_se,
                     y = gate_label, yend = gate_label),
                 linewidth = line_width, colour = "#172033") +
    geom_point(position = position_jitter(width = 0, height = 0.08, seed = 2026),
               shape = 21, fill = "white", size = 1.3, stroke = point_stroke, alpha = 0.72) +
    geom_point(data = summary_data, aes(method_mean, gate_label), inherit.aes = FALSE,
               shape = 23, fill = "white", colour = "#172033", size = 2.0, stroke = point_stroke) +
    scale_colour_manual(values = gate_colours, drop = FALSE) +
    labs(title = paste(spec[[1]], spec[[3]]), x = spec[[4]], y = NULL) +
    theme_panel(9.2) +
    theme(legend.position = "none", panel.grid.major.y = element_blank())
  save_panel(plot, paste0("S2", spec[[1]], "_", spec[[2]]), 85, 55, panel_data)
}

# S9: categorical batch trajectory ------------------------------------------------------
trajectory <- read_csv(file.path(source_root, "SF_categorical_batch_r2_trajectory.csv"),
                       show_col_types = FALSE) |>
  mutate(dataset = factor(dataset, levels = unname(dataset_labels)),
         stage = factor(stage, levels = c("Raw", "Outlier shrinkage", "Selective drift",
                                          "Selective batch", "PQN")))
for (spec in list(list("A", "simulation", "Simulation"), list("B", "mtbls79", "MTBLS79"),
                  list("C", "clinical_fiams", "Clinical FIA-MS"),
                  list("D", "batchcorr_set1", "BatchCorr Set 1"),
                  list("E", "sacurine", "Sacurine"), list("F", "waveica", "WaveICA"))) {
  panel_data <- trajectory |> filter(dataset_key == spec[[2]])
  plot <- ggplot(panel_data, aes(stage, value, group = 1)) +
    geom_line(linewidth = line_width, colour = "#0072B2") +
    geom_point(size = 1.6, colour = "#0072B2", shape = 16, stroke = 0) +
    scale_x_discrete(labels = c("Raw" = "Raw", "Outlier shrinkage" = "MAD",
                                "Selective drift" = "Drift", "Selective batch" = "Batch",
                                "PQN" = "PQN")) +
    labs(title = paste(spec[[1]], spec[[3]]), x = NULL, y = "Batch weighted-PC R²") +
    theme_panel(9.2) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  save_panel(plot, paste0("S9", spec[[1]], "_batch_r2_", spec[[2]]), 85, 50, panel_data)
}
