#!/usr/bin/env Rscript

# This script is intentionally downstream-only. It reads saved cumulative metric
# tables and cannot invoke a correction method, tuning routine, or selector.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
source(file.path(repo_root, "analysis", "ablation_stage_helpers.R"), local = TRUE)

required_packages <- c("ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)

table_dir <- file.path(repo_root, "results", "final", "tables")
figure_dir <- file.path(repo_root, "results", "final", "figures")
source_dir <- file.path(figure_dir, "source_data")
source_paths <- list.files(table_dir, pattern = "^ablation_cumulative_values_.*[.]csv$", full.names = TRUE)
if (length(source_paths) != 6L) {
  stop("Expected six cumulative ablation tables; found ", length(source_paths), ".", call. = FALSE)
}

cumulative <- do.call(rbind, lapply(source_paths, function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}))
rownames(cumulative) <- NULL
cumulative <- normalize_ablation_cumulative_metric_names(cumulative)
impacts <- derive_ablation_stage_impacts(cumulative)
validate_ablation_stage_impacts(cumulative, impacts)

write.csv(impacts, file.path(table_dir, "ablation_stage_impacts.csv"), row.names = FALSE, na = "")
write.csv(impacts, file.path(table_dir, "ablation_stage_impacts_table.csv"), row.names = FALSE, na = "")
for (table_id in unique(impacts$table_id)) {
  write.csv(
    impacts[impacts$table_id == table_id, , drop = FALSE],
    file.path(table_dir, paste0("ablation_step_impacts_", table_id, ".csv")),
    row.names = FALSE, na = ""
  )
}
write.csv(impacts, file.path(source_dir, "ablation_step_impacts.csv"), row.names = FALSE, na = "")

dataset_levels <- c(
  "Simulation", "MTBLS79", "Clinical FIA-MS", "BatchCorr Set 1", "Sacurine", "WaveICA"
)
stage_levels <- ablation_stage_spec()$stage
plot_source <- impacts[
  impacts$table_id %in% c("heldout_qc_cv", "gam_deviance", "batch_r2"),
  , drop = FALSE
]
plot_source$metric_label <- factor(
  plot_source$table_id,
  levels = c("heldout_qc_cv", "gam_deviance", "batch_r2"),
  labels = c("Held-out QC CV", "Residual GAM deviance", "Batch weighted-PC R²")
)
plot_source$dataset <- factor(plot_source$dataset, levels = dataset_levels)
plot_source$stage <- factor(plot_source$stage, levels = stage_levels)
plot_source$stage_label <- factor(
  plot_source$stage_label,
  levels = ablation_stage_spec()$stage_label
)
plot_source$within_metric_score <- NA_real_
for (metric_value in levels(plot_source$metric_label)) {
  selected <- plot_source$metric_label == metric_value & is.finite(plot_source$direction_normalized_benefit)
  maximum <- max(abs(plot_source$direction_normalized_benefit[selected]), na.rm = TRUE)
  if (is.finite(maximum) && maximum > 0) {
    plot_source$within_metric_score[selected] <-
      100 * plot_source$direction_normalized_benefit[selected] / maximum
  }
}

write.csv(plot_source, file.path(source_dir, "ablation_stage_impacts.csv"), row.names = FALSE, na = "")

theme_publication <- function(base_size = 9) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "#343A40"),
      axis.title = ggplot2::element_text(colour = "#202428", face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#EEF1F4", colour = NA),
      strip.text = ggplot2::element_text(colour = "#202428", face = "bold"),
      plot.title = ggplot2::element_text(colour = "#202428", face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "#4F5963"),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold")
    )
}

impact_plot <- ggplot2::ggplot(plot_source, ggplot2::aes(stage_label, dataset, fill = within_metric_score)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
  ggplot2::facet_wrap(~metric_label, scales = "free", ncol = 1) +
  ggplot2::scale_fill_gradient2(
    low = "#B04A5A", mid = "#F3F0EA", high = "#2B6F77", midpoint = 0,
    limits = c(-100, 100), oob = scales::squish
  ) +
  ggplot2::labs(
    title = "Incremental contribution of WiNN stages",
    subtitle = "Changes are consecutive; positive colors indicate improvement and are scaled within each metric",
    x = NULL, y = NULL, fill = "Direction-normalized\nwithin-metric score"
  ) +
  theme_publication(9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

ggplot2::ggsave(
  file.path(figure_dir, "ablation_stage_impacts.pdf"), impact_plot,
  width = 8, height = 8, units = "in"
)
ggplot2::ggsave(
  file.path(figure_dir, "ablation_stage_impacts.png"), impact_plot,
  width = 8, height = 8, units = "in", dpi = 320
)

figure_manifest_path <- file.path(figure_dir, "figure_manifest.csv")
figure_manifest <- read.csv(figure_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
figure_manifest$source_csv[figure_manifest$figure == "ablation_stage_impacts"] <-
  "figures/source_data/ablation_stage_impacts.csv"
write.csv(figure_manifest, figure_manifest_path, row.names = FALSE, na = "")

cat("Repaired ablation derivative tables and figure without executing a correction method.\n")
