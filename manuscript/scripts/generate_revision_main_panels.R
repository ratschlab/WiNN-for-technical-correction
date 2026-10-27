#!/usr/bin/env Rscript

# Regenerate revised main panels and S5, from included frozen inputs.
# Corrected QC-RFSC simulation endpoints are integrated from validated Euler outputs.
# No correction, tuning, or model selection is run by this script.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript or R --file.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
package_root <- dirname(dirname(script_path))
project_root <- dirname(package_root)
revision_inputs <- file.path(package_root, "generated", "revision_panel_inputs")

required <- c("dplyr", "ggplot2", "ggrepel", "readr", "scales", "tidyr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(scales)
  library(tidyr)
})

final_root <- revision_inputs
source_root <- file.path(package_root, "generated", "source_data")
main_dir <- file.path(package_root, "figures", "main_panels")
supp_dir <- file.path(package_root, "figures", "supplementary_panels")
panel_source_dir <- file.path(package_root, "figures", "source_data")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_source_dir, recursive = TRUE, showWarnings = FALSE)
pt_to_mm <- 25.4 / 72.27
line_width <- 0.45 * pt_to_mm
point_stroke <- 0.35 * pt_to_mm

main_method_ids <- c(
  "raw", "combat", "qc_rlsc", "tiger", "serrf",
  "winn_fixed_no_qc"
)
simulation_method_ids <- c("raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf", "winn_fixed_no_qc")
all_method_ids <- c(
  "raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf",
  "winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc"
)
method_labels <- c(
  raw = "Raw", combat = "ComBat", qc_rlsc = "QC-RLSC",
  qc_rfsc = "QC-RFSC", tiger = "TIGER", serrf = "SERRF",
  winn_auto_qc = "WiNN (QC-tuned)",
  winn_auto_batch_qc = "WiNN (inferred segments)",
  winn_fixed_no_qc = "WiNN"
)
method_colours <- c(
  "Raw" = "#6B7280", "ComBat" = "#0072B2", "QC-RLSC" = "#009E73",
  "QC-RFSC" = "#56B4E9", "TIGER" = "#E69F00", "SERRF" = "#D55E00",
  "WiNN (QC-tuned)" = "#CC79A7", "WiNN (inferred segments)" = "#7A5195",
  "WiNN" = "#B58B00"
)
gate_colours <- c(
  "Selective / selective" = "#0072B2", "Forced drift" = "#D55E00",
  "Forced batch" = "#009E73", "Forced drift + batch" = "#7A5195"
)
confounding_colours <- c(
  "Batch" = "#377EB8", "Run order" = "#1B9E77", "Both" = "#E64B35"
)
dataset_labels <- c(
  simulation = "Simulation", mtbls79 = "MTBLS79",
  clinical_fiams = "Clinical FIA-MS", batchcorr_set1 = "BatchCorr Set 1",
  sacurine = "Sacurine", waveica = "WaveICA"
)

theme_panel <- function(base_size = 9.5) {
  # Preserve a true 9 pt floor after the relative reductions below.
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
      plot.title.position = "plot",
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

panel_manifest <- list()
save_panel <- function(plot, directory, stem, width_mm, height_mm, source_data,
                       figure_id, panel_id, title) {
  pdf_path <- file.path(directory, paste0(stem, ".pdf"))
  png_path <- file.path(directory, paste0(stem, ".png"))
  csv_path <- file.path(panel_source_dir, paste0(stem, ".csv"))
  ggsave(
    pdf_path, plot, width = width_mm, height = height_mm, units = "mm",
    device = grDevices::pdf, useDingbats = FALSE, bg = "white"
  )
  source_data <- as.data.frame(source_data)
  if (all(c("method", "method_label") %in% names(source_data))) {
    source_data$method <- as.character(source_data$method_label)
  } else if ("method" %in% names(source_data) && "method_id" %in% names(source_data)) {
    display <- unname(method_labels[as.character(source_data$method_id)])
    source_data$method[!is.na(display)] <- display[!is.na(display)]
  }
  write_csv(source_data, csv_path, na = "")
  panel_manifest[[length(panel_manifest) + 1L]] <<- data.frame(
    figure = figure_id, panel = panel_id, title = title,
    pdf = file.path(basename(directory), basename(pdf_path)),
    png = NA_character_,
    source_csv = file.path("source_data", basename(csv_path)),
    width_mm = width_mm, height_mm = height_mm,
    stringsAsFactors = FALSE
  )
}

add_xy_se <- function(plot, data, x, x_se, y, y_se) {
  plot +
    geom_segment(
      data = data,
      aes(x = {{ x }} - {{ x_se }}, xend = {{ x }} + {{ x_se }},
          y = {{ y }}, yend = {{ y }}, colour = method_label),
      linewidth = line_width, alpha = 0.8
    ) +
    geom_segment(
      data = data,
      aes(x = {{ x }}, xend = {{ x }},
          y = {{ y }} - {{ y_se }}, yend = {{ y }} + {{ y_se }},
          colour = method_label),
      linewidth = line_width, alpha = 0.8
    )
}

# -----------------------------------------------------------------------------
# Main Figure 2: simulation truth recovery.
# -----------------------------------------------------------------------------

simulation_file <- file.path(revision_inputs, "simulation_seed_method_metrics.csv")
simulation <- read_csv(simulation_file, show_col_types = FALSE) |>
  filter(method_id %in% all_method_ids) |>
  mutate(
    method_label = unname(method_labels[method_id]),
    method_label = factor(method_label, levels = rev(unname(method_labels[all_method_ids])))
  )

simulation_specs <- list(
  list("A", "heldout_qc_cv_mean", "QC CV", "Held-out QC CV (%)", "lower"),
  list("B", "truth_log1p_rmse", "Clean-truth RMSE", "Clean-truth RMSE", "lower"),
  list("C", "truth_metabolite_profile_icc_mean", "Feature ICC(A,1)", "ICC(A,1) to clean truth", "higher"),
  list("D", "truth_sample_profile_pearson_mean", "Sample correlation", "Pearson r to clean truth", "higher")
)

for (spec in simulation_specs[2:4]) {
  panel_id <- spec[[1]]
  metric_id <- spec[[2]]
  panel_title <- spec[[3]]
  axis_title <- spec[[4]]
  panel_data <- simulation |>
    filter(metric == metric_id, method_id %in% simulation_method_ids) |>
    group_by(method_id, method_label) |>
    mutate(
      method_mean = mean(value),
      method_se = sd(value) / sqrt(n())
    ) |>
    ungroup()
  summary_data <- panel_data |>
    distinct(method_id, method_label, method_mean, method_se)

  plot <- ggplot(panel_data, aes(value, method_label, colour = method_label)) +
    geom_segment(
      data = summary_data,
      aes(x = method_mean - method_se, xend = method_mean + method_se,
          y = method_label, yend = method_label),
      linewidth = line_width, colour = "#172033"
    ) +
    geom_point(
      position = position_jitter(width = 0, height = 0.08, seed = 2026),
      shape = 21, fill = "white", size = 1.35, stroke = point_stroke, alpha = 0.72
    ) +
    geom_point(
      data = summary_data, aes(x = method_mean, y = method_label),
      inherit.aes = FALSE, shape = 23, fill = "white", colour = "#172033",
      size = 2.1, stroke = point_stroke
    ) +
    scale_colour_manual(values = method_colours, drop = FALSE) +
    scale_x_continuous(n.breaks = 4, expand = expansion(mult = c(0.06, 0.08))) +
    labs(title = paste(panel_id, panel_title), x = axis_title, y = NULL) +
    theme_panel(9.5) +
    theme(
      legend.position = "none", panel.grid.major.y = element_blank(),
      plot.title = element_text(size = 9.5)
    )

  save_panel(
    plot, main_dir, paste0("F2", panel_id, "_", metric_id), 85, 55,
    panel_data, "F2", panel_id, panel_title
  )
}


# Supplementary Figure S3 uses the same revised per-seed matrices as Figure 2.
s3_specs <- list(
  list("A", "heldout_qc_cv_mean", "Held-out QC CV", "Held-out QC CV (%)"),
  list("B", "residual_run_order_gam_mean_deviance", "Run-order structure", "GAM deviance explained"),
  list("C", "batch_weighted_pc_r2_categorical", "Residual batch structure", "Batch weighted-PC R2"),
  list("D", "truth_log1p_rmse", "Clean-truth RMSE", "Clean-truth RMSE"),
  list("E", "truth_metabolite_profile_icc_mean", "Feature ICC(A,1)", "ICC(A,1) to clean truth"),
  list("F", "truth_sample_profile_pearson_mean", "Sample correlation", "Pearson r to clean truth")
)
for (spec in s3_specs) {
  panel_data <- simulation |> filter(metric == spec[[2]]) |>
    group_by(method_id, method_label) |> mutate(method_mean = mean(value), method_se = sd(value)/sqrt(n())) |> ungroup()
  summary_data <- panel_data |> distinct(method_id, method_label, method_mean, method_se)
  plot <- ggplot(panel_data, aes(value, method_label, colour=method_label)) +
    geom_segment(data=summary_data, aes(x=method_mean-method_se, xend=method_mean+method_se,
      y=method_label,yend=method_label), colour="#172033", linewidth=line_width) +
    geom_point(position=position_jitter(width=0,height=.08,seed=2026), shape=21,
      fill="white",size=1.35,stroke=point_stroke,alpha=.72) +
    geom_point(data=summary_data,aes(x=method_mean),shape=23,fill="white",colour="#172033",size=2.1,stroke=point_stroke) +
    scale_colour_manual(values=method_colours,drop=FALSE) +
    scale_x_continuous(n.breaks=4,expand=expansion(mult=c(.06,.08))) +
    labs(title=paste(spec[[1]],spec[[3]]),x=spec[[4]],y=NULL) +
    theme_panel() + theme(legend.position="none",panel.grid.major.y=element_blank())
  save_panel(plot,supp_dir,paste0("S3",spec[[1]],"_simulation_",spec[[2]]),85,65,panel_data,"S3",spec[[1]],spec[[3]])
}

gate_truth <- read_csv(file.path(revision_inputs, "selective_forced_truth_rmse.csv"), show_col_types = FALSE) |>
  mutate(policy = factor(method_id, levels = c("G_AA", "G_SA", "G_AS", "G_SS"),
                         labels = c("Force both", "Force batch", "Force drift", "Selective")))
gate_summary <- gate_truth |> distinct(policy, method_mean, method_se)
p <- ggplot(gate_truth, aes(value, policy)) +
  geom_segment(data = gate_summary, aes(x = method_mean - method_se, xend = method_mean + method_se,
                                       y = policy, yend = policy), linewidth = line_width) +
  geom_point(position = position_jitter(width = 0, height = .075, seed = 2026),
             shape = 21, fill = "white", size = 1.35, stroke = point_stroke) +
  geom_point(data = gate_summary, aes(x = method_mean), shape = 23, fill = "white", size = 2.1,
             stroke = point_stroke) +
  scale_x_continuous(limits = c(.070, .080), breaks = c(.070, .075, .080),
                     expand = expansion(mult = c(.04, .12)), labels = scales::label_number(accuracy = .001)) +
  labs(title = "A Selective versus forced gates", x = "Clean-truth RMSE", y = NULL) +
  theme_panel() + theme(panel.grid.major.y = element_blank())
save_panel(p, main_dir, "F2A_selective_forced_truth_rmse", 85, 55, gate_truth,
           "F2", "A", "Selective versus forced gates")

# -----------------------------------------------------------------------------
# Main Figure 3: empirical technical and repeated-identity trade-offs.
# -----------------------------------------------------------------------------

identity_file <- file.path(revision_inputs, "eccb2026_identity_tradeoffs_v2.csv")
identity <- read_csv(identity_file, show_col_types = FALSE) |>
  filter(method_id %in% main_method_ids) |>
  mutate(
    method_label = unname(method_labels[method_id]),
    dataset_label = factor(dataset_label, levels = c("BatchCorr Set 1", "Clinical FIA-MS"))
  )

make_identity_scatter <- function(panel_data, panel_id, title, kind, stem) {
  if (kind == "qcv_gam") {
    base <- ggplot(panel_data)
    plot <- add_xy_se(
      base, panel_data, qc_cv, qc_cv_se,
      residual_run_order_gam_deviance_percent,
      residual_run_order_gam_deviance_percent_se
    ) +
      geom_point(
        aes(qc_cv, residual_run_order_gam_deviance_percent, colour = method_label),
        size = 1.7, shape = 21, fill = "white", stroke = point_stroke
      ) +
      geom_text_repel(
        aes(qc_cv, residual_run_order_gam_deviance_percent,
            label = method_label, colour = method_label),
        size = 3.2, seed = 2026, max.overlaps = Inf, force = 1.2,
        min.segment.length = 0, box.padding = 0.16, point.padding = 0.12,
        segment.size = line_width, show.legend = FALSE
      ) +
      labs(x = "Held-out reference CV (%)", y = "Residual run-order deviance (%)")
  } else {
    base <- ggplot(panel_data)
    plot <- add_xy_se(
      base, panel_data, batch_pc_variance_percent, batch_pc_variance_percent_se,
      point_estimate_genuine_replicate_icc_median,
      bootstrap_se_genuine_replicate_icc_median
    ) +
      geom_point(
        aes(batch_pc_variance_percent, point_estimate_genuine_replicate_icc_median,
            colour = method_label),
        size = 1.7, shape = 21, fill = "white", stroke = point_stroke
      ) +
      (if (panel_id == "C") {
        # Five corrected methods cluster near zero batch R2. Direct labels are
        # deliberately separated, with leaders linking every label to its point.
        label_data <- panel_data |>
          mutate(label_x = ifelse(method_id == "raw", 58, 17),
                 label_y = unname(c(raw = .886, combat = .941,
                   qc_rlsc = .952, tiger = .930, serrf = .919,
                   winn_fixed_no_qc = .963)[method_id]))
        list(
          geom_segment(data = label_data,
            aes(x = batch_pc_variance_percent,
                y = point_estimate_genuine_replicate_icc_median,
                xend = label_x - 1, yend = label_y, colour = method_label),
            linewidth = line_width),
          geom_text(data = label_data,
            aes(x = label_x, y = label_y, label = method_label,
                colour = method_label), hjust = 0, size = 3.2,
            show.legend = FALSE)
        )
      } else {
        geom_text_repel(
          aes(batch_pc_variance_percent, point_estimate_genuine_replicate_icc_median,
              label = method_label, colour = method_label),
          size = 3.2, seed = 2026, max.overlaps = Inf, force = 1.2,
          min.segment.length = 0, box.padding = 0.16, point.padding = 0.12,
          segment.size = line_width, show.legend = FALSE
        )
      }) +
      labs(
        x = ifelse(
          unique(as.character(panel_data$dataset_label)) == "Clinical FIA-MS",
          "Plate weighted-PC R² (%)",
          "Batch weighted-PC R² (%)"
        ),
        y = "Replicate-profile ICC(A,1)"
      )
  }

  plot <- plot +
    scale_colour_manual(values = method_colours, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.13, 0.27))) +
    scale_y_continuous(expand = expansion(mult = c(0.16, 0.24))) +
    coord_cartesian(clip = "off") +
    labs(title = paste(panel_id, title)) +
    theme_panel(9.4) +
    theme(legend.position = "none")

  save_panel(
    plot, main_dir, stem, 85, 55, panel_data, "F3", panel_id, title
  )
}

make_identity_scatter(
  filter(identity, dataset_key == "batchcorr_set1"), "A",
  "BatchCorr: reference precision", "qcv_gam",
  "F3A_batchcorr_qc_cv_vs_run_order_deviance"
)
make_identity_scatter(
  filter(identity, dataset_key == "clinical_fiams"), "B",
  "Clinical: reference precision", "qcv_gam",
  "F3B_clinical_qc_cv_vs_run_order_deviance"
)
make_identity_scatter(
  filter(identity, dataset_key == "batchcorr_set1"), "C",
  "BatchCorr: batch vs ICC", "batch_icc",
  "F3C_batchcorr_batch_variance_vs_replicate_icc"
)
make_identity_scatter(
  filter(identity, dataset_key == "clinical_fiams"), "D",
  "Clinical: plate vs ICC", "batch_icc",
  "F3D_clinical_plate_variance_vs_replicate_icc"
)

# -----------------------------------------------------------------------------
# Main Figure 2E-F: identifiability boundary.
# -----------------------------------------------------------------------------

confounding_file <- file.path(revision_inputs, "confounding_fixed_winn.csv")
confounding <- read_csv(confounding_file, show_col_types = FALSE) |>
  filter(metric %in% c("median_attenuation_ratio_clean_responsive", "empirical_fdr")) |>
  mutate(
    axis = recode(confounding_type, batch = "Batch", run_order = "Run order", joint = "Both"),
    endpoint = recode(
      metric,
      median_attenuation_ratio_clean_responsive = "Case-control effect recovery",
      empirical_fdr = "Empirical FDR"
    ),
    endpoint = factor(endpoint, levels = c("Case-control effect recovery", "Empirical FDR"))
  )

confounding_main_specs <- list(
  list(
    "E", "median_attenuation_ratio_clean_responsive",
    "WiNN: effect recovery", "Clean-reference recovery"
  ),
  list(
    "F", "empirical_fdr",
    "WiNN: empirical FDR", "Empirical FDR"
  )
)
for (spec in confounding_main_specs) {
  panel_data <- confounding |> filter(metric == spec[[2]])
  plot <- ggplot(panel_data, aes(alignment_percent, mean, colour = axis, group = axis)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 1.5, linewidth = line_width) +
    geom_line(linewidth = line_width) +
    geom_point(size = 1.6, shape = 16, stroke = 0) +
    scale_colour_manual(values = confounding_colours, breaks = c("Batch", "Run order", "Both")) +
    scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_y_continuous(limits = if (spec[[2]] == "empirical_fdr") c(0, 0.30) else c(0, 1.03),
                       breaks = if (spec[[2]] == "empirical_fdr") c(0, 0.05, 0.15, 0.25) else c(0, 0.5, 1)) +
    labs(
      title = paste(spec[[1]], spec[[3]]),
      x = "Alignment strength (%)", y = spec[[4]], colour = NULL
    ) +
    theme_panel(9.4) +
    theme(legend.position = "bottom", legend.direction = "horizontal")
  if (spec[[2]] == "empirical_fdr") plot <- plot + geom_hline(yintercept = .05, linetype = "dashed", linewidth = line_width)
  save_panel(
    plot, main_dir, paste0("F2", spec[[1]], "_", spec[[2]]), 85, 55,
    panel_data, "F2", spec[[1]], spec[[3]]
  )
}

public_effect_file <- file.path(revision_inputs, "eccb2026_public_effect_tradeoffs_v2.csv")
public_effect <- read_csv(public_effect_file, show_col_types = FALSE) |>
  filter(method_id %in% main_method_ids) |>
  mutate(
    method_label = unname(method_labels[method_id]),
    dataset_label = factor(
      dataset_label,
      levels = c("WaveICA: group 0 vs 1", "MTBLS79: cow vs sheep")
    )
  )

public_specs <- list(
  list("A", "waveica", "WaveICA group 0 vs 1", "F4A_waveica_batch_vs_group"),
  list("B", "mtbls79", "MTBLS79 cow vs sheep", "F4B_mtbls79_batch_vs_group")
)
for (spec in public_specs) {
  panel_data <- public_effect |> filter(dataset_key == spec[[2]])
  base <- ggplot(panel_data)
  plot <- add_xy_se(
    base, panel_data,
    batch_pc_variance_percent, batch_pc_variance_percent_se,
    biological_pc_variance_percent, biological_pc_variance_percent_se
  ) +
    geom_point(
      aes(batch_pc_variance_percent, biological_pc_variance_percent, colour = method_label),
      size = 1.7, shape = 21, fill = "white", stroke = point_stroke
    ) +
    geom_text_repel(
      aes(batch_pc_variance_percent, biological_pc_variance_percent,
          label = method_label, colour = method_label),
      size = 3.2, seed = 2026, max.overlaps = Inf, force = 1.25,
      min.segment.length = 0, box.padding = 0.16, point.padding = 0.12,
      segment.size = line_width, show.legend = FALSE
    ) +
    scale_colour_manual(values = method_colours, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.13, 0.28))) +
    scale_y_continuous(expand = expansion(mult = c(0.16, 0.24))) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste(spec[[1]], spec[[3]]),
      x = "Batch weighted-PC R² (%)",
      y = "Group weighted-PC R² (%)"
    ) +
    theme_panel(9.4) +
    theme(legend.position = "none")
  save_panel(
    plot, main_dir, spec[[4]], 85, 55, panel_data,
    "F4", spec[[1]], spec[[3]]
  )
}

# -----------------------------------------------------------------------------
partial_specs <- list(
  list("A", "heldout_qc_cv_mean", "Held-out QC CV", "Held-out QC CV (%)"),
  list("B", "residual_run_order_gam_mean_deviance", "Run-order structure", "GAM deviance explained"),
  list("C", "batch_weighted_pc_r2_categorical", "Residual batch structure", "Weighted-PC R2"),
  list("D", "clean_sample_profile_pearson_mean", "Clean sample-profile recovery", "Pearson r to clean data"),
  list("E", "median_attenuation_ratio_clean_responsive", "Clean-reference effect recovery", "Corrected / clean effect ratio"),
  list("F", "empirical_fdr", "Empirical false-discovery rate", "Empirical FDR")
)
for (spec in partial_specs) {
  panel_data <- read_csv(file.path(revision_inputs, paste0("S5", spec[[1]], "_partial_", spec[[2]], ".csv")), show_col_types = FALSE)
  panel_data$method_label <- ifelse(panel_data$method == "Raw", "Raw", "WiNN")
  plot <- ggplot(panel_data, aes(alignment_percent, mean, colour = method_label, group = method_label)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 1.5, linewidth = line_width) +
    geom_line(linewidth = line_width) +
    geom_point(size = 1.5, shape = 16, stroke = 0) +
    facet_wrap(vars(axis), nrow = 1) +
    scale_colour_manual(values = method_colours[c("Raw", "WiNN")], drop = FALSE) +
    scale_x_continuous(breaks = c(0, 50, 100)) +
    labs(
      title = paste(spec[[1]], spec[[3]]), x = "Alignment strength (%)",
      y = spec[[4]]
    ) +
    theme_panel(9.2) +
    theme(legend.position = "bottom", panel.spacing.x = unit(4, "mm"),
          plot.margin = margin(3, 9, 3, 4))
  if (spec[[2]] == "empirical_fdr") plot <- plot + geom_hline(yintercept = .05, linetype = "dashed", linewidth = line_width)
  save_panel(
    plot, supp_dir, paste0("S5", spec[[1]], "_partial_", spec[[2]]),
    85, 65, panel_data, "S5", spec[[1]], spec[[3]]
  )
}

# -----------------------------------------------------------------------------

write_csv(bind_rows(panel_manifest), file.path(package_root, "generated", "revision_main_panel_manifest.csv"))
