#!/usr/bin/env Rscript
# Rebuild S4, S8 and S10 from the original aggregate panel inputs.
# QC-RFSC is omitted because these analyses have not been repeated for its revised adapter.
# Every other estimate, uncertainty interval and hidden-reference split is unchanged.
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)))
root <- dirname(dirname(script))
suppressPackageStartupMessages({library(ggplot2); library(ggrepel)})
input <- file.path(root, "generated", "revision_panel_inputs")
out <- file.path(root, "figures", "supplementary_panels")
source_out <- file.path(root, "figures", "source_data")
ids <- c("raw", "combat", "qc_rlsc", "tiger", "serrf", "winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc")
labels <- setNames(c("Raw", "ComBat", "QC-RLSC", "TIGER", "SERRF", "WiNN (QC-tuned)", "WiNN (inferred segments)", "WiNN"), ids)
colors <- setNames(c("#6B7280", "#0072B2", "#009E73", "#E69F00", "#D55E00", "#CC79A7", "#7A5195", "#B58B00"), labels)
line_width <- 0.45 * 25.4 / 72.27
point_stroke <- 0.35 * 25.4 / 72.27
panel_theme <- theme_classic(base_size = 9.5, base_family = "sans") + theme(
 axis.line = element_line(colour = "#596273", linewidth = line_width),
 axis.ticks = element_line(colour = "#596273", linewidth = line_width),
 panel.grid.major = element_line(colour = "#E2E7ED", linewidth = line_width),
 panel.grid.minor = element_blank(),
 panel.background = element_rect(fill = "white", colour = NA),
 plot.background = element_rect(fill = "white", colour = NA),
 plot.title = element_text(face = "bold", colour = "#172033", size = 9.5),
 plot.title.position = "plot",
 axis.title = element_text(colour = "#263247", size = 9.3),
 axis.text = element_text(colour = "#263247", size = 9.1),
 legend.position = "none", panel.grid.major.y = element_blank(),
 plot.margin = margin(3, 4, 3, 4))
read_input <- function(stem, rows) {
 d <- read.csv(file.path(input, paste0(stem, ".csv")), check.names = FALSE)
 d <- d[d$method_id %in% ids, , drop = FALSE]
 stopifnot(nrow(d) == rows, setequal(unique(d$method_id), ids))
 d$method_label <- factor(unname(labels[d$method_id]), levels = rev(unname(labels)))
 if ("method" %in% names(d)) d$method <- unname(labels[d$method_id])
 d
}
save_panel <- function(p, stem, data, width = 85, height = 55) {
 ggsave(file.path(out, paste0(stem, ".pdf")), p, width = width, height = height,
        units = "mm", device = "pdf", useDingbats = FALSE, bg = "white")
 ggsave(file.path(out, paste0(stem, ".png")), p, width = width, height = height,
        units = "mm", dpi = 180, bg = "white")
 write.csv(data, file.path(source_out, paste0(stem, ".csv")), row.names = FALSE)
}
specs4 <- list(c("A","simulation","Simulation"), c("B","mtbls79","MTBLS79"),
 c("C","clinical_fiams","Clinical FIA-MS"), c("D","batchcorr_set1","BatchCorr Set 1"),
 c("E","sacurine","Sacurine"), c("F","waveica","WaveICA"))
for (spec in specs4) {
 stem <- paste0("S4", spec[1], "_reference_", spec[2]); d <- read_input(stem, 80)
 means <- unique(d[c("method_id", "method_label", "method_mean")])
 stopifnot(all(vapply(ids, function(id) abs(mean(d$value[d$method_id == id]) - means$method_mean[means$method_id == id]) < 1e-10, logical(1))))
 p <- ggplot(d, aes(value, method_label, colour = method_label)) +
  geom_point(position = position_jitter(width = 0, height = .08, seed = 2026),
             shape = 21, fill = "white", size = 1.3, stroke = point_stroke, alpha = .72) +
  geom_point(data = means, aes(method_mean, method_label), inherit.aes = FALSE,
             shape = 23, fill = "white", colour = "#172033", size = 2, stroke = point_stroke) +
  scale_colour_manual(values = colors, drop = FALSE) +
  labs(title = paste(spec[1], spec[3]), x = "Held-out QC CV (%)", y = NULL) + panel_theme
 save_panel(p, stem, d)
}
specs8 <- list(
 c("A", "batchcorr_set1", "batch_weighted_pc_r2_categorical", "BatchCorr batch structure", "Weighted-PC R2"),
 c("B", "batchcorr_set1", "sample_replicate_pearson_median", "BatchCorr replicate profiles", "Replicate-profile Pearson r"),
 c("C", "batchcorr_set1", "genuine_replicate_icc_median", "BatchCorr replicate agreement", "Replicate ICC(A,1)"),
 c("D", "clinical_fiams", "batch_weighted_pc_r2_categorical", "Clinical plate structure", "Weighted-PC R2"),
 c("E", "clinical_fiams", "sample_replicate_pearson_median", "Clinical replicate profiles", "Replicate-profile Pearson r"),
 c("F", "clinical_fiams", "genuine_replicate_icc_median", "Clinical replicate agreement", "Replicate ICC(A,1)"))
for (spec in specs8) {
 stem <- paste0("S8", spec[1], "_", spec[2], "_", spec[3]); d <- read_input(stem, 8)
 stopifnot(all(d$valid_bootstrap == 2000), all(d$ci_lower <= d$ci_upper))
 p <- ggplot(d, aes(point_estimate, method_label, colour = method_label)) +
  geom_segment(aes(x = ci_lower, xend = ci_upper, y = method_label, yend = method_label), linewidth = line_width) +
  geom_point(size = 1.8, shape = 16, stroke = 0) + scale_colour_manual(values = colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(.08, .16))) +
  labs(title = paste(spec[1], spec[4]), x = spec[5], y = NULL) + panel_theme +
  theme(plot.margin = margin(3, 12, 3, 4))
 save_panel(p, stem, d)
}
specs10 <- list(c("A", "age_weighted_pc_r2", "Sacurine age"),
 c("B", "bmi_weighted_pc_r2", "Sacurine BMI"), c("C", "gender_weighted_pc_r2", "Sacurine gender"))
for (spec in specs10) {
 stem <- paste0("S10", spec[1], "_sacurine_", spec[2]); d <- read_input(stem, 8)
 stopifnot(all(is.finite(as.matrix(d[c("x", "x_se", "y", "y_se")]))))
 # Put labels in one column, with leaders to the unchanged estimates.
 # Acquisition-variance removal clusters seven methods near x=0.
 label_data <- d
 selected <- which(label_data$method_id != "raw")
 selected <- selected[order(label_data$y[selected])]
 yr <- range(c(d$y-d$y_se, d$y+d$y_se))
 label_data$label_x <- label_data$x + .3
 label_data$label_y <- label_data$y
 label_data$label_x[selected] <- 2.2
 label_data$label_y[selected] <- seq(yr[1]+.15*diff(yr), yr[2]-.15*diff(yr), length.out=length(selected))
 p <- ggplot(d) +
  geom_segment(aes(x = x-x_se, xend = x+x_se, y = y, yend = y, colour = method_label), linewidth = line_width) +
  geom_segment(aes(x = x, xend = x, y = y-y_se, yend = y+y_se, colour = method_label), linewidth = line_width) +
  geom_point(aes(x, y, colour = method_label), shape = 21, fill = "white", size = 1.8, stroke = point_stroke) +
  geom_segment(data=label_data, aes(x=x, y=y, xend=label_x-.1, yend=label_y,
    colour=method_label), linewidth=line_width) +
  geom_text(data=label_data, aes(label_x, label_y, label=method_label, colour=method_label),
    size=3.2, hjust=0, show.legend=FALSE) +
  scale_colour_manual(values = colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(.12, .20))) +
  scale_y_continuous(expand = expansion(mult = c(.08, .10))) +
  labs(title = paste(spec[1], spec[3]), x = "Batch weighted-PC R² (%)", y = "Covariate weighted-PC R² (%)") +
  panel_theme + theme(panel.grid.major.y = element_line(colour = "#E2E7ED", linewidth = line_width),
                      plot.margin = margin(3, 10, 3, 4))
 save_panel(p, stem, d, width = 174, height = 60)
}
cat("Rebuilt 15 supplementary panels with unchanged non-RFSC inputs.\n")
