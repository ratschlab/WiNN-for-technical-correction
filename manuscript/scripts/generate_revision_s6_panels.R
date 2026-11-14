#!/usr/bin/env Rscript
# Descriptive permutation-tail and absolute-effect plots from frozen outputs.
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)))
root <- dirname(dirname(script))
suppressPackageStartupMessages(library(ggplot2))
source(file.path(dirname(script), "method_palette.R"))
input <- file.path(root, "generated", "revision_panel_inputs")
counts <- read.csv(gzfile(file.path(input, "permutation_null.csv.gz")), stringsAsFactors = FALSE)
summary <- read.csv(file.path(input, "permutation_summary.csv"), stringsAsFactors = FALSE)
effects <- read.csv(file.path(input, "biological_association_summary.csv"), stringsAsFactors = FALSE)
ids <- c("raw", "combat", "qc_rlsc", "tiger", "serrf", "winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc")
labels <- c("Raw", "ComBat", "QC-RLSC", "TIGER", "SERRF", "WiNN (QC-tuned)", "WiNN (inferred segments)", "WiNN")
names(labels) <- ids
colors <- palette_colours_by_id[ids]
specs <- list(
  c("mtbls79", "cow_vs_sheep", "A", "B", "MTBLS79 cow vs sheep", "mtbls79_cow_vs_sheep"),
  c("waveica", "biological_group", "C", "D", "WaveICA biological groups", "waveica_biological_group"),
  c("sacurine", "age", "E", "F", "Sacurine age", "sacurine_age"),
  c("sacurine", "bmi", "G", "H", "Sacurine BMI", "sacurine_bmi"),
  c("sacurine", "gender", "I", "J", "Sacurine gender", "sacurine_gender")
)
out <- file.path(root, "figures", "supplementary_panels")
source_out <- file.path(root, "figures", "source_data")
dir.create(source_out, showWarnings = FALSE, recursive = TRUE)
base_theme <- theme_bw(base_size = 9.5) + theme(panel.grid.minor = element_blank(),
  panel.border = element_rect(linewidth = .3), axis.text = element_text(size = 9),
  axis.title = element_text(size = 9.3), plot.title = element_text(size = 10, face = "bold"),
  plot.subtitle = element_text(size = 9), plot.caption = element_text(size = 9),
  legend.text = element_text(size = 9), legend.title = element_blank(),
  plot.title.position = "plot", plot.caption.position = "plot",
  plot.margin = margin(4, 10, 4, 3))
log_breaks <- function(limits) {
  lo <- ceiling(log10(limits[1])); hi <- floor(log10(limits[2]))
  value <- if (lo <= hi) 10^seq(lo, hi) else numeric()
  if (length(value) < 2L) value <- signif(exp(seq(log(limits[1]), log(limits[2]), length.out = 4))[2:3], 2)
  if (length(value) > 2L) value <- value[c(1L, length(value))]
  value[value > 0]
}
short_scientific <- function(value) vapply(value, function(x) {
  if (!is.finite(x) || x == 0) return(as.character(x))
  exponent <- floor(log10(abs(x)))
  if (exponent >= -2 && exponent <= 1) return(format(signif(x, 2), scientific = FALSE, trim = TRUE))
  mantissa <- signif(x / 10^exponent, 2)
  paste0(format(mantissa, scientific = FALSE, trim = TRUE), "e", exponent)
}, character(1))
save_plot <- function(p, name) {
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = 174, height = 70, units = "mm", device = "pdf", useDingbats = FALSE)
  ggsave(file.path(out, paste0(name, ".png")), p, width = 174, height = 70, units = "mm", dpi = 180)
}
for (spec in specs) {
  d <- counts[counts$dataset_key == spec[1] & counts$endpoint == spec[2], ]
  obs <- summary[summary$dataset_key == spec[1] & summary$endpoint == spec[2], ]
  tails <- do.call(rbind, lapply(ids, function(id) {
    values <- d$significant_features[d$method_id == id]
    observed <- obs$observed_count[obs$method_id == id]
    stopifnot(length(values) == 5000L, length(observed) == 1L)
    p <- (1 + sum(values >= observed)) / (length(values) + 1)
    stopifnot(abs(p - obs$empirical_p_value[obs$method_id == id]) < 1e-12)
    data.frame(method_id = id, method = labels[id], observed = observed,
      null_median = median(values), null_q95 = unname(quantile(values, .95)),
      null_q99 = unname(quantile(values, .99)), null_q999 = unname(quantile(values, .999)),
      null_max = max(values), probability_any_null_discovery = mean(values > 0),
      exceedances = sum(values >= observed), permutations = length(values), empirical_p = p)
  }))
  tails$method_label <- factor(tails$method_id, levels = rev(ids), labels = labels[rev(ids)])
  tails$p_label <- sprintf("p=%.4f", tails$empirical_p)
  max_count <- max(c(tails$observed, tails$null_max))
  p_x <- expm1(log1p(max_count) + .22)
  upper <- expm1(log1p(max_count) + 2.25)
  breaks <- if (max_count <= 10) c(0,1,2,5,10) else if (max_count <= 100) c(0,1,5,20,100) else c(0,1,10,100,1000,10000)
  breaks <- breaks[breaks <= max_count]
  null_plot <- ggplot(tails, aes(y = method_label)) +
    geom_segment(aes(x = 0, xend = null_max, yend = method_label), color = "#B8BEC6", linewidth = .35) +
    geom_point(aes(x = null_q99, shape = "Null 99%"), color = "#4D4D4D", size = 1.8) +
    geom_point(aes(x = null_max, shape = "Null max"), color = "#4D4D4D", size = 1.8) +
    geom_point(aes(x = observed, shape = "Observed", color = method_id), size = 2.4) +
    scale_color_manual(values = colors, guide = "none") +
    geom_text(aes(x = p_x, label = p_label), hjust = 0, size = 3.17, color = "#333333") +
    scale_shape_manual(values = c("Null 99%" = 1, "Null max" = 17, "Observed" = 18),
      breaks = c("Null 99%", "Null max", "Observed")) +
    scale_x_continuous(trans = "log1p", breaks = breaks, limits = c(0, upper), expand = c(0, .02)) +
    labs(title = paste(spec[3], spec[5]), subtitle = "Null tails and observed counts",
      x = "Significant count (log1p axis)", y = NULL,
      caption = "5,000 constrained permutations.\np = (1 + null counts >= observed) / 5,001.") +
    base_theme + theme(legend.position = "top", legend.box.margin = margin(0,0,0,0),
      legend.margin = margin(0,0,0,0), panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 9), legend.key.width = unit(3, "mm"), legend.key.height = unit(3, "mm"))
  name <- paste0("S6", spec[3], "_permutation_", spec[6])
  write.csv(tails, file.path(source_out, paste0(name, ".csv")), row.names = FALSE)
  save_plot(null_plot, name)

  e <- effects[effects$dataset_key == spec[1] & effects$endpoint == spec[2] & effects$method_id %in% ids, ]
  if ("method" %in% names(e)) e$method <- unname(labels[e$method_id])
  e$short_label <- unname(labels[e$method_id])
  e$method_label <- e$short_label
  stopifnot(nrow(e) == length(ids), all(e$median_absolute_effect > 0), all(e$median_within_group_or_residual_variance > 0))
  long <- rbind(data.frame(method_id = e$method_id, statistic = "Absolute effect", value = e$median_absolute_effect),
    data.frame(method_id = e$method_id, statistic = "Variance", value = e$median_within_group_or_residual_variance))
  long$method_label <- factor(long$method_id, levels = rev(ids), labels = labels[rev(ids)])
  long$statistic <- factor(long$statistic, levels = c("Absolute effect", "Variance"))
  effect_plot <- ggplot(long, aes(x = value, y = method_label, color = method_id)) +
    geom_point(size = 2) + facet_wrap(~statistic, nrow = 1, scales = "free_x") +
    scale_color_manual(values = colors) +
    scale_x_log10(labels = short_scientific, breaks = log_breaks,
      expand = expansion(mult = c(.15, .15))) +
    labs(title = paste(spec[4], spec[5]), subtitle = "Absolute effects and variances",
      x = "Median value (log axis)", y = NULL,
      caption = if (spec[1] == "sacurine") "Effects: log1p after osmolality scaling." else
        "Effects: log1p intensity scale.") +
    base_theme + theme(legend.position = "none", panel.grid.major.y = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(size = 9), strip.background = element_rect(fill = "grey96", linewidth = .3))
  name <- paste0("S6", spec[4], "_effect_variance_", spec[6])
  write.csv(e, file.path(source_out, paste0(name, ".csv")), row.names = FALSE)
  save_plot(effect_plot, name)
}
cat("Regenerated all ten S6 panels from full frozen permutation counts and absolute effects.\n")
