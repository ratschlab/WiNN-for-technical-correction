#!/usr/bin/env Rscript
# Render actual fold-specific conditional bootstrap estimates; no correction rerun.
args_all <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", grep("^--file=", args_all, value = TRUE)))
root <- dirname(dirname(script))
stopifnot(requireNamespace("ggplot2", quietly = TRUE))
library(ggplot2)
d <- read.csv(file.path(root, "generated", "supplementary_table_data",
                        "leave_one_batch_out_bootstrap_summary.csv"), check.names = FALSE)
method_ids <- c("raw", "combat", "qc_rlsc", "tiger", "serrf",
                "winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc")
method_labels <- c("Raw", "ComBat", "QC-RLSC", "TIGER", "SERRF",
                   "WiNN (QC-tuned)", "WiNN (inferred segments)", "WiNN")
d <- d[d$method_id %in% method_ids, ]
d$method_label <- factor(d$method_id, levels = rev(method_ids), labels = rev(method_labels))
d$metric <- factor(d$metric, levels = c("pearson", "spearman", "sign_concordance"),
                   labels = c("Pearson effect\ncorrelation", "Spearman effect\ncorrelation", "Sign\nconcordance"))
d <- d[!is.na(d$metric), ]
panels <- list(
  list("mtbls79", "cow_vs_sheep", "S7A_lobo_mtbls79_cow_vs_sheep", "A  MTBLS79 cow versus sheep"),
  list("waveica", "biological_group", "S7B_lobo_waveica_biological_group", "B  WaveICA deposited groups"),
  list("sacurine", "age", "S7C_lobo_sacurine_age", "C  Sacurine age"),
  list("sacurine", "bmi", "S7D_lobo_sacurine_bmi", "D  Sacurine BMI"),
  list("sacurine", "gender", "S7E_lobo_sacurine_gender", "E  Sacurine gender")
)
shapes <- c(16, 17, 15, 18, 0, 1, 2, 5)
colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00", "#333333", "#882255")
out <- file.path(root, "figures", "supplementary_panels")
for (spec in panels) {
  x <- d[d$dataset_key == spec[[1]] & d$endpoint == spec[[2]], ]
  stopifnot(nrow(x) > 0)
  x$heldout_batch <- factor(x$heldout_batch, levels = unique(x$heldout_batch))
  x$y <- as.integer(x$method_label) +
    (as.integer(x$heldout_batch) - (nlevels(x$heldout_batch) + 1) / 2) *
    min(0.075, 0.5 / nlevels(x$heldout_batch))
  p <- ggplot(x, aes(x = bootstrap_mean, y = y, color = heldout_batch, shape = heldout_batch)) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.2) +
    geom_segment(aes(x = ci_lower, xend = ci_upper, yend = y), linewidth = 0.23) +
    geom_point(size = 1.6, stroke = 0.35) +
    facet_wrap(~metric, nrow = 1, scales = "free_x") +
    scale_x_continuous(breaks = scales::breaks_extended(n = 3)) +
    scale_color_manual(values = colors, name = "Held-out batch") +
    scale_shape_manual(values = shapes, name = "Held-out batch") +
    scale_y_continuous(breaks = seq_along(method_ids), labels = rev(method_labels),
                       limits = c(0.5, length(method_ids) + 0.5)) +
    labs(title = spec[[4]], x = "Fold bootstrap mean and 95% percentile interval", y = NULL,
         subtitle = "Effects estimated separately by batch after joint correction",
         caption = "Intervals describe individual folds after joint correction.") +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.spacing.x = unit(8, "mm"),
          panel.border = element_rect(linewidth = 0.35), strip.background = element_rect(fill = "grey96", linewidth = 0.3),
          strip.text = element_text(size = 10), legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 10),
          plot.caption = element_text(size = 9.6), axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 10), legend.text = element_text(size = 10))
  write.csv(x, file.path(root, "figures", "source_data", paste0(spec[[3]], ".csv")), row.names = FALSE)
  ggsave(file.path(out, paste0(spec[[3]], ".pdf")), p, width = 6.8, height = 7.0,
         device = "pdf", useDingbats = FALSE)
}
