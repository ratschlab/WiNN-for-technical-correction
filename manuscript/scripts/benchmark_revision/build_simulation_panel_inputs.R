#!/usr/bin/env Rscript
# Builds two Figure 2 panel inputs from the release aggregate tables:
#   selective_forced_truth_rmse.csv  (F2A; per-seed gate-factorial truth RMSE)
#   confounding_fixed_winn.csv       (F2E/F; WiNN fixed partial-confounding summary)
# Usage: Rscript build_simulation_panel_inputs.R --tables=<release results/final/tables> --out=<dir>
suppressPackageStartupMessages(library(readr))
args <- commandArgs(TRUE)
arg <- function(name) {
  hit <- sub(paste0("^--", name, "="), "", grep(paste0("^--", name, "="), args, value = TRUE))
  if (length(hit) != 1L) stop("Supply --", name, "=", call. = FALSE)
  hit
}
tables <- arg("tables")
out <- arg("out")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# F2A ----------------------------------------------------------------------
abl <- read_csv(file.path(tables, "simulation_seed_ablation_metrics.csv"), show_col_types = FALSE,
                col_types = cols(.default = col_character()))
gate_labels <- c(G_SS = "Selective / selective", G_AS = "Forced drift",
                 G_SA = "Forced batch", G_AA = "Forced drift + batch")
gate <- abl[abl$metric == "truth_log1p_rmse" & abl$method_id %in% names(gate_labels), , drop = FALSE]
gate$value_num <- as.numeric(gate$value)
gate$gate_label <- unname(gate_labels[gate$method_id])
stats <- do.call(rbind, lapply(split(gate$value_num, gate$method_id), function(v) {
  data.frame(method_mean = mean(v), method_se = stats::sd(v) / sqrt(length(v)), n = length(v))
}))
stats$method_id <- rownames(stats)
if (any(stats$n != 10L)) stop("F2A: expected 10 seeds per gate policy.", call. = FALSE)
gate <- merge(gate, stats[, c("method_id", "method_mean", "method_se")], by = "method_id", sort = FALSE)
gate <- gate[order(gate$seed_id, match(gate$method_id, names(gate_labels))), , drop = FALSE]
gate_out <- data.frame(
  metric = gate$metric, value = gate$value_num, direction = gate$direction,
  denominator = as.numeric(gate$denominator), units = gate$units,
  notes = ifelse(is.na(gate$notes), "", gate$notes), dataset_key = gate$dataset_key,
  method_id = gate$method_id, method = gate$method, seed_id = gate$seed_id,
  gate_label = gate$gate_label, method_mean = gate$method_mean, method_se = gate$method_se
)
write_csv(gate_out, file.path(out, "selective_forced_truth_rmse.csv"), na = "")

# F2E/F --------------------------------------------------------------------
conf <- read_csv(file.path(tables, "partial_confounding_method_metrics.csv"), show_col_types = FALSE)
conf <- conf[conf$method == "WiNN fixed default (no QC)" &
               conf$metric %in% c("median_attenuation_ratio_clean_responsive",
                                  "true_positive_rate", "empirical_fdr"), , drop = FALSE]
types <- c("batch", "run_order", "joint")
type_labels <- c(batch = "Batch-aligned biology", run_order = "Run-order-aligned biology",
                 joint = "Jointly aligned biology")
endpoint_labels <- c(batch = "Batch", run_order = "Run order", joint = "Joint")
levels <- c(0, 0.25, 0.5, 0.75, 0.9, 1)
rows <- list()
for (type in types) for (metric in c("median_attenuation_ratio_clean_responsive",
                                     "true_positive_rate", "empirical_fdr")) for (level in levels) {
  source_type <- if (level == 0) "none" else type
  v <- conf$value[conf$confounding_type == source_type & abs(conf$level - level) < 1e-12 &
                    conf$metric == metric]
  v <- v[is.finite(v)]
  m <- mean(v); s <- stats::sd(v); se <- s / sqrt(length(v))
  labelled <- metric != "empirical_fdr"
  rows[[length(rows) + 1L]] <- data.frame(
    confounding_type = type, level = level, metric = metric, mean = m, sd = s, se = se,
    n_seeds = length(v),
    confounding_label = if (labelled) type_labels[[type]] else "",
    alignment_percent = 100 * level,
    endpoint_label = if (labelled) endpoint_labels[[type]] else "",
    lower = m - se, upper = m + se
  )
}
conf_out <- do.call(rbind, rows)
# Original layout: labelled endpoints first (type alphabetical, level, metric), then FDR rows.
conf_out <- conf_out[order(conf_out$metric == "empirical_fdr", conf_out$confounding_type,
                           conf_out$level, conf_out$metric), ]
write_csv(conf_out, file.path(out, "confounding_fixed_winn.csv"), na = "")
cat("wrote", nrow(gate_out), "F2A rows and", nrow(conf_out), "F2E/F rows to", out, "\n")
