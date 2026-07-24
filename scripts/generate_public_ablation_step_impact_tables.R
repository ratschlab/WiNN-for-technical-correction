#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript or R --file.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly = TRUE)
combined_arg <- grep("^--combined-dir=", args, value = TRUE)

combined_dir <- if (length(combined_arg) == 1L) {
  normalizePath(sub("^--combined-dir=", "", combined_arg), mustWork = TRUE)
} else {
  file.path(repo_root, "results", "winn_ablations", "combined")
}
input_path <- file.path(combined_dir, "cumulative_primary_metrics.csv")
if (!file.exists(input_path)) stop("Run scripts/summarize_public_winn_ablations.R first.")
metrics <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
metrics$value <- as.numeric(metrics$value)

datasets <- c("simulation", "mtbls79", "batchcorr_set1", "sacurine", "waveica")
dataset_labels <- c(
  simulation = "Simulation",
  mtbls79 = "MTBLS79",
  batchcorr_set1 = "BatchCorr Set 1",
  sacurine = "Sacurine",
  waveica = "WaveICA"
)
stages <- c("C0_RAW", "C1_OUTLIER", "C2_SELECTIVE_DRIFT", "C3_SELECTIVE_BATCH", "C4_FULL_FIXED")
stage_labels <- c("Raw", "Outlier", "+ selective drift", "+ selective batch", "+ PQN")

endpoint_specs <- list(
  heldout_cv = list(
    title = "Held-out QC/reference CV mean (%)",
    direction = "lower is better",
    units = "percent",
    stem = "heldout_cv",
    mapping = stats::setNames(rep("heldout_qc_cv_mean", length(datasets)), datasets)
  ),
  gam_deviance = list(
    title = "Residual run-order GAM deviance explained (mean)",
    direction = "lower is better",
    units = "proportion",
    stem = "gam_deviance",
    mapping = stats::setNames(rep("residual_gam_deviance_mean", length(datasets)), datasets)
  ),
  batch_r2 = list(
    title = "Batch/plate weighted-PC R-squared",
    direction = "lower is better",
    units = "weighted R-squared",
    stem = "batch_r2",
    mapping = stats::setNames(rep("batch_weighted_pc_r2", length(datasets)), datasets)
  ),
  metabolite_icc = list(
    title = "Metabolite ICC(A,1)",
    direction = "higher is better",
    units = "ICC(A,1)",
    stem = "metabolite_icc_a1",
    mapping = c(
      simulation = "truth_metabolite_profile_icc_mean",
      mtbls79 = "cross_batch_metabolite_icc_median",
      batchcorr_set1 = NA_character_,
      sacurine = NA_character_,
      waveica = NA_character_
    )
  ),
  sample_pearson = list(
    title = "Sample Pearson agreement",
    direction = "higher is better",
    units = "correlation",
    stem = "sample_pearson",
    mapping = c(
      simulation = "truth_sample_profile_pearson_mean",
      mtbls79 = "sample_replicate_pearson_median",
      batchcorr_set1 = "genuine_replicate_pearson_median",
      sacurine = NA_character_,
      waveica = NA_character_
    )
  )
)

unavailable_reason <- c(
  simulation = "",
  mtbls79 = "",
  batchcorr_set1 = "No per-metabolite ICC(A,1) is defined in the saved schema; feature repeatability is not relabeled as ICC.",
  sacurine = "No genuine repeated biological sample measurements.",
  waveica = "No genuine repeated biological sample measurements."
)

out_dir <- file.path(combined_dir, "step_impact_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fmt <- function(value) {
  if (!is.finite(value)) return("NA")
  if (abs(value) >= 10) sprintf("%.2f", value) else sprintf("%.4f", value)
}
fmt_delta <- function(value) {
  if (!is.finite(value)) return("NA")
  paste0(if (value >= 0) "+" else "", fmt(value))
}

tables <- list()
long_rows <- list()
mapping_rows <- list()
for (endpoint_name in names(endpoint_specs)) {
  spec <- endpoint_specs[[endpoint_name]]
  rows <- lapply(datasets, function(dataset) {
    metric_name <- unname(spec$mapping[[dataset]])
    values <- rep(NA_real_, length(stages))
    if (!is.na(metric_name)) {
      selected <- metrics[
        metrics$dataset == dataset & metrics$metric == metric_name & metrics$variant_id %in% stages,
        c("variant_id", "value"), drop = FALSE
      ]
      values <- selected$value[match(stages, selected$variant_id)]
    }
    deltas <- c(NA_real_, diff(values))
    data.frame(
      dataset = dataset,
      dataset_label = unname(dataset_labels[[dataset]]),
      endpoint_metric = if (is.na(metric_name)) NA_character_ else metric_name,
      endpoint_definition = spec$title,
      direction = spec$direction,
      raw_value = values[1],
      outlier_value = values[2], outlier_delta = deltas[2],
      selective_drift_value = values[3], selective_drift_delta = deltas[3],
      selective_batch_value = values[4], selective_batch_delta = deltas[4],
      pqn_value = values[5], pqn_delta = deltas[5],
      unavailable_reason = if (is.na(metric_name)) unname(unavailable_reason[[dataset]]) else "",
      stringsAsFactors = FALSE
    )
  })
  table_numeric <- do.call(rbind, rows)
  write.csv(
    table_numeric,
    file.path(out_dir, paste0(spec$stem, "_step_impact.csv")),
    row.names = FALSE, quote = TRUE, na = ""
  )
  tables[[endpoint_name]] <- table_numeric

  long_rows[[endpoint_name]] <- do.call(rbind, lapply(seq_len(nrow(table_numeric)), function(i) {
    row <- table_numeric[i, ]
    values <- as.numeric(row[c("raw_value", "outlier_value", "selective_drift_value", "selective_batch_value", "pqn_value")])
    data.frame(
      endpoint = endpoint_name,
      endpoint_definition = spec$title,
      direction = spec$direction,
      units = spec$units,
      dataset = row$dataset,
      dataset_label = row$dataset_label,
      endpoint_metric = row$endpoint_metric,
      variant_id = stages,
      stage = stage_labels,
      value = values,
      previous_value = c(NA_real_, values[-length(values)]),
      absolute_delta = c(NA_real_, diff(values)),
      improvement = if (spec$direction == "lower is better") c(NA_real_, -diff(values)) else c(NA_real_, diff(values)),
      unavailable_reason = row$unavailable_reason,
      stringsAsFactors = FALSE
    )
  }))

  mapping_rows[[endpoint_name]] <- data.frame(
    endpoint = endpoint_name,
    endpoint_definition = spec$title,
    direction = spec$direction,
    units = spec$units,
    dataset = datasets,
    dataset_label = unname(dataset_labels[datasets]),
    source_metric = unname(spec$mapping[datasets]),
    unavailable_reason = ifelse(is.na(unname(spec$mapping[datasets])), unname(unavailable_reason[datasets]), ""),
    stringsAsFactors = FALSE
  )
}

write.csv(do.call(rbind, long_rows), file.path(out_dir, "all_step_impacts_long.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(do.call(rbind, mapping_rows), file.path(out_dir, "endpoint_mapping.csv"), row.names = FALSE, quote = TRUE, na = "")

md <- c(
  "# WiNN cumulative-stage impact tables",
  "",
  "Each corrected-stage cell reports the absolute endpoint value followed by the change from the immediately preceding stage. Positive `improvement` in the long source table always means movement in the endpoint's preferred direction. `NA` is retained when a design does not support a genuine endpoint; no proxy is silently relabeled.",
  ""
)
for (endpoint_name in names(endpoint_specs)) {
  spec <- endpoint_specs[[endpoint_name]]
  tab <- tables[[endpoint_name]]
  md <- c(md, paste0("## ", spec$title), "", paste0("Direction: ", spec$direction, "."), "")
  header <- "| Dataset | Raw | Outlier | + selective drift | + selective batch | + PQN |"
  separator <- "|---|---:|---:|---:|---:|---:|"
  body <- vapply(seq_len(nrow(tab)), function(i) {
    row <- tab[i, ]
    values <- as.numeric(row[c("raw_value", "outlier_value", "selective_drift_value", "selective_batch_value", "pqn_value")])
    deltas <- c(NA_real_, as.numeric(row[c("outlier_delta", "selective_drift_delta", "selective_batch_delta", "pqn_delta")]))
    cells <- c(fmt(values[1]), vapply(2:5, function(j) {
      if (!is.finite(values[j])) "NA" else paste0(fmt(values[j]), " (Delta ", fmt_delta(deltas[j]), ")")
    }, character(1)))
    paste0("| ", row$dataset_label, " | ", paste(cells, collapse = " | "), " |")
  }, character(1))
  md <- c(md, header, separator, body, "")
  unavailable <- tab[nzchar(tab$unavailable_reason), c("dataset_label", "unavailable_reason"), drop = FALSE]
  if (nrow(unavailable)) {
    md <- c(md, "Unavailable endpoints:", "", paste0("- ", unavailable$dataset_label, ": ", unavailable$unavailable_reason), "")
  }
}
writeLines(md, file.path(out_dir, "WINN_ABLATION_STEP_IMPACT_TABLES.md"))

cat("Five cumulative-stage impact tables written to:", normalizePath(out_dir), "\n")
