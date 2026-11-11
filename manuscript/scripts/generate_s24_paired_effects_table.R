#!/usr/bin/env Rscript
# Table S24: same-feature effects relative to input abundance.
# Renders generated/supplementary_tables/table_s24_paired_feature_effects.tex from
# generated/review_validation/paired_feature_effects_descriptive.csv (written by
# experiments/run_targeted_biology_audit.R). Raw (the reference) and the previous
# QC-RFSC adapter are excluded, as in the published table. No analysis is run here.
#
# Usage (from the manuscript directory):
#   Rscript scripts/generate_s24_paired_effects_table.R [--input=<csv>] [--output=<tex>]
args <- commandArgs(TRUE)
opt <- function(name, default) {
  hit <- sub(paste0("^--", name, "="), "", grep(paste0("^--", name, "="), args, value = TRUE))
  if (length(hit)) hit[1] else default
}
file_arg <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE)
root <- dirname(dirname(normalizePath(file_arg)))
input <- opt("input", file.path(root, "generated", "review_validation", "paired_feature_effects_descriptive.csv"))
output <- opt("output", file.path(root, "generated", "supplementary_tables", "table_s24_paired_feature_effects.tex"))

effects <- read.csv(input, stringsAsFactors = FALSE)
endpoint_order <- data.frame(
  dataset_key = c("mtbls79", "waveica", "sacurine", "sacurine", "sacurine"),
  endpoint = c("cow_vs_sheep", "biological_group", "age", "bmi", "gender"),
  label = c("MTBLS79 cow--sheep", "WaveICA groups", "Sacurine age", "Sacurine bmi", "Sacurine gender"),
  stringsAsFactors = FALSE
)
method_labels <- c(
  combat = "ComBat", qc_rlsc = "QC-RLSC", tiger = "TIGER", serrf = "SERRF",
  winn_auto_qc = "WiNN (QC-tuned)", winn_auto_batch_qc = "WiNN (inferred segments)",
  winn_fixed_no_qc = "WiNN"
)
rows <- character()
for (i in seq_len(nrow(endpoint_order))) {
  part <- effects[effects$dataset_key == endpoint_order$dataset_key[i] &
                    effects$endpoint == endpoint_order$endpoint[i], , drop = FALSE]
  for (method_id in names(method_labels)) {
    row <- part[part$method_id == method_id, , drop = FALSE]
    if (nrow(row) != 1L) stop("Missing or duplicated S24 row: ", endpoint_order$label[i], " / ", method_id)
    rows <- c(rows, sprintf(
      "%s & %s & %.3f & %.3f & %.3f \\\\", endpoint_order$label[i], method_labels[[method_id]],
      row$slope_through_origin_vs_raw, row$sign_concordance_vs_raw, row$pearson_vs_raw
    ))
  }
}
header <- c(
  "\\scriptsize",
  "\\begin{longtable}{llrrr}",
  paste0(
    "\\caption{Same-feature effects relative to input abundance. Slope is the through-origin ",
    "least-squares slope of corrected on input feature effects. Sign concordance and Pearson ",
    "correlation use matched finite features. Effects use the log1p scale and the dataset-specific ",
    "osmolality adjustment.}\\label{tab:supp_paired_effects}\\\\"
  ),
  "\\toprule",
  "Dataset / endpoint & Method & Slope & Sign agreement & Pearson \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\toprule",
  "Dataset / endpoint & Method & Slope & Sign agreement & Pearson \\\\",
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\endlastfoot"
)
writeLines(c(header, rows, "\\end{longtable}"), output, useBytes = TRUE)
cat("wrote", output, "with", length(rows), "rows\n")
