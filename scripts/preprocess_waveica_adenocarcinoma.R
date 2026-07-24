#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript or R --file.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
source(file.path(repo_root, "scripts", "public_benchmark_audit_helpers.R"))
require_packages(c("jsonlite", "openssl"))

force <- "--force" %in% commandArgs(trailingOnly = TRUE)
source_dir <- file.path(repo_root, "data", "public", "waveica_adenocarcinoma", "source")
processed_dir <- file.path(repo_root, "data", "public", "waveica_adenocarcinoma", "processed")
audit_dir <- file.path(repo_root, "data", "public", "waveica_adenocarcinoma", "audit")
report_dir <- file.path(repo_root, "reports")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

outputs <- c(
  matrix = file.path(processed_dir, "waveica_intensity_matrix.rds"),
  metadata = file.path(processed_dir, "waveica_metadata.csv"),
  features = file.path(processed_dir, "waveica_feature_metadata.csv"),
  summary = file.path(processed_dir, "waveica_preprocessing_summary.json"),
  removed = file.path(processed_dir, "waveica_removed_features.csv"),
  report_md = file.path(report_dir, "waveica_adenocarcinoma_suitability.md"),
  report_csv = file.path(report_dir, "waveica_adenocarcinoma_suitability.csv")
)
if (!force && all(file.exists(outputs))) {
  message("All WaveICA audit/preprocessing outputs exist; reusing them.")
  quit(save = "no", status = 0)
}

source_path <- file.path(source_dir, "Amide_data.rda")
if (!file.exists(source_path)) stop("Missing Amide_data.rda. Run scripts/download_waveica_adenocarcinoma.py.")
loaded <- new.env(parent = emptyenv())
loaded_names <- load(source_path, envir = loaded)
if (!identical(loaded_names, "Amide_data")) stop("Expected only the Amide_data object in Amide_data.rda; found: ", paste(loaded_names, collapse = ", "))
source_data <- loaded$Amide_data
if (!is.data.frame(source_data) || !all(c("Injection_order", "group", "batch") %in% names(source_data))) stop("Amide_data is not the documented data frame.")
if (ncol(source_data) <= 3L) stop("Amide_data contains no feature columns.")

sample_ids <- rownames(source_data)
if (is.null(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) sample_ids <- sprintf("WaveICA_%04d", seq_len(nrow(source_data)))
feature_source_ids <- names(source_data)[-(1:3)]
intensity <- t(as.matrix(source_data[, -(1:3), drop = FALSE]))
storage.mode(intensity) <- "double"
rownames(intensity) <- make.unique(feature_source_ids, sep = "__duplicate_")
colnames(intensity) <- sample_ids
if (anyDuplicated(rownames(intensity)) || anyDuplicated(colnames(intensity))) stop("WaveICA has unresolved duplicate feature or sample IDs.")
if (any(!is.finite(intensity)) || any(intensity < 0)) stop("WaveICA contains missing, nonfinite, or negative intensities.")
if (anyNA(source_data$Injection_order) || anyDuplicated(source_data$Injection_order)) stop("WaveICA injection order is missing or not globally unique.")
if (!identical(as.integer(source_data$Injection_order), seq_len(nrow(source_data)))) stop("WaveICA rows are not in globally increasing deposited injection order.")

is_qc <- as.character(source_data$group) == "QC"
group_neutral <- ifelse(is_qc, NA_character_, paste0("group_", as.character(source_data$group)))
# The source object has exactly 497 group-1 and 71 group-0 study samples, matching
# the two published cohort totals. Counts alone do not prove the row-level label
# direction, so disease names are deliberately not assigned.
disease_label <- rep(NA_character_, nrow(source_data))
batch_chr <- as.character(source_data$batch)
within_order <- ave(source_data$Injection_order, batch_chr, FUN = function(x) rank(x, ties.method = "first"))
metadata <- data.frame(
  sample_id = sample_ids,
  original_row_name = rownames(source_data),
  original_row_index = seq_len(nrow(source_data)),
  sample_type = ifelse(is_qc, "QC", "Study"),
  class = ifelse(is_qc, "QC", "Sample"),
  is_qc = is_qc,
  biological_group = group_neutral,
  deposited_group = as.character(source_data$group),
  disease_label = disease_label,
  disease_label_mapping_evidence = ifelse(is_qc, NA_character_, "not proven: deposited 0/1 counts match published disease totals, but no source maps group values to diagnoses"),
  batch = batch_chr,
  deposited_injection_order = as.integer(source_data$Injection_order),
  within_batch_order = as.integer(within_order),
  global_run_order = as.integer(source_data$Injection_order),
  run_order = as.integer(source_data$Injection_order),
  assay = "XCMS-derived untargeted LC-MS plasma peak table",
  source_file = "Amide_data.rda",
  stringsAsFactors = FALSE
)
if (!identical(colnames(intensity), metadata$sample_id)) stop("WaveICA matrix/metadata alignment failed.")

detected_fraction <- rowMeans(is.finite(intensity))
feature_sd <- apply(intensity, 1, stats::sd)
keep <- detected_fraction >= 0.80 & is.finite(feature_sd) & feature_sd > 0
feature_metadata <- data.frame(
  feature_id = rownames(intensity),
  original_feature_id = feature_source_ids,
  original_column_index = seq_along(feature_source_ids) + 3L,
  mz = NA_real_, retention_time = NA_real_,
  detected_fraction = detected_fraction,
  missing_fraction = 1 - detected_fraction,
  retained = keep,
  removal_reason = ifelse(detected_fraction < 0.80, "observed_in_less_than_80_percent_of_injections", ifelse(!is.finite(feature_sd) | feature_sd <= 0, "zero_or_undefined_variance", "retained")),
  annotation_note = "The public object exposes Var1...Var6461 without m/z or retention-time annotation.",
  stringsAsFactors = FALSE
)
intensity <- intensity[keep, , drop = FALSE]
if (!all(is.finite(intensity)) || any(apply(intensity, 1, stats::sd) <= 0)) stop("Invalid values or zero-variance features remain in WaveICA.")

saveRDS(intensity, outputs[["matrix"]], compress = "xz")
write.csv(metadata, outputs[["metadata"]], row.names = FALSE, quote = TRUE, na = "")
write.csv(feature_metadata, outputs[["features"]], row.names = FALSE, quote = TRUE, na = "")
write.csv(feature_metadata[!feature_metadata$retained, , drop = FALSE], outputs[["removed"]], row.names = FALSE, quote = TRUE, na = "")

study <- metadata[!metadata$is_qc, , drop = FALSE]
group_batch <- table(study$batch, study$biological_group)
group_test <- suppressWarnings(stats::chisq.test(group_batch, correct = FALSE))
randomization <- list(
  data.frame(diagnostic = "group_vs_batch_chisq", batch = "all", estimate = unname(group_test$statistic), p_value = group_test$p.value, detail = paste(capture.output(print(group_batch)), collapse = " ")),
  data.frame(diagnostic = "group_vs_batch_cramers_v", batch = "all", estimate = cramers_v(group_batch), p_value = NA_real_, detail = "study samples only")
)
for (batch_id in sort(unique(study$batch))) {
  d <- study[study$batch == batch_id, , drop = FALSE]
  d <- d[order(d$within_batch_order), , drop = FALSE]
  run_test <- binary_runs_test(d$biological_group)
  randomization[[length(randomization) + 1L]] <- data.frame(diagnostic = "within_batch_runs_test", batch = batch_id, estimate = unname(run_test["z"]), p_value = unname(run_test["p_value"]), detail = sprintf("runs=%d; expected=%.2f", as.integer(run_test["runs"]), run_test["expected"]))
  randomization[[length(randomization) + 1L]] <- data.frame(diagnostic = "group_vs_within_order_rank_biserial", batch = batch_id, estimate = rank_biserial_order(d$biological_group, d$within_batch_order), p_value = NA_real_, detail = "group_1 versus group_0")
}
logistic <- stats::glm(I(biological_group == "group_1") ~ factor(batch) + within_batch_order, data = study, family = stats::binomial())
logistic_table <- as.data.frame(summary(logistic)$coefficients)
logistic_table$term <- rownames(logistic_table)
rownames(logistic_table) <- NULL
names(logistic_table)[1:4] <- c("estimate", "standard_error", "z_value", "p_value")
write.csv(do.call(rbind, randomization), file.path(audit_dir, "randomization_diagnostics.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(as.data.frame.matrix(group_batch), file.path(audit_dir, "group_by_batch.csv"), quote = TRUE)
write.csv(logistic_table, file.path(audit_dir, "group_logistic_batch_order.csv"), row.names = FALSE, quote = TRUE)
write.csv(qc_spacing_by_batch(metadata), file.path(audit_dir, "qc_spacing_by_batch.csv"), row.names = FALSE, quote = TRUE)
write.csv(metadata[, c("sample_id", "sample_type", "biological_group", "disease_label", "batch", "within_batch_order", "global_run_order")], file.path(audit_dir, "acquisition_sequence.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(histogram_source(intensity, "WaveICA adenocarcinoma"), file.path(audit_dir, "intensity_histogram_source.csv"), row.names = FALSE, quote = TRUE)

manifest <- read.csv(file.path(repo_root, "data", "public", "waveica_adenocarcinoma", "manifests", "download_manifest.csv"), stringsAsFactors = FALSE)
audit <- intensity_audit(intensity)
summary_value <- list(
  dataset = "WaveICA 2.0 human plasma adenocarcinoma demonstration dataset",
  source_repository = unique(manifest$source_repository), source_commit = unique(manifest$source_commit),
  source_files = manifest[, c("source_url", "local_path", "sha256", "role")],
  source_object = list(name = "Amide_data", dimensions = dim(source_data), metadata_columns = names(source_data)[1:3], feature_columns = length(feature_source_ids)),
  input = audit,
  preprocessing = list(
    orientation = "transposed from deposited injections_by_features to features_by_injections",
    retained_features = nrow(intensity), removed_features = sum(!keep),
    missingness_rule = "retain features observed in at least 80% of injections",
    imputation = "not applied because no values were missing",
    imputation_seed_reserved = 42,
    zero_handling = "48,820 deposited zeros retained as measured/deposited values; no source documentation establishes that they are missing",
    correction_or_normalization = "none"
  ),
  design = list(
    total_injections = nrow(metadata), study_injections = sum(!metadata$is_qc), qc_injections = sum(metadata$is_qc),
    batches = as.list(table(metadata$batch)), batch_study = as.list(table(metadata$batch[!metadata$is_qc])), batch_qc = as.list(table(metadata$batch[metadata$is_qc])),
    group_counts = as.list(table(metadata$biological_group, useNA = "no")),
    disease_mapping = "not assigned; deposited 0/1 counts match published disease totals but do not prove label direction",
    order_evidence = "Injection_order is globally unique 1..642, matches row order, and batches occupy contiguous intervals 1..217, 218..434, and 435..642"
  )
)
write_json_pretty(summary_value, outputs[["summary"]])

criteria <- data.frame(
  criterion_id = 1:14,
  criterion = c("Numerical feature-by-injection matrix", "Uncorrected or minimally processed intensities", "No prior batch/drift/QC correction", "Study and pooled QC in same matrix", "Unambiguous matrix-metadata matching", "Reliable supplied batch labels", "Reliable acquisition order", "Identifiable QCs", "Usable biological labels", "Sufficient study samples per batch", "Sufficient training QCs per batch", "Boundary QCs can remain in training", "No severe biological-batch/order confounding", "Safe input scale for WiNN"),
  status = c("PASS", "PASS WITH LIMITATION", "PASS WITH LIMITATION", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS WITH LIMITATION", "PASS", "PASS", "PASS", "PASS", "PASS"),
  evidence = c(
    "6461 numeric peak rows by 642 aligned injections after orientation.",
    "The paper/repository describe an XCMS-derived peak table; m/z and retention-time annotations and the complete XCMS script are not included.",
    "The package example treats Amide_data as original input to WaveICA; values are raw-intensity-like with zeros and strong batch/drift structure, with no evidence of prior correction.",
    "568 study samples and 74 pooled QCs coexist in Amide_data.",
    "Metadata and features are columns of the same deposited data frame; row names are unique sample IDs.",
    "Three supplied batches occupy contiguous acquisition intervals.",
    "Injection_order is globally unique 1..642 and equals deposited row order.",
    "group=QC identifies 25, 25, and 24 QCs by batch.",
    "Neutral groups are direct. Disease labels are not assigned: group counts match the two published cohort totals but do not prove which numeric value denotes which diagnosis.",
    "192, 192, and 184 study samples in batches 1, 2, and 3.",
    "Five hidden QCs per batch leave 20, 20, and 19 training QCs.",
    "The fixed holdout excludes the first and last QC of each batch.",
    sprintf("Both groups occur in every batch in virtually identical proportions; Cramer's V=%.4f and chi-squared p=%.4g.", cramers_v(group_batch), group_test$p.value),
    "All values are finite and nonnegative; no inverse transformation is required."
  ), stringsAsFactors = FALSE
)
write_suitability_csv(criteria, outputs[["report_csv"]])
report_lines <- c(
  "# WaveICA adenocarcinoma suitability audit",
  "",
  "**Verdict: SUITABLE WITH LIMITATIONS.** `Amide_data` provides an aligned XCMS-derived, pre-correction intensity matrix with global injection order, three batches, pooled QCs, and two balanced-within-batch biological groups. The benchmark may proceed.",
  "",
  "## Source and provenance",
  "",
  sprintf("Pinned source: `%s` at commit `%s`; object `Amide_data`.", unique(manifest$source_repository), unique(manifest$source_commit)),
  "Only the 33 MB R data object and small package source/documentation were downloaded. No mzXML, mzML, vendor, or raw spectra were downloaded.",
  "",
  "## Design and scale",
  "",
  sprintf("- Matrix: %d peaks × %d injections; %d study samples and %d QCs.", nrow(intensity), ncol(intensity), sum(!metadata$is_qc), sum(metadata$is_qc)),
  "- Batches: 192 study + 25 QC, 192 + 25, and 184 + 24.",
  "- Study groups by batch: group 1 = 168/168/161; group 0 = 24/24/23. All batches contain both groups in essentially identical proportions.",
  "- Injection order is globally unique 1–642 and batches are contiguous (1–217, 218–434, 435–642).",
  sprintf("- Intensities range from %.3g to %.3g (median %.3g); %.3f%% of entries are zero and no values are missing or nonfinite.", audit$minimum, audit$maximum, audit$median, 100 * audit$zero_fraction),
  "- Zeros are retained because the deposit does not establish that zero is an NA placeholder. No missingness filtering or imputation was needed.",
  "",
  "## Pre-correction evidence and label mapping",
  "",
  "The package source labels `Amide_data` as original data, removes only its first three metadata columns, and passes the remaining peak table directly to WaveICA 2.0. The associated paper describes XCMS preprocessing followed by batch-effect correction. Raw-intensity magnitude, zeros, and visible batch/run-order structure are consistent with this provenance; there is no evidence of prior drift or batch correction.",
  "The object stores groups only as `0` and `1`. It contains exactly 497 group-1 and 71 group-0 study samples, while the publication reports cohort totals of 497 colorectal-cancer and 71 chronic-enteritis patients. This count agreement does not prove the numeric label direction, so disease names are not assigned and all analyses use the neutral deposited groups.",
  "",
  "## Limitations",
  "",
  sprintf("The source lacks m/z/retention-time annotations, a complete XCMS parameter record, and a source-backed mapping from numeric groups to diagnoses. Group-batch Cramer's V is %.4f (p=%.4g); sequence diagnostics and logistic models are saved under `data/public/waveica_adenocarcinoma/audit/`.", cramers_v(group_batch), group_test$p.value),
  "Patients are not repeated across batches, so no replicate correlation or patient ICC will be calculated.",
  "",
  "## Criterion ledger",
  "",
  paste0("- ", criteria$criterion_id, ". **", criteria$status, "** — ", criteria$criterion, ": ", criteria$evidence)
)
writeLines(report_lines, outputs[["report_md"]], useBytes = TRUE)
message("WaveICA suitability audit and preprocessing completed: SUITABLE WITH LIMITATIONS.")
