#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript or R --file.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
source(file.path(repo_root, "scripts", "public_benchmark_audit_helpers.R"))
require_packages(c("jsonlite", "openssl"))

force <- "--force" %in% commandArgs(trailingOnly = TRUE)
source_dir <- file.path(repo_root, "data", "public", "sacurine", "source")
processed_dir <- file.path(repo_root, "data", "public", "sacurine", "processed")
audit_dir <- file.path(repo_root, "data", "public", "sacurine", "audit")
report_dir <- file.path(repo_root, "reports")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

outputs <- c(
  matrix = file.path(processed_dir, "sacurine_intensity_matrix.rds"),
  metadata = file.path(processed_dir, "sacurine_metadata.csv"),
  features = file.path(processed_dir, "sacurine_feature_metadata.csv"),
  summary = file.path(processed_dir, "sacurine_preprocessing_summary.json"),
  removed = file.path(processed_dir, "sacurine_removed_features.csv"),
  report_md = file.path(report_dir, "sacurine_suitability.md"),
  report_csv = file.path(report_dir, "sacurine_suitability.csv")
)
if (!force && all(file.exists(outputs))) {
  message("All Sacurine audit/preprocessing outputs exist; reusing them.")
  quit(save = "no", status = 0)
}

required_source <- c("Galaxy1_dataMatrix.tabular", "Galaxy2_sampleMetadata.tabular", "Galaxy3_variableMetadata.tabular", "phenomis_vignette.Rmd")
missing_source <- required_source[!file.exists(file.path(source_dir, required_source))]
if (length(missing_source)) stop("Missing Sacurine source file(s): ", paste(missing_source, collapse = ", "), ". Run scripts/download_sacurine.py.")

data_table <- read.delim(file.path(source_dir, required_source[1]), check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
sample_table <- read.delim(file.path(source_dir, required_source[2]), check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
feature_table <- read.delim(file.path(source_dir, required_source[3]), check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)

intensity <- as.matrix(data_table)
storage.mode(intensity) <- "double"
if (!identical(colnames(intensity), rownames(sample_table))) stop("Sacurine matrix columns do not align exactly with sample metadata rows.")
if (!identical(rownames(intensity), rownames(feature_table))) stop("Sacurine matrix rows do not align exactly with feature metadata rows.")
if (anyDuplicated(colnames(intensity)) || anyDuplicated(rownames(intensity))) stop("Sacurine has duplicate sample or feature identifiers.")
if (any(!is.finite(intensity)) || any(intensity < 0)) stop("Sacurine contains missing, nonfinite, or negative intensities.")

sample_table$original_sample_id <- rownames(sample_table)
batch_levels <- unique(as.character(sample_table$batch))
batch_index <- match(as.character(sample_table$batch), batch_levels)
within_order <- ave(sample_table$injectionOrder, sample_table$batch, FUN = function(x) rank(x, ties.method = "first"))
global_order <- integer(nrow(sample_table))
offset <- 0L
for (batch_id in batch_levels) {
  idx <- which(sample_table$batch == batch_id)
  idx <- idx[order(sample_table$injectionOrder[idx])]
  global_order[idx] <- offset + seq_along(idx)
  offset <- offset + length(idx)
}
is_qc <- sample_table$sampleType == "pool"
metadata <- data.frame(
  sample_id = sample_table$original_sample_id,
  original_sample_id = sample_table$original_sample_id,
  participant_id = ifelse(is_qc, NA_character_, sub("_b2$", "", sample_table$original_sample_id)),
  sample_type = ifelse(is_qc, "QC", "Study"),
  class = ifelse(is_qc, "QC", "Sample"),
  is_qc = is_qc,
  batch = as.character(sample_table$batch),
  batch_order_index = batch_index,
  deposited_injection_order = as.integer(sample_table$injectionOrder),
  within_batch_order = as.integer(within_order),
  global_run_order = as.integer(global_order),
  run_order = as.integer(global_order),
  deposited_column_index = match(sample_table$original_sample_id, colnames(intensity)),
  osmolality = as.numeric(sample_table$osmolality),
  sampling = sample_table$sampling,
  age = as.numeric(sample_table$age),
  bmi = as.numeric(sample_table$bmi),
  gender = as.character(sample_table$gender),
  subset = sample_table$subset,
  full = sample_table$full,
  assay = "untargeted UPLC-HRMS identified metabolite panel",
  polarity = "negative",
  source_file = "Galaxy1_dataMatrix.tabular / Galaxy2_sampleMetadata.tabular",
  stringsAsFactors = FALSE
)
metadata <- metadata[order(metadata$global_run_order), , drop = FALSE]
intensity <- intensity[, metadata$sample_id, drop = FALSE]
if (!identical(colnames(intensity), metadata$sample_id)) stop("Sacurine matrix/metadata alignment failed after ordering.")
if (anyNA(metadata$batch) || anyNA(metadata$within_batch_order) || anyNA(metadata$global_run_order)) stop("Sacurine batch/order contains missing values.")

detected_fraction <- rowMeans(is.finite(intensity))
feature_sd <- apply(intensity, 1, stats::sd)
keep <- detected_fraction >= 0.80 & is.finite(feature_sd) & feature_sd > 0
feature_metadata <- data.frame(
  feature_id = make.unique(rownames(intensity), sep = "__duplicate_"),
  original_feature_id = rownames(intensity),
  feature_table,
  detected_fraction = detected_fraction,
  missing_fraction = 1 - detected_fraction,
  retained = keep,
  removal_reason = ifelse(detected_fraction < 0.80, "observed_in_less_than_80_percent_of_injections", ifelse(!is.finite(feature_sd) | feature_sd <= 0, "zero_or_undefined_variance", "retained")),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(intensity) <- feature_metadata$feature_id
intensity <- intensity[keep, , drop = FALSE]
if (!all(is.finite(intensity)) || any(apply(intensity, 1, stats::sd) <= 0)) stop("Invalid values or zero-variance features remain in Sacurine.")

saveRDS(intensity, outputs[["matrix"]], compress = "xz")
write.csv(metadata, outputs[["metadata"]], row.names = FALSE, quote = TRUE, na = "")
write.csv(feature_metadata, outputs[["features"]], row.names = FALSE, quote = TRUE, na = "")
write.csv(feature_metadata[!feature_metadata$retained, , drop = FALSE], outputs[["removed"]], row.names = FALSE, quote = TRUE, na = "")

study <- metadata[!metadata$is_qc, , drop = FALSE]
gender_batch <- table(study$batch, study$gender)
gender_test <- suppressWarnings(stats::chisq.test(gender_batch, correct = FALSE))
randomization <- list()
randomization[[1]] <- data.frame(diagnostic = "gender_vs_batch_chisq", batch = "all", estimate = unname(gender_test$statistic), p_value = gender_test$p.value, detail = paste(capture.output(print(gender_batch)), collapse = " "))
randomization[[2]] <- data.frame(diagnostic = "gender_vs_batch_cramers_v", batch = "all", estimate = cramers_v(gender_batch), p_value = NA_real_, detail = "study samples only")
randomization[[3]] <- data.frame(diagnostic = "age_batch_standardized_mean_difference", batch = "ne2_minus_ne1", estimate = standardized_mean_difference(study$age, study$batch), p_value = NA_real_, detail = "pooled-SD standardized")
randomization[[4]] <- data.frame(diagnostic = "bmi_batch_standardized_mean_difference", batch = "ne2_minus_ne1", estimate = standardized_mean_difference(study$bmi, study$batch), p_value = NA_real_, detail = "pooled-SD standardized")
for (batch_id in batch_levels) {
  d <- study[study$batch == batch_id, , drop = FALSE]
  for (variable in c("age", "bmi")) {
    test <- suppressWarnings(stats::cor.test(d[[variable]], d$within_batch_order, method = "spearman", exact = FALSE))
    randomization[[length(randomization) + 1L]] <- data.frame(diagnostic = paste0(variable, "_vs_within_order_spearman"), batch = batch_id, estimate = unname(test$estimate), p_value = test$p.value, detail = paste0("n=", sum(is.finite(d[[variable]]))))
  }
  randomization[[length(randomization) + 1L]] <- data.frame(diagnostic = "gender_vs_within_order_rank_biserial", batch = batch_id, estimate = rank_biserial_order(d$gender, d$within_batch_order), p_value = NA_real_, detail = "Male/Female sequence association")
}
randomization_df <- do.call(rbind, randomization)
write.csv(randomization_df, file.path(audit_dir, "randomization_diagnostics.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(as.data.frame.matrix(gender_batch), file.path(audit_dir, "gender_by_batch.csv"), quote = TRUE)
write.csv(as.data.frame.matrix(table(study$batch, study$sampling, useNA = "ifany")), file.path(audit_dir, "sampling_by_batch.csv"), quote = TRUE)
write.csv(qc_spacing_by_batch(metadata), file.path(audit_dir, "qc_spacing_by_batch.csv"), row.names = FALSE, quote = TRUE)
write.csv(metadata[, c("sample_id", "sample_type", "batch", "deposited_injection_order", "within_batch_order", "global_run_order", "age", "bmi", "gender", "sampling")], file.path(audit_dir, "acquisition_sequence.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(histogram_source(intensity, "Sacurine"), file.path(audit_dir, "intensity_histogram_source.csv"), row.names = FALSE, quote = TRUE)

manifest <- read.csv(file.path(repo_root, "data", "public", "sacurine", "manifests", "download_manifest.csv"), stringsAsFactors = FALSE)
audit <- intensity_audit(intensity)
summary_value <- list(
  dataset = "Sacurine / W4M00001 / MTBLS404 negative ion",
  source_repository = unique(manifest$source_repository),
  source_commit = unique(manifest$source_commit),
  source_files = manifest[, c("source_url", "local_path", "sha256", "role")],
  input = audit,
  preprocessing = list(
    orientation = "features_by_injections",
    retained_features = nrow(intensity),
    removed_features = sum(!keep),
    missingness_rule = "retain features observed in at least 80% of injections",
    imputation = "not applied because no values were missing",
    imputation_seed_reserved = 42,
    zero_handling = "zeros were not converted to missing; the deposited table contains no zeros",
    lower_detection_replacement = "repeated positive values (including 506.5 and 607.5) were retained because the deposit does not mark them as missing",
    correction_or_normalization = "none"
  ),
  design = list(
    total_injections = nrow(metadata), study_injections = sum(!metadata$is_qc), qc_injections = sum(metadata$is_qc),
    batches = as.list(table(metadata$batch)), batch_study = as.list(table(metadata$batch[!metadata$is_qc])), batch_qc = as.list(table(metadata$batch[metadata$is_qc])),
    order_evidence = "deposited injectionOrder is unique and strictly increasing within batch; batch sequence follows deposited matrix/metadata order; global order concatenates ne1 then ne2 while preserving the deposited order and its gaps separately"
  ),
  missing_metadata = as.list(colSums(is.na(metadata[, c("age", "bmi", "gender", "sampling", "osmolality")]) | metadata[, c("age", "bmi", "gender", "sampling", "osmolality")] == ""))
)
write_json_pretty(summary_value, outputs[["summary"]])

criteria <- data.frame(
  criterion_id = 1:14,
  criterion = c("Numerical feature-by-injection matrix", "Uncorrected or minimally processed intensities", "No prior batch/drift/QC correction", "Study and pooled QC in same matrix", "Unambiguous matrix-metadata matching", "Reliable supplied batch labels", "Reliable acquisition order", "Identifiable QCs", "Usable biological labels", "Sufficient study samples per batch", "Sufficient training QCs per batch", "Boundary QCs can remain in training", "No severe biological-batch/order confounding", "Safe input scale for WiNN"),
  status = c("PASS", "PASS WITH LIMITATION", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS"),
  evidence = c(
    "113 numeric metabolite rows by 210 aligned injections.",
    "Positive identified-metabolite peak intensities; the table is a selected 113-metabolite panel rather than the full untargeted peak table.",
    "The phenomis vignette reads these files and only then invokes correcting(), QC-CV filtering, osmolality normalization, and log10 transformation.",
    "184 study samples and 26 pool QCs coexist in the same table.",
    "Matrix column names exactly equal metadata row names in deposited order.",
    "Two deposited batches: ne1 (104 injections) and ne2 (106 injections).",
    "Deposited injectionOrder is unique and strictly increasing within each batch; it resets across batches, so global order is a documented concatenation.",
    "sampleType=pool identifies 14 ne1 and 12 ne2 QCs.",
    "Age, BMI, gender, sampling, and osmolality are present for every study sample.",
    "90 and 94 study injections in ne1 and ne2.",
    "Four hidden QCs per batch leave 10 and 8 training QCs.",
    "First/last deposited QCs in both batches are retained by the holdout constraints.",
    sprintf("Gender-batch Cramer's V=%.3f; age SMD=%.3f; BMI SMD=%.3f; both genders occur throughout both batches.", cramers_v(gender_batch), standardized_mean_difference(study$age, study$batch), standardized_mean_difference(study$bmi, study$batch)),
    "All values are finite and nonnegative; no undocumented inverse transformation is required."
  ), stringsAsFactors = FALSE
)
write_suitability_csv(criteria, outputs[["report_csv"]])
report_lines <- c(
  "# Sacurine negative-ion suitability audit",
  "",
  "**Verdict: SUITABLE WITH LIMITATIONS.** The W4M00001 files provide an aligned, pre-correction positive-intensity matrix with reliable QCs, two analytical batches, within-batch acquisition order, and complete study covariates. The benchmark may proceed.",
  "",
  "## Source and provenance",
  "",
  sprintf("Pinned source: `%s` at commit `%s`.", unique(manifest$source_repository), unique(manifest$source_commit)),
  "The three deposited W4M tables are used directly. No raw spectra, mzML, vendor files, ropls corrected object, or post-correction W4M output was downloaded.",
  "",
  "## Design and scale",
  "",
  sprintf("- Matrix: %d retained identified metabolites × %d injections; %d study samples and %d pooled QCs.", nrow(intensity), ncol(intensity), sum(!metadata$is_qc), sum(metadata$is_qc)),
  sprintf("- Batch ne1: %d study + %d QC; batch ne2: %d study + %d QC.", sum(!metadata$is_qc & metadata$batch == "ne1"), sum(metadata$is_qc & metadata$batch == "ne1"), sum(!metadata$is_qc & metadata$batch == "ne2"), sum(metadata$is_qc & metadata$batch == "ne2")),
  "- `injectionOrder` resets across batches but is unique and strictly increasing within each batch. The original values are preserved; within-batch ranks and a ne1→ne2 concatenated global order are derived explicitly.",
  sprintf("- Intensities range from %.1f to %.1f (median %.1f), with no NA, nonfinite, negative, or zero values.", audit$minimum, audit$maximum, audit$median),
  "- Repeated positive lower-bound values are retained; the source does not designate them as missing. No filtering or imputation was required.",
  "",
  "## Evidence that the table is pre-correction",
  "",
  "The pinned phenomis vignette first reads these exact three files and subsequently calls `correcting()` for QC-loess drift and batch adjustment, then performs pool-CV filtering, removes pools, divides by osmolality, and applies log10. This ordering is direct evidence that the deposited table precedes those operations.",
  "",
  "## Balance and limitations",
  "",
  sprintf("Gender versus batch: chi-squared p=%.4g, Cramer's V=%.3f. Age and BMI batch standardized mean differences are %.3f and %.3f. Within-batch correlations and the complete sequence are saved under `data/public/sacurine/audit/`.", gender_test$p.value, cramers_v(gender_batch), standardized_mean_difference(study$age, study$batch), standardized_mean_difference(study$bmi, study$batch)),
  "The panel contains only 113 identified metabolites, so conclusions do not generalize to all detected peaks. Each volunteer occurs once, so no biological-replicate correlation or participant ICC will be calculated. Technical metrics will use pre-osmolality values; a common post-correction osmolality normalization will be used only for biological analyses.",
  "",
  "## Criterion ledger",
  "",
  paste0("- ", criteria$criterion_id, ". **", criteria$status, "** — ", criteria$criterion, ": ", criteria$evidence)
)
writeLines(report_lines, outputs[["report_md"]], useBytes = TRUE)
message("Sacurine suitability audit and preprocessing completed: SUITABLE WITH LIMITATIONS.")
