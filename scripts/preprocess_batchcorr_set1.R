#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript or R --file.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args

required_packages <- c("jsonlite", "missForest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Use the R 4.5 benchmark runtime documented in README.md."
  )
}

source_path <- file.path(repo_root, "data", "public", "batchcorr_set1", "source", "BC.RData")
processed_dir <- file.path(repo_root, "data", "public", "batchcorr_set1", "processed")
audit_dir <- file.path(repo_root, "data", "public", "batchcorr_set1", "audit")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

output_paths <- c(
  metadata = file.path(processed_dir, "BatchCorr_set1_metadata.csv"),
  unfiltered = file.path(processed_dir, "BatchCorr_set1_data_unfiltered_log.csv"),
  filtered = file.path(processed_dir, "BatchCorr_set1_data_filtered_log.csv"),
  imputed_log = file.path(processed_dir, "BatchCorr_set1_imputed_log.csv"),
  imputed_linear = file.path(processed_dir, "BatchCorr_set1_imputed_data.csv"),
  feature_metadata = file.path(processed_dir, "BatchCorr_set1_feature_metadata.csv"),
  removed_features = file.path(processed_dir, "BatchCorr_set1_removed_features.csv"),
  summary = file.path(processed_dir, "BatchCorr_set1_preprocessing_summary.json")
)

if (!force && all(file.exists(output_paths))) {
  message("All BatchCorr Set 1 preprocessing outputs already exist; reusing them.")
  quit(save = "no", status = 0)
}

if (!file.exists(source_path)) {
  stop("Missing ", source_path, ". Run scripts/download_batchcorr_set1.py first.")
}

write_matrix_csv <- function(mat, path) {
  out <- data.frame(feature_id = rownames(mat), mat, check.names = FALSE)
  write.csv(out, path, row.names = FALSE, quote = TRUE, na = "NA")
}

sha256_file <- function(path) {
  unname(tools::md5sum(path))
}

loaded <- new.env(parent = emptyenv())
load(source_path, envir = loaded)
if (!all(c("set.1", "set.1.Y") %in% ls(loaded))) {
  stop("BC.RData does not contain both set.1 and set.1.Y.")
}

x_source <- loaded$set.1
y_source <- loaded$set.1.Y
if (!is.matrix(x_source) || !is.numeric(x_source)) stop("set.1 is not a numeric matrix.")
if (!is.data.frame(y_source)) stop("set.1.Y is not a data frame.")
if (!identical(rownames(x_source), rownames(y_source))) {
  stop("set.1 rows do not align exactly with set.1.Y rows.")
}
if (!all(c("SeqNr", "Batch", "SCode") %in% names(y_source))) {
  stop("set.1.Y is missing SeqNr, Batch, or SCode.")
}
if (anyNA(y_source$SeqNr) || anyNA(y_source$Batch) || anyNA(y_source$SCode)) {
  stop("set.1.Y contains missing sequence, batch, or sample-code values.")
}
if (anyDuplicated(y_source$SeqNr)) stop("SeqNr is not globally unique.")

sample_id <- rownames(y_source)
is_qc <- as.character(y_source$SCode) == "ref"
batch_chr <- as.character(y_source$Batch)
run_order <- as.integer(y_source$SeqNr)
batch_first_order <- tapply(run_order, batch_chr, min)
chronological_batches <- names(sort(batch_first_order))
batch_order_index <- match(batch_chr, chronological_batches)
within_batch_order <- integer(length(run_order))
for (batch_id in unique(batch_chr)) {
  idx <- which(batch_chr == batch_id)
  within_batch_order[idx] <- rank(run_order[idx], ties.method = "first")
}

accession_id <- ifelse(is_qc, NA_character_, as.character(y_source$SCode))
biological_replicate <- rep(NA_integer_, length(sample_id))
for (accession in unique(accession_id[!is.na(accession_id)])) {
  idx <- which(accession_id == accession)
  idx <- idx[order(run_order[idx])]
  biological_replicate[idx] <- seq_along(idx)
}

metadata <- data.frame(
  sample_id = sample_id,
  original_sample_name = sample_id,
  sample_type = ifelse(is_qc, "QC", "Study"),
  class = ifelse(is_qc, "QC", "Sample"),
  accession_id = accession_id,
  genotype = accession_id,
  biological_replicate = biological_replicate,
  replicate_group = accession_id,
  batch = batch_chr,
  batch_order_index = batch_order_index,
  run_order = run_order,
  global_run_order = run_order,
  within_batch_order = within_batch_order,
  is_qc = is_qc,
  is_reference = is_qc,
  is_study_sample = !is_qc,
  assay = "untargeted LC-MS",
  polarity = "negative",
  source_object = "set.1 / set.1.Y",
  stringsAsFactors = FALSE
)
metadata <- metadata[order(metadata$run_order), , drop = FALSE]

feature_ids <- sprintf("Feature_%04d", seq_len(ncol(x_source)))
original_indices <- seq_len(ncol(x_source))
detected_count <- colSums(!is.na(x_source))
detected_fraction <- colMeans(!is.na(x_source))
feature_variance <- apply(x_source, 2, stats::var, na.rm = TRUE)
keep <- detected_fraction >= 0.80 & is.finite(feature_variance) & feature_variance > 0
removal_reason <- ifelse(
  detected_fraction < 0.80,
  "observed_in_less_than_80_percent_of_injections",
  ifelse(!is.finite(feature_variance) | feature_variance <= 0, "zero_or_undefined_variance", "retained")
)

feature_metadata <- data.frame(
  feature_id = feature_ids,
  original_column_index = original_indices,
  original_feature_name = NA_character_,
  detected_count = detected_count,
  detected_fraction = detected_fraction,
  missing_fraction = 1 - detected_fraction,
  variance_on_deposited_log_scale = feature_variance,
  retained = keep,
  removal_reason = removal_reason,
  note = "The anonymized package matrix has no column names or m/z/retention-time annotations.",
  stringsAsFactors = FALSE
)

x_feature_by_sample <- t(x_source)
rownames(x_feature_by_sample) <- feature_ids
colnames(x_feature_by_sample) <- sample_id
x_filtered_log <- x_feature_by_sample[keep, metadata$sample_id, drop = FALSE]

set.seed(42)
if (anyNA(x_filtered_log)) {
  missforest_input <- as.data.frame(t(x_filtered_log), check.names = FALSE)
  imputation <- missForest::missForest(
    missforest_input,
    maxiter = 10,
    ntree = 100,
    verbose = TRUE,
    parallelize = "no"
  )
  x_imputed_log <- t(as.matrix(imputation$ximp))
  oob_error <- unname(imputation$OOBerror)
} else {
  x_imputed_log <- x_filtered_log
  oob_error <- NA_real_
}
rownames(x_imputed_log) <- rownames(x_filtered_log)
colnames(x_imputed_log) <- colnames(x_filtered_log)

# The package demo explicitly states that set.1 is already log-scaled but does
# not state the logarithm base. expm1 is therefore used as a reversible
# canonical linearization for multiplicative correction methods. Every metric
# is evaluated after log1p, which recovers the deposited/imputed analysis scale.
x_imputed_linear <- expm1(x_imputed_log)
if (any(!is.finite(x_imputed_linear)) || any(x_imputed_linear < 0)) {
  stop("Canonical expm1 linearization produced invalid values.")
}

if (anyNA(x_imputed_log) || anyNA(x_imputed_linear)) stop("Imputation left missing values.")
if (anyDuplicated(metadata$sample_id)) stop("Duplicate sample IDs remain.")
if (!identical(colnames(x_imputed_linear), metadata$sample_id)) stop("Matrix and metadata orders differ.")
if (any(apply(x_imputed_linear, 1, stats::sd) <= 0)) stop("Zero-variance features remain.")

write.csv(metadata, output_paths[["metadata"]], row.names = FALSE, quote = TRUE, na = "")
write_matrix_csv(x_feature_by_sample, output_paths[["unfiltered"]])
write_matrix_csv(x_filtered_log, output_paths[["filtered"]])
write_matrix_csv(x_imputed_log, output_paths[["imputed_log"]])
write_matrix_csv(x_imputed_linear, output_paths[["imputed_linear"]])
write.csv(feature_metadata, output_paths[["feature_metadata"]], row.names = FALSE, quote = TRUE, na = "")
write.csv(feature_metadata[!feature_metadata$retained, , drop = FALSE], output_paths[["removed_features"]], row.names = FALSE, quote = TRUE, na = "")

batch_summary <- do.call(rbind, lapply(chronological_batches, function(batch_id) {
  d <- metadata[metadata$batch == batch_id, , drop = FALSE]
  data.frame(
    batch = batch_id,
    chronological_index = unique(d$batch_order_index),
    n_total = nrow(d),
    n_study = sum(d$is_study_sample),
    n_qc = sum(d$is_qc),
    first_global_order = min(d$run_order),
    last_global_order = max(d$run_order),
    sequence_span = max(d$run_order) - min(d$run_order) + 1L,
    omitted_sequence_numbers = max(d$run_order) - min(d$run_order) + 1L - nrow(d),
    first_qc_order = min(d$run_order[d$is_qc]),
    last_qc_order = max(d$run_order[d$is_qc]),
    max_study_injections_between_qcs = {
      q <- sort(d$within_batch_order[d$is_qc])
      max(diff(c(0L, q, nrow(d) + 1L)) - 1L)
    },
    stringsAsFactors = FALSE
  )
}))
write.csv(batch_summary, file.path(audit_dir, "batch_and_qc_counts.csv"), row.names = FALSE, quote = TRUE)

qc_positions <- metadata[metadata$is_qc, c("sample_id", "batch", "run_order", "within_batch_order")]
write.csv(qc_positions, file.path(audit_dir, "qc_positions.csv"), row.names = FALSE, quote = TRUE)

study_metadata <- metadata[metadata$is_study_sample, , drop = FALSE]
accession_groups <- split(study_metadata, study_metadata$accession_id)
accession_balance <- do.call(rbind, lapply(names(accession_groups), function(accession) {
  d <- accession_groups[[accession]]
  data.frame(
    accession_id = accession,
    n_biological_replicates = nrow(d),
    n_batches = length(unique(d$batch)),
    batches = paste(sort(unique(d$batch)), collapse = ";"),
    cross_batch_replicates = length(unique(d$batch)) > 1L,
    stringsAsFactors = FALSE
  )
}))
write.csv(accession_balance, file.path(audit_dir, "accession_batch_balance.csv"), row.names = FALSE, quote = TRUE)

source_manifest <- read.csv(
  file.path(repo_root, "data", "public", "batchcorr_set1", "source", "download_manifest.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
source_sha <- source_manifest$sha256[source_manifest$local_path == "data/public/batchcorr_set1/source/BC.RData"]

summary <- list(
  source = list(
    repository = "rwehrens/BatchCorrMetabolomics",
    commit = "e0c7668140e206dcdae9afa602dd2e1b337ac4f6",
    object = "set.1 and set.1.Y",
    bc_rdata_sha256 = if (length(source_sha)) source_sha[[1]] else NA_character_
  ),
  deposited_matrix = list(
    samples = nrow(x_source),
    features = ncol(x_source),
    missing_values = sum(is.na(x_source)),
    missing_fraction = mean(is.na(x_source)),
    zeros = sum(x_source == 0, na.rm = TRUE),
    scale = "already log-scaled according to the package demo; base not documented"
  ),
  filtering = list(
    policy = "retain features observed in at least 80% of all injections and remove zero-variance features",
    retained_features = sum(keep),
    removed_features = sum(!keep),
    retained_missing_fraction_before_imputation = mean(is.na(x_filtered_log))
  ),
  imputation = list(
    method = "missForest",
    seed = 42,
    maxiter = 10,
    ntree = 100,
    oob_error = oob_error,
    scale = "deposited log scale"
  ),
  analysis_matrix = list(
    features = nrow(x_imputed_linear),
    injections = ncol(x_imputed_linear),
    transformation = "expm1 canonical linearization; log1p is the exact evaluation inverse",
    minimum = min(x_imputed_linear),
    maximum = max(x_imputed_linear)
  ),
  design = list(
    study_injections = sum(metadata$is_study_sample),
    qc_injections = sum(metadata$is_qc),
    batches = length(unique(metadata$batch)),
    accessions = length(unique(stats::na.omit(metadata$accession_id))),
    accessions_with_at_least_two_replicates = sum(accession_balance$n_biological_replicates >= 2L),
    accessions_with_cross_batch_replicates = sum(accession_balance$cross_batch_replicates),
    chronological_batch_order = chronological_batches
  ),
  zero_handling = "No zeros occur in set.1. NAs are documented non-detects; no zero-to-NA conversion was performed."
)
jsonlite::write_json(summary, output_paths[["summary"]], pretty = TRUE, auto_unbox = TRUE, na = "null")

session_lines <- capture.output(sessionInfo())
writeLines(session_lines, file.path(processed_dir, "preprocessing_sessionInfo.txt"))
message(
  "Completed BatchCorr Set 1 preprocessing: ", nrow(x_imputed_linear), " features x ",
  ncol(x_imputed_linear), " injections (", sum(metadata$is_qc), " QCs)."
)
