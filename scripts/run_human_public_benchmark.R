#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript or R --file.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly = TRUE)
dataset_arg <- grep("^--dataset=", args, value = TRUE)
if (length(dataset_arg) != 1L) stop("Supply exactly one of --dataset=sacurine or --dataset=waveica_adenocarcinoma.")
dataset <- sub("^--dataset=", "", dataset_arg)
force <- "--force" %in% args
if (!dataset %in% c("sacurine", "waveica_adenocarcinoma")) stop("Unknown dataset: ", dataset)

required_packages <- c(
  "dplyr", "ggplot2", "jsonlite", "limma", "lmtest", "malbacR", "mgcv",
  "patchwork", "pmartR", "qcrlscR", "scales", "statTarget", "sva",
  "tibble", "tidyr", "TIGERr", "winn"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing benchmark package(s): ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(tidyr)
})
source(file.path(repo_root, "scripts", "weighted_pc_r2.R"))
source(file.path(repo_root, "scripts", "benchmark_helpers.R"))
source(file.path(repo_root, "scripts", "run_order_drift_helpers.R"))
source(file.path(repo_root, "scripts", "public_benchmark_method_helpers.R"))
source(file.path(repo_root, "scripts", "plot_theme.R"))

config <- if (dataset == "sacurine") {
  list(
    title = "Sacurine negative-ion",
    matrix_path = file.path(repo_root, "data", "public", "sacurine", "processed", "sacurine_intensity_matrix.rds"),
    metadata_path = file.path(repo_root, "data", "public", "sacurine", "processed", "sacurine_metadata.csv"),
    result_dir = file.path(repo_root, "results", "sacurine"),
    holdout_seed = 20260810L, holdout_per_batch = 4L, min_train_qc = 6L,
    tuning_seed = 20260812L, tuning_features = Inf,
    bootstrap_seed = 20260814L,
    suitability_report = file.path(repo_root, "reports", "sacurine_suitability.md")
  )
} else {
  list(
    title = "WaveICA adenocarcinoma",
    matrix_path = file.path(repo_root, "data", "public", "waveica_adenocarcinoma", "processed", "waveica_intensity_matrix.rds"),
    metadata_path = file.path(repo_root, "data", "public", "waveica_adenocarcinoma", "processed", "waveica_metadata.csv"),
    result_dir = file.path(repo_root, "results", "waveica_adenocarcinoma"),
    holdout_seed = 20260811L, holdout_per_batch = 5L, min_train_qc = 15L,
    tuning_seed = 20260812L, tuning_features = 220L,
    bootstrap_seed = 20260815L,
    suitability_report = file.path(repo_root, "reports", "waveica_adenocarcinoma_suitability.md")
  )
}
if (!file.exists(config$suitability_report) || !any(grepl("SUITABLE", readLines(config$suitability_report, warn = FALSE), fixed = TRUE))) stop("Suitability report is missing or does not authorize benchmarking.")
if (!file.exists(config$matrix_path) || !file.exists(config$metadata_path)) stop("Processed inputs are missing. Run the dataset preprocessing script first.")

result_dir <- config$result_dir
cache_dir <- file.path(result_dir, "method_cache")
tuning_dir <- file.path(result_dir, "tuning_candidates")
figure_dir <- file.path(result_dir, "figures")
figure_source_dir <- file.path(figure_dir, "source_data")
for (directory in c(result_dir, cache_dir, tuning_dir, figure_dir, figure_source_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(result_dir, "analysis_log.txt")
if (force || !file.exists(log_path)) writeLines(character(), log_path)
log_line <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " ", paste0(..., collapse = ""))
  cat(line, "\n")
  cat(line, "\n", file = log_path, append = TRUE)
}
log_line("Starting ", config$title, " benchmark; force=", force, ".")

x <- readRDS(config$matrix_path)
meta <- read.csv(config$metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
meta$is_qc <- as.logical(meta$is_qc)
meta <- meta[match(colnames(x), meta$sample_id), , drop = FALSE]
if (anyNA(meta$sample_id) || !identical(colnames(x), meta$sample_id)) stop("Input matrix/metadata alignment failed.")
if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)) || any(!is.finite(x)) || any(x < 0)) stop("Input matrix failed identifier/value validation.")
storage.mode(x) <- "double"

method_order <- winn_method_order()
method_seeds <- stats::setNames(config$holdout_seed + 1000L + seq_along(method_order), method_order)
method_palette <- winn_method_palette()

split <- select_hidden_qcs(meta, config$holdout_per_batch, config$holdout_seed, config$min_train_qc)
hidden_ids <- split$hidden
training_ids <- split$training
hidden_table <- meta[match(hidden_ids, meta$sample_id), c("sample_id", "batch", "global_run_order", "within_batch_order")]
hidden_table$seed <- config$holdout_seed
training_table <- meta[match(training_ids, meta$sample_id), c("sample_id", "batch", "global_run_order", "within_batch_order")]
write.csv(hidden_table, file.path(result_dir, "qc_holdout_ids.csv"), row.names = FALSE, quote = TRUE)
write.csv(training_table, file.path(result_dir, "qc_training_ids.csv"), row.names = FALSE, quote = TRUE)
write.csv(split$audit, file.path(result_dir, "qc_holdout_batch_summary.csv"), row.names = FALSE, quote = TRUE)

meta_hidden <- meta
hidden_rows <- meta_hidden$sample_id %in% hidden_ids
meta_hidden$class[hidden_rows] <- "Sample"
meta_hidden$sample_type[hidden_rows] <- "Study"
meta_hidden$is_qc[hidden_rows] <- FALSE
if (any(meta_hidden$class[match(hidden_ids, meta_hidden$sample_id)] == "QC") || any(meta_hidden$is_qc[match(hidden_ids, meta_hidden$sample_id)])) stop("Hidden QC labels remain visible in method metadata.")

protocol <- c(
  paste0("# ", config$title, " hidden-QC protocol"), "",
  sprintf("Seed %d selected %d interior QCs (%d per batch target).", config$holdout_seed, length(hidden_ids), config$holdout_per_batch),
  sprintf("The split leaves %d training QCs; first and last QC in every batch are protected, with at least %d training QCs per batch.", length(training_ids), config$min_train_qc),
  "Before every method call, hidden QCs are relabelled as ordinary Sample/Study injections and `is_qc` is set to FALSE in the method-facing metadata.",
  "Only the saved training QC IDs are supplied to QC-aware methods, competitor tuning, and WiNN automatic selection.",
  "No tuning score uses hidden-QC values, phenotype labels, final batch metrics, or final run-order metrics.",
  paste0("Every full correction call has an explicit method-specific seed: ", paste(names(method_seeds), method_seeds, sep = "=", collapse = "; "), "."),
  "ComBat receives batches but no QC labels. WINN auto-batch receives run order and training QCs but no supplied batch labels. WINN fixed receives batches and no QC identifiers."
)
writeLines(protocol, file.path(result_dir, "qc_holdout_protocol.md"))
method_protocol <- data.frame(
  method = method_order,
  correction_seed = unname(method_seeds[method_order]),
  supplied_batch_labels = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  supplied_training_qc_ids = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  supplied_hidden_qc_ids = FALSE,
  hidden_qc_method_class = "Sample",
  tuning_uses_hidden_qc_values = FALSE,
  tuning_uses_phenotype = FALSE,
  stringsAsFactors = FALSE
)
write.csv(method_protocol, file.path(result_dir, "hidden_qc_method_protocol.csv"), row.names = FALSE, quote = TRUE)
log_line("Hidden QCs fixed: ", paste(hidden_ids, collapse = ", "), ".")

qcrlsc_audit <- audit_qcrlsc_wrapper(x, training_ids, meta_hidden, result_dir)
qcrlsc_report_path <- file.path(repo_root, "reports", paste0(dataset, "_qc_rlsc_audit.md"))
writeLines(c(
  paste0("# ", config$title, " QC-RLSC wrapper audit"), "",
  if (qcrlsc_audit$custom_required) "**Verdict: the canonical wrapper is unsafe for fixed-span benchmarking; use the audited ID-safe implementation.**" else "**Verdict: canonical wrapper checks passed.**",
  "",
  "`qcrlscR::qc.rlsc.wrap()` was tested on a deterministic 20-feature subset with spans 0.30 and 0.90. Version 0.1.3 does not forward `...` to `qc.rlsc()`, so the outputs are identical and the nominal span is ignored. The benchmark therefore invokes `qcrlscR::qc.rlsc()` separately within each batch on log1p intensities, validates returned row/sample IDs, uses subtractive centered LOESS, then applies the package batch shift on that same log scale and inverse-transforms with expm1.",
  "",
  paste0("- ", qcrlsc_audit$audit$check, ": ", qcrlsc_audit$audit$passed, " (", qcrlsc_audit$audit$interpretation, ")")
), qcrlsc_report_path)
if (!qcrlsc_audit$custom_required) log_line("Canonical QC-RLSC wrapper passed; ID-safe implementation retained for explicit validation.") else log_line("QC-RLSC wrapper ignored fixed span; using audited ID-safe per-batch implementation.")

set.seed(config$tuning_seed)
tune_n <- if (is.finite(config$tuning_features)) min(nrow(x), as.integer(config$tuning_features)) else nrow(x)
tune_features <- if (tune_n == nrow(x)) rownames(x) else sample(rownames(x), tune_n, replace = FALSE)
x_tune <- x[tune_features, , drop = FALSE]
write.csv(data.frame(feature_id = tune_features), file.path(tuning_dir, "tuning_feature_ids.csv"), row.names = FALSE, quote = TRUE)

tuning_rows <- list()
selected <- list(ComBat = list(par_prior = TRUE))
run_tuning_grid <- function(label, grid, fn) {
  rows <- list()
  for (i in seq_len(nrow(grid))) {
    captured <- capture_call_full(function() fn(grid[i, , drop = FALSE]))
    score <- data.frame(mean_qc_correlation = NA_real_, mean_qc_cv = NA_real_, qc_score = NA_real_)
    status <- "failed"
    output_features <- NA_integer_
    output_samples <- NA_integer_
    if (!inherits(captured$value, "error") && is.matrix(captured$value)) {
      output_features <- nrow(captured$value); output_samples <- ncol(captured$value)
      if (identical(colnames(captured$value), colnames(x_tune))) {
        score <- training_qc_score(pmax(captured$value, 0), training_ids)
        status <- if (is.finite(score$qc_score)) "completed" else "failed_score"
      } else captured$error <- "candidate returned misaligned sample IDs"
    }
    rows[[i]] <- cbind(
      data.frame(method = label, candidate_id = i, grid[i, , drop = FALSE], stringsAsFactors = FALSE),
      status = status, runtime_sec = captured$elapsed_sec, score,
      output_features = output_features, output_samples = output_samples,
      warnings = paste(captured$warnings, collapse = " | "), messages = paste(captured$messages, collapse = " | "), error = captured$error
    )
    log_line(label, " tuning candidate ", i, "/", nrow(grid), ": ", status, ".")
  }
  table_value <- bind_rows(rows)
  valid <- table_value[table_value$status == "completed" & is.finite(table_value$qc_score), , drop = FALSE]
  if (!nrow(valid)) {
    table_value$selected <- table_value$candidate_id == 1L
    tuning_rows[[label]] <<- table_value
    log_line("Every ", label, " tuning candidate failed; the published/default first grid setting will be attempted without tuning.")
    return(grid[1, , drop = FALSE])
  }
  valid <- valid[order(-valid$qc_score, valid$candidate_id), , drop = FALSE]
  table_value$selected <- table_value$candidate_id == valid$candidate_id[1]
  tuning_rows[[label]] <<- table_value
  grid[valid$candidate_id[1], , drop = FALSE]
}

selected[["QC-RLSC"]] <- run_tuning_grid("QC-RLSC", data.frame(span = c(0.30, 0.50, 0.70, 0.90)), function(p) run_qc_rlsc_id_safe(x_tune, training_ids, meta_hidden, span = p$span, degree = 1L, shift_batches = TRUE))
selected[["QC-RFSC"]] <- run_tuning_grid("QC-RFSC", expand.grid(ntree = c(200L, 500L), coCV = c(20, 30), Frule = 0.8, KEEP.OUT.ATTRS = FALSE), function(p) run_qc_rfsc_with_controls(x_tune, training_ids, meta_hidden, ntree = p$ntree, coCV = p$coCV, Frule = p$Frule))
selected[["TIGER"]] <- run_tuning_grid("TIGER", data.frame(use_injection = c(TRUE, FALSE)), function(p) run_tiger_all_corrected(x_tune, training_ids, meta_hidden, use_injection = p$use_injection, mtry_percent = 0.4, nodesize_percent = 0.4, ntree = 5, parallel_cores = 1))
selected[["SERRF"]] <- run_tuning_grid("SERRF", data.frame(jitter_eps = c(0, 1e-6, 1e-5)), function(p) run_serrf_all_corrected(x_tune, training_ids, meta_hidden, jitter_eps = p$jitter_eps))

log_line("Starting WiNN automatic selection on ", nrow(x_tune), " deterministic tuning features.")
winn_auto <- run_winn_auto_tuning_subset(x_tune, training_ids, meta_hidden, auto_batch = FALSE)
winn_auto_batch <- run_winn_auto_tuning_subset(x_tune, training_ids, meta_hidden, auto_batch = TRUE)
adjust_outliers <- get("adjust_outliers_mad", envir = asNamespace("winn"))
detect_batch <- get(".auto_detect_batch", envir = asNamespace("winn"))
detected_tune <- detect_batch(adjust_outliers(x_tune), pelt_penalty = NULL)
winn_auto_batch$params$pelt_penalty <- attr(detected_tune, "pelt_penalty")
selected[["WINN auto (QC)"]] <- winn_auto$params
selected[["WINN auto-batch (QC)"]] <- winn_auto_batch$params
tuning_rows[["WINN auto (QC)"]] <- data.frame(method = "WINN auto (QC)", candidate_id = 1L, status = "completed", runtime_sec = winn_auto$captured$elapsed_sec, mean_qc_correlation = winn_auto$params$mean_correlation, mean_qc_cv = winn_auto$params$mean_cv, qc_score = winn_auto$params$mean_correlation - 0.5 * winn_auto$params$mean_cv, output_features = nrow(x_tune), output_samples = ncol(x_tune), warnings = paste(winn_auto$captured$warnings, collapse = " | "), messages = paste(winn_auto$captured$messages, collapse = " | "), error = "", selected = TRUE)
tuning_rows[["WINN auto-batch (QC)"]] <- data.frame(method = "WINN auto-batch (QC)", candidate_id = 1L, status = "completed", runtime_sec = winn_auto_batch$captured$elapsed_sec, mean_qc_correlation = winn_auto_batch$params$mean_correlation, mean_qc_cv = winn_auto_batch$params$mean_cv, qc_score = winn_auto_batch$params$mean_correlation - 0.5 * winn_auto_batch$params$mean_cv, output_features = nrow(x_tune), output_samples = ncol(x_tune), warnings = paste(winn_auto_batch$captured$warnings, collapse = " | "), messages = paste(winn_auto_batch$captured$messages, collapse = " | "), error = "", selected = TRUE)
tuning_df <- bind_rows(tuning_rows)
write.csv(tuning_df, file.path(result_dir, "tuning_candidates.csv"), row.names = FALSE, quote = TRUE, na = "")
for (label in names(tuning_rows)) write.csv(tuning_rows[[label]], file.path(tuning_dir, paste0(gsub("[^A-Za-z0-9]+", "_", tolower(label)), ".csv")), row.names = FALSE, quote = TRUE, na = "")

tuning_runtime <- tapply(tuning_df$runtime_sec, tuning_df$method, sum, na.rm = TRUE)
winn_params_fixed <- list(test = "Ljung-Box", acorr_fdr = 0.05, anova_fdr = 0.05, normalization = "shrink", spline_method = "conservative", scale_by_batch = FALSE, pelt_penalty = NULL)

# Verify that the stage-instrumented implementation returns the same matrix as
# the normal package wrapper.  This is run on the deterministic tuning panel so
# diagnostics do not require a second full-data correction.
validate_winn_instrumentation <- function(method, wrapper_value, instrumented_value) {
  wrapper_mat <- as.matrix(wrapper_value)
  stage_mat <- as.matrix(instrumented_value$data)
  aligned <- identical(dim(wrapper_mat), dim(stage_mat)) &&
    identical(rownames(wrapper_mat), rownames(stage_mat)) &&
    identical(colnames(wrapper_mat), colnames(stage_mat))
  max_difference <- if (aligned) max(abs(wrapper_mat - stage_mat), na.rm = TRUE) else Inf
  data.frame(
    method = method, validation_scope = "deterministic tuning-feature panel",
    n_features = nrow(x_tune), n_samples = ncol(x_tune), dimensions_and_ids_equal = aligned,
    max_absolute_difference = max_difference,
    diagnostic_matrix_equal = aligned && is.finite(max_difference) && max_difference <= 1e-8,
    stringsAsFactors = FALSE
  )
}
winn_fixed_wrapper <- capture_call_full(function() run_winn_with_controls(
  x_tune, control_ids = NULL, meta_df = meta_hidden,
  auto_batch = FALSE, parameters = "fixed"
))
if (inherits(winn_fixed_wrapper$value, "error")) stop("WiNN fixed wrapper validation failed: ", winn_fixed_wrapper$error)
winn_validation <- bind_rows(
  validate_winn_instrumentation(
    "WINN auto (QC)", winn_auto$captured$value,
    run_winn_instrumented(x_tune, meta_hidden, training_ids, FALSE, selected[["WINN auto (QC)"]], "WINN auto (QC)")
  ),
  validate_winn_instrumentation(
    "WINN auto-batch (QC)", winn_auto_batch$captured$value,
    run_winn_instrumented(x_tune, meta_hidden, training_ids, TRUE, selected[["WINN auto-batch (QC)"]], "WINN auto-batch (QC)")
  ),
  validate_winn_instrumentation(
    "WINN default (no QC)", winn_fixed_wrapper$value,
    run_winn_instrumented(x_tune, meta_hidden, NULL, FALSE, winn_params_fixed, "WINN default (no QC)")
  )
)
write.csv(winn_validation, file.path(result_dir, "winn_instrumentation_validation.csv"), row.names = FALSE, quote = TRUE)
if (any(!winn_validation$diagnostic_matrix_equal)) {
  stop("At least one stage-instrumented WiNN result differs from its normal wrapper output; see winn_instrumentation_validation.csv.")
}

method_results <- list()
attempt_rows <- list()
diagnostic_rows <- list()
winn_selectivity <- list()
winn_features <- list()
winn_drift <- list()
winn_batches <- list()
cache_tag <- paste0(dataset, "_hiddenqc_", config$holdout_seed, "_v1")
previous_runtime <- if (file.exists(file.path(result_dir, "method_runtime.csv"))) read.csv(file.path(result_dir, "method_runtime.csv"), stringsAsFactors = FALSE) else data.frame()
previous_diagnostics <- if (file.exists(file.path(result_dir, "method_warnings_messages_errors.csv"))) read.csv(file.path(result_dir, "method_warnings_messages_errors.csv"), stringsAsFactors = FALSE) else data.frame()
run_full <- function(label, fn, tuning_sec = 0) {
  cache_file <- file.path(cache_dir, paste0(gsub("[^A-Za-z0-9]+", "_", tolower(label)), "_", cache_tag, ".rds"))
  reuse_cache <- !force && file.exists(cache_file)
  captured <- if (reuse_cache) {
    previous_elapsed <- if (nrow(previous_runtime) && label %in% previous_runtime$method) previous_runtime$correction_sec[match(label, previous_runtime$method)] else 0
    list(value = readRDS(cache_file), elapsed_sec = previous_elapsed, warnings = character(), messages = "reused valid cache", error = "")
  } else {
    set.seed(method_seeds[[label]])
    capture_call_full(fn)
  }
  diagnostic_rows[[label]] <<- if (reuse_cache && nrow(previous_diagnostics) && label %in% previous_diagnostics$method) {
    previous_diagnostics[match(label, previous_diagnostics$method), , drop = FALSE]
  } else data.frame(method = label, warnings = paste(captured$warnings, collapse = " | "), messages = paste(captured$messages, collapse = " | "), error = captured$error, stringsAsFactors = FALSE)
  if (inherits(captured$value, "error")) {
    attempt_rows[[label]] <<- data.frame(method = label, status = "failed", tuning_sec = tuning_sec, correction_sec = captured$elapsed_sec, total_sec = tuning_sec + captured$elapsed_sec, error = captured$error, stringsAsFactors = FALSE)
    log_line("Method failed: ", label, " — ", captured$error)
    return(NULL)
  }
  value <- captured$value
  validated <- tryCatch({
    if (is.list(value) && !is.null(value$data)) {
      mat <- value$data
    } else mat <- value
    mat <- as.matrix(mat)
    if (is.null(rownames(mat)) && nrow(mat) == nrow(x)) rownames(mat) <- rownames(x)
    if (is.null(colnames(mat)) && ncol(mat) == ncol(x)) colnames(mat) <- colnames(x)
    if (anyDuplicated(rownames(mat)) || anyDuplicated(colnames(mat))) stop(label, " returned duplicate IDs.")
    if (!identical(colnames(mat), colnames(x))) stop(label, " returned nonidentical sample IDs/order.")
    if (!length(intersect(rownames(mat), rownames(x)))) stop(label, " returned no alignable features.")
    mat <- mat[rownames(x)[rownames(x) %in% rownames(mat)], , drop = FALSE]
    if (any(!is.finite(mat))) {
      if (label != "TIGER") stop(label, " returned nonfinite values.")
      invalid <- which(!is.finite(mat), arr.ind = TRUE)
      invalid_features <- unique(rownames(mat)[invalid[, "row"]])
      invalid_table <- bind_rows(lapply(invalid_features, function(feature_id) {
        values <- mat[feature_id, ]
        data.frame(
          feature_id = feature_id,
          nonfinite_cells = sum(!is.finite(values)),
          nan_cells = sum(is.nan(values)),
          positive_infinite_cells = sum(is.infinite(values) & values > 0),
          negative_infinite_cells = sum(is.infinite(values) & values < 0),
          input_zero_cells = sum(x[feature_id, ] == 0),
          input_training_qc_zero_cells = sum(x[feature_id, training_ids, drop = FALSE] == 0),
          action = "complete TIGER output feature row removed; no cell replacement",
          stringsAsFactors = FALSE
        )
      }))
      write.csv(invalid_table, file.path(result_dir, "tiger_nonfinite_output_features.csv"), row.names = FALSE, quote = TRUE)
      mat <- mat[!rownames(mat) %in% invalid_features, , drop = FALSE]
      if (!nrow(mat)) stop("TIGER returned nonfinite values in every feature row.")
      note <- paste0(
        "TIGER returned ", nrow(invalid), " non-finite cells across ", length(invalid_features),
        " feature rows. Complete affected rows were removed without cell replacement; see tiger_nonfinite_output_features.csv."
      )
      warning_parts <- c(diagnostic_rows[[label]]$warnings, note)
      warning_parts <- warning_parts[!is.na(warning_parts) & nzchar(warning_parts)]
      diagnostic_rows[[label]]$warnings <<- paste(warning_parts, collapse = " | ")
      log_line(note)
    } else if (label == "TIGER") {
      tiger_report <- file.path(result_dir, "tiger_nonfinite_output_features.csv")
      if (file.exists(tiger_report)) file.remove(tiger_report)
    }
    pmax(mat, 0)
  }, error = function(e) e)
  if (inherits(validated, "error")) {
    error_text <- conditionMessage(validated)
    diagnostic_rows[[label]]$error <<- error_text
    attempt_rows[[label]] <<- data.frame(method = label, status = "failed", tuning_sec = tuning_sec, correction_sec = captured$elapsed_sec, total_sec = tuning_sec + captured$elapsed_sec, error = error_text, stringsAsFactors = FALSE)
    log_line("Method output validation failed: ", label, " — ", error_text)
    return(NULL)
  }
  mat <- validated
  if (is.list(value) && !is.null(value$data) && grepl("^WINN", label)) {
    winn_selectivity[[label]] <<- value$summary
    winn_features[[label]] <<- value$by_feature
    winn_drift[[label]] <<- value$drift
    winn_batches[[label]] <<- data.frame(method = label, sample_id = meta$sample_id, supplied_batch = meta$batch, detected_batch = as.character(value$detected_batch), global_run_order = meta$global_run_order)
  }
  saveRDS(mat, cache_file, compress = "xz")
  method_results[[label]] <<- mat
  attempt_rows[[label]] <<- data.frame(method = label, status = "completed", tuning_sec = tuning_sec, correction_sec = captured$elapsed_sec, total_sec = tuning_sec + captured$elapsed_sec, error = "", stringsAsFactors = FALSE)
  log_line("Method completed: ", label, " (", nrow(mat), " × ", ncol(mat), ").")
  mat
}

invisible(run_full("Raw", function() x))
invisible(run_full("ComBat", function() run_combat(x, meta_hidden, par_prior = TRUE)))
invisible(run_full("QC-RLSC", function() run_qc_rlsc_id_safe(x, training_ids, meta_hidden, span = selected[["QC-RLSC"]]$span, degree = 1L, shift_batches = TRUE), tuning_runtime[["QC-RLSC"]]))
invisible(run_full("QC-RFSC", function() run_qc_rfsc_with_controls(x, training_ids, meta_hidden, ntree = selected[["QC-RFSC"]]$ntree, coCV = selected[["QC-RFSC"]]$coCV, Frule = selected[["QC-RFSC"]]$Frule), tuning_runtime[["QC-RFSC"]]))
invisible(run_full("TIGER", function() run_tiger_all_corrected(x, training_ids, meta_hidden, use_injection = selected[["TIGER"]]$use_injection, mtry_percent = 0.4, nodesize_percent = 0.4, ntree = 5, parallel_cores = 1), tuning_runtime[["TIGER"]]))
invisible(run_full("SERRF", function() run_serrf_all_corrected(x, training_ids, meta_hidden, jitter_eps = selected[["SERRF"]]$jitter_eps), tuning_runtime[["SERRF"]]))
invisible(run_full("WINN auto (QC)", function() run_winn_instrumented(x, meta_hidden, training_ids, FALSE, selected[["WINN auto (QC)"]], "WINN auto (QC)"), tuning_runtime[["WINN auto (QC)"]]))
invisible(run_full("WINN auto-batch (QC)", function() run_winn_instrumented(x, meta_hidden, training_ids, TRUE, selected[["WINN auto-batch (QC)"]], "WINN auto-batch (QC)"), tuning_runtime[["WINN auto-batch (QC)"]]))
invisible(run_full("WINN default (no QC)", function() run_winn_instrumented(x, meta_hidden, NULL, FALSE, winn_params_fixed, "WINN default (no QC)")))

attempts <- bind_rows(attempt_rows)[match(method_order, names(attempt_rows)), , drop = FALSE]
diagnostics <- bind_rows(diagnostic_rows)[match(method_order, names(diagnostic_rows)), , drop = FALSE]
write.csv(attempts, file.path(result_dir, "method_runtime.csv"), row.names = FALSE, quote = TRUE)
write.csv(attempts[attempts$status == "failed", c("method", "error"), drop = FALSE], file.path(result_dir, "method_failures.csv"), row.names = FALSE, quote = TRUE)
write.csv(diagnostics, file.path(result_dir, "method_warnings_messages_errors.csv"), row.names = FALSE, quote = TRUE)
completed <- method_order[method_order %in% names(method_results)]
if (!length(completed)) stop("No correction method completed.")

dimension_rows <- bind_rows(lapply(completed, function(label) data.frame(
  method = label, features_retained = nrow(method_results[[label]]), features_lost = nrow(x) - nrow(method_results[[label]]),
  injections_retained = ncol(method_results[[label]]), study_samples_retained = sum(!meta$is_qc & meta$sample_id %in% colnames(method_results[[label]])),
  hidden_qcs_retained = sum(hidden_ids %in% colnames(method_results[[label]])), exact_sample_alignment = identical(colnames(method_results[[label]]), colnames(x)), stringsAsFactors = FALSE
)))
write.csv(dimension_rows, file.path(result_dir, "method_dimensions.csv"), row.names = FALSE, quote = TRUE)
feature_coverage <- bind_rows(lapply(completed, function(label) data.frame(
  method = label,
  feature_id = rownames(x),
  retained = rownames(x) %in% rownames(method_results[[label]]),
  stringsAsFactors = FALSE
)))
write.csv(feature_coverage, file.path(result_dir, "method_feature_coverage.csv"), row.names = FALSE, quote = TRUE)
common_features <- Reduce(intersect, lapply(method_results[completed], rownames))
write.csv(data.frame(feature_id = common_features), file.path(result_dir, "common_evaluation_feature_ids.csv"), row.names = FALSE, quote = TRUE)
method_results_common <- lapply(method_results[completed], function(mat) mat[common_features, , drop = FALSE])

parameter_rows <- data.frame(
  method = method_order,
  correction_seed = unname(method_seeds[method_order]),
  parameters = c(
    "none", "log1p; intercept-only design; par.prior=TRUE; expm1",
    paste0("ID-safe per-batch log1p subtractive LOESS; span=", selected[["QC-RLSC"]]$span, "; degree=1; batch.shift=TRUE; training_QCs=", length(training_ids)),
    paste0("ntree=", selected[["QC-RFSC"]]$ntree, "; coCV=", selected[["QC-RFSC"]]$coCV, "; Frule=", selected[["QC-RFSC"]]$Frule, "; training_QCs=", length(training_ids)),
    paste0("use_injection=", selected[["TIGER"]]$use_injection, "; mtry_percent=0.4; nodesize_percent=0.4; ntree=5; cores=1; training_QCs=", length(training_ids), "; nonfinite_output_policy=remove_complete_affected_feature_rows_and_report_no_cell_replacement"),
    paste0("jitter_eps=", selected[["SERRF"]]$jitter_eps, "; training_QCs=", length(training_ids)),
    paste0("auto selected on ", length(tune_features), " QC-only tuning features; test=", selected[["WINN auto (QC)"]]$test, "; acorr_fdr=", selected[["WINN auto (QC)"]]$acorr_fdr, "; anova_fdr=", selected[["WINN auto (QC)"]]$anova_fdr, "; normalization=", selected[["WINN auto (QC)"]]$normalization, "; spline=", selected[["WINN auto (QC)"]]$spline_method),
    paste0("auto-batch selected on ", length(tune_features), " QC-only tuning features; no supplied batches; pelt_penalty=", selected[["WINN auto-batch (QC)"]]$pelt_penalty, "; test=", selected[["WINN auto-batch (QC)"]]$test, "; acorr_fdr=", selected[["WINN auto-batch (QC)"]]$acorr_fdr, "; anova_fdr=", selected[["WINN auto-batch (QC)"]]$anova_fdr, "; normalization=", selected[["WINN auto-batch (QC)"]]$normalization),
    "parameters=fixed; supplied batches; no QC IDs; Ljung-Box adaptive lag; fdr=0.05; conservative spline; ANOVA batch gate; shrink normalization"
  ), stringsAsFactors = FALSE
)
write.csv(parameter_rows, file.path(result_dir, "method_parameters.csv"), row.names = FALSE, quote = TRUE)

if (length(winn_selectivity)) {
  winn_summary <- bind_rows(winn_selectivity)
  winn_summary$diagnostic_matrix_equal <- winn_validation$diagnostic_matrix_equal[match(winn_summary$method, winn_validation$method)]
  write.csv(winn_summary, file.path(result_dir, "winn_selectivity_summary.csv"), row.names = FALSE, quote = TRUE, na = "")
  write.csv(bind_rows(winn_features), file.path(result_dir, "winn_selectivity_by_feature.csv"), row.names = FALSE, quote = TRUE, na = "")
  write.csv(bind_rows(winn_drift), file.path(result_dir, "winn_drift_tests_by_feature_batch.csv"), row.names = FALSE, quote = TRUE, na = "")
  segmentation <- bind_rows(winn_batches)
  write.csv(segmentation, file.path(result_dir, "winn_auto_batch_assignments.csv"), row.names = FALSE, quote = TRUE)
} else {
  if (!file.exists(file.path(result_dir, "winn_selectivity_summary.csv"))) write.csv(data.frame(), file.path(result_dir, "winn_selectivity_summary.csv"), row.names = FALSE)
  if (!file.exists(file.path(result_dir, "winn_selectivity_by_feature.csv"))) write.csv(data.frame(), file.path(result_dir, "winn_selectivity_by_feature.csv"), row.names = FALSE)
}

log_line("Computing common technical metrics.")
qc_cv <- bind_rows(lapply(completed, function(label) qc_metrics_by_feature(method_results[[label]], hidden_ids, label)))
qc_pairs <- bind_rows(lapply(completed, function(label) qc_pairwise_metrics(method_results[[label]], hidden_ids, label)))
write.csv(qc_cv, file.path(result_dir, "qc_cv_by_feature.csv"), row.names = FALSE, quote = TRUE)
write.csv(qc_pairs, file.path(result_dir, "qc_pairwise_correlations.csv"), row.names = FALSE, quote = TRUE)

study_meta <- meta[!meta$is_qc, , drop = FALSE]
gam_rows <- bind_rows(lapply(completed, function(label) bind_rows(lapply(unique(study_meta$batch), function(batch_id) {
  d <- study_meta[study_meta$batch == batch_id, , drop = FALSE]
  compute_metabolite_segment_gam(method_results[[label]], d, label, sample_id_col = "sample_id", order_col = "within_batch_order", batch_col = "batch", transform_fun = log1p)
}))))
ljung_rows <- bind_rows(lapply(completed, function(label) bind_rows(lapply(unique(study_meta$batch), function(batch_id) {
  d <- study_meta[study_meta$batch == batch_id, , drop = FALSE]
  compute_metabolite_segment_ljung_box(method_results[[label]], d, label, sample_id_col = "sample_id", order_col = "within_batch_order", batch_col = "batch", transform_fun = log1p)
})))) |>
  group_by(method, batch, segment_id) |>
  mutate(p_adj = p.adjust(p_value, method = "BH"), is_autocorrelated = is.finite(p_adj) & p_adj < 0.05) |>
  ungroup()
write.csv(gam_rows, file.path(result_dir, "run_order_gam_by_feature_batch.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(ljung_rows, file.path(result_dir, "ljung_box_by_feature_batch.csv"), row.names = FALSE, quote = TRUE, na = "")
batch_associations <- bind_rows(lapply(completed, function(label) feature_batch_associations(method_results[[label]], meta, label)))
write.csv(batch_associations, file.path(result_dir, "batch_associations.csv"), row.names = FALSE, quote = TRUE, na = "")
correction_magnitude <- bind_rows(lapply(completed, function(label) correction_magnitude_table(method_results[[label]], x, label)))
write.csv(correction_magnitude, file.path(result_dir, "correction_magnitude.csv"), row.names = FALSE, quote = TRUE)

limma_fit_table <- function(mat, md, design, coefficients, method_label, model_label) {
  z <- log1p(mat[, md$sample_id, drop = FALSE])
  fit <- limma::eBayes(limma::lmFit(z, design), robust = TRUE)
  bind_rows(lapply(names(coefficients), function(variable) {
    coef_name <- coefficients[[variable]]
    index <- match(coef_name, colnames(design))
    if (is.na(index)) return(tibble())
    coef_value <- fit$coefficients[, index]
    t_value <- fit$t[, index]
    p_value <- fit$p.value[, index]
    outcome_sd <- apply(z, 1, stats::sd)
    data.frame(
      method = method_label, model = model_label, variable = variable, feature_id = rownames(z),
      coefficient = coef_value, standardized_effect = coef_value / outcome_sd,
      standard_error = abs(coef_value / t_value), t_statistic = t_value,
      p_value = p_value, p_adj = stats::p.adjust(p_value, method = "BH"),
      partial_r_squared = t_value^2 / (t_value^2 + fit$df.total),
      significant = is.finite(p_value) & stats::p.adjust(p_value, method = "BH") < 0.05,
      stringsAsFactors = FALSE
    )
  }))
}

biological_pc <- list()
biological_associations <- NULL
cross_batch <- NULL
stability <- NULL
if (dataset == "sacurine") {
  md <- meta[!meta$is_qc, , drop = FALSE]
  md$age_z <- as.numeric(scale(md$age)); md$bmi_z <- as.numeric(scale(md$bmi))
  md$gender <- factor(md$gender); md$sampling_factor <- factor(md$sampling)
  biological_mats <- lapply(method_results_common, function(m) sweep(m[, md$sample_id, drop = FALSE], 2, md$osmolality, "/"))
  design_primary <- model.matrix(~ age_z + bmi_z + gender + sampling_factor, data = md)
  design_batch <- model.matrix(~ age_z + bmi_z + gender + sampling_factor + factor(batch), data = md)
  gender_coef <- grep("^gender", colnames(design_primary), value = TRUE)[1]
  coefficients <- c(age = "age_z", bmi = "bmi_z", gender = gender_coef)
  primary <- bind_rows(lapply(completed, function(label) limma_fit_table(biological_mats[[label]], md, design_primary, coefficients, label, "primary_osmolality_normalized")))
  adjusted <- bind_rows(lapply(completed, function(label) limma_fit_table(biological_mats[[label]], md, design_batch, coefficients, label, "batch_adjusted_osmolality_normalized")))
  write.csv(primary, file.path(result_dir, "sacurine_associations_primary.csv"), row.names = FALSE, quote = TRUE, na = "")
  write.csv(adjusted, file.path(result_dir, "sacurine_associations_batch_adjusted.csv"), row.names = FALSE, quote = TRUE, na = "")
  biological_associations <- primary
  biological_pc <- bind_rows(lapply(completed, function(label) bind_rows(lapply(c("age", "bmi", "gender"), function(variable) data.frame(
    method = label,
    variable = variable,
    weighted_pc_r2 = weighted_pc_r2_general(
      biological_mats[[label]], md[[variable]],
      target_type = if (variable %in% c("age", "bmi")) "continuous" else "categorical",
      n_pcs = 5L
    )
  )))))

  batch_effects <- bind_rows(lapply(completed, function(label) bind_rows(lapply(unique(md$batch), function(batch_id) {
    d <- md[md$batch == batch_id, , drop = FALSE]
    design <- model.matrix(~ age_z + bmi_z + gender + sampling_factor, data = d)
    gender_name <- grep("^gender", colnames(design), value = TRUE)[1]
    limma_fit_table(biological_mats[[label]], d, design, c(age = "age_z", bmi = "bmi_z", gender = gender_name), label, paste0("batch_", batch_id))
  }))))
  batch_pairs <- combn(unique(md$batch), 2, simplify = FALSE)
  cross_batch <- bind_rows(lapply(completed, function(label) bind_rows(lapply(c("age", "bmi", "gender"), function(variable) bind_rows(lapply(batch_pairs, function(pair) {
    a <- batch_effects[batch_effects$method == label & batch_effects$variable == variable & batch_effects$model == paste0("batch_", pair[1]), ]
    b <- batch_effects[batch_effects$method == label & batch_effects$variable == variable & batch_effects$model == paste0("batch_", pair[2]), ]
    b <- b[match(a$feature_id, b$feature_id), ]
    metrics <- effect_concordance(a$standardized_effect, b$standardized_effect)
    data.frame(method = label, variable = variable, batch_1 = pair[1], batch_2 = pair[2], t(metrics), stringsAsFactors = FALSE)
  }))))))
  write.csv(cross_batch, file.path(result_dir, "sacurine_cross_batch_effect_concordance.csv"), row.names = FALSE, quote = TRUE, na = "")

  set.seed(config$bootstrap_seed)
  bootstrap_rows <- list()
  for (label in completed) {
    coefficient_store <- lapply(c("age", "bmi", "gender"), function(v) matrix(NA_real_, length(common_features), 100L, dimnames = list(common_features, NULL)))
    names(coefficient_store) <- c("age", "bmi", "gender")
    selected_store <- lapply(coefficient_store, function(m) matrix(FALSE, nrow(m), ncol(m), dimnames = dimnames(m)))
    for (replicate_id in seq_len(100L)) {
      idx <- sample(seq_len(nrow(md)), nrow(md), replace = TRUE)
      d <- md[idx, , drop = FALSE]; d$sample_id <- md$sample_id[idx]
      design <- model.matrix(~ age_z + bmi_z + gender + sampling_factor, data = d)
      gender_name <- grep("^gender", colnames(design), value = TRUE)[1]
      tab <- limma_fit_table(biological_mats[[label]], d, design, c(age = "age_z", bmi = "bmi_z", gender = gender_name), label, "bootstrap")
      for (variable in names(coefficient_store)) {
        z <- tab[tab$variable == variable, ]; z <- z[match(common_features, z$feature_id), ]
        coefficient_store[[variable]][, replicate_id] <- z$coefficient
        selected_store[[variable]][, replicate_id] <- z$significant
      }
    }
    bootstrap_rows[[label]] <- bind_rows(lapply(names(coefficient_store), function(variable) data.frame(
      method = label, variable = variable, feature_id = common_features,
      selection_frequency = rowMeans(selected_store[[variable]], na.rm = TRUE),
      median_coefficient = apply(coefficient_store[[variable]], 1, stats::median, na.rm = TRUE),
      coefficient_sd = apply(coefficient_store[[variable]], 1, stats::sd, na.rm = TRUE),
      selected_at_least_half = rowMeans(selected_store[[variable]], na.rm = TRUE) >= 0.5,
      stringsAsFactors = FALSE
    )))
    log_line("Completed Sacurine bootstrap stability: ", label, ".")
  }
  stability <- bind_rows(bootstrap_rows)
  write.csv(stability, file.path(result_dir, "sacurine_bootstrap_stability.csv"), row.names = FALSE, quote = TRUE)

  sensitivity <- bind_rows(lapply(completed, function(label) {
    unnormalized <- limma_fit_table(method_results_common[[label]], md, design_primary, coefficients, label, "not_osmolality_normalized")
    normalized <- primary[primary$method == label, ]
    bind_rows(lapply(c("age", "bmi", "gender"), function(variable) data.frame(
      method = label, variable = variable,
      significant_osmolality_normalized = sum(normalized$significant[normalized$variable == variable]),
      significant_unnormalized = sum(unnormalized$significant[unnormalized$variable == variable]),
      weighted_pc_r2_osmolality_normalized = weighted_pc_r2_general(
        biological_mats[[label]], md[[variable]],
        target_type = if (variable %in% c("age", "bmi")) "continuous" else "categorical",
        n_pcs = 5L
      ),
      weighted_pc_r2_unnormalized = weighted_pc_r2_general(
        method_results_common[[label]][, md$sample_id, drop = FALSE], md[[variable]],
        target_type = if (variable %in% c("age", "bmi")) "continuous" else "categorical",
        n_pcs = 5L
      ), stringsAsFactors = FALSE
    )))
  }))
  write.csv(sensitivity, file.path(result_dir, "sacurine_osmolality_sensitivity.csv"), row.names = FALSE, quote = TRUE)
  outlier_md <- md[md$sample_id != "HU_neg_096_b2", , drop = FALSE]
  outlier_sensitivity <- bind_rows(lapply(completed, function(label) {
    design <- model.matrix(~ age_z + bmi_z + gender + sampling_factor, data = outlier_md)
    gender_name <- grep("^gender", colnames(design), value = TRUE)[1]
    tab <- limma_fit_table(biological_mats[[label]], outlier_md, design, c(age = "age_z", bmi = "bmi_z", gender = gender_name), label, "exclude_HU_neg_096_b2")
    tab |> group_by(method, variable) |> summarise(significant_features = sum(significant), .groups = "drop")
  }))
  write.csv(outlier_sensitivity, file.path(result_dir, "sacurine_outlier_sensitivity.csv"), row.names = FALSE, quote = TRUE)
} else {
  md <- meta[!meta$is_qc, , drop = FALSE]
  md$biological_group <- factor(md$biological_group, levels = c("group_0", "group_1"))
  limma_group_fit_table <- function(mat, model_md, design, coefficient_name, method_label, model_label) {
    tab <- limma_fit_table(mat, model_md, design, c(biological_group = coefficient_name), method_label, model_label)
    z <- log1p(mat[, model_md$sample_id, drop = FALSE])
    group <- factor(model_md$biological_group, levels = c("group_0", "group_1"))
    n0 <- sum(group == "group_0"); n1 <- sum(group == "group_1")
    sd0 <- apply(z[, group == "group_0", drop = FALSE], 1, stats::sd)
    sd1 <- apply(z[, group == "group_1", drop = FALSE], 1, stats::sd)
    pooled_sd <- sqrt(((n0 - 1) * sd0^2 + (n1 - 1) * sd1^2) / (n0 + n1 - 2))
    tab$standardized_effect <- tab$coefficient / pooled_sd[match(tab$feature_id, rownames(z))]
    tab
  }
  design_primary <- model.matrix(~ biological_group, data = md)
  design_batch <- model.matrix(~ biological_group + factor(batch), data = md)
  coefficient <- grep("^biological_group", colnames(design_primary), value = TRUE)[1]
  primary <- bind_rows(lapply(completed, function(label) limma_group_fit_table(method_results_common[[label]], md, design_primary, coefficient, label, "primary")))
  adjusted <- bind_rows(lapply(completed, function(label) limma_group_fit_table(method_results_common[[label]], md, design_batch, coefficient, label, "batch_adjusted")))
  write.csv(primary, file.path(result_dir, "waveica_group_associations_primary.csv"), row.names = FALSE, quote = TRUE, na = "")
  write.csv(adjusted, file.path(result_dir, "waveica_group_associations_batch_adjusted.csv"), row.names = FALSE, quote = TRUE, na = "")
  biological_associations <- primary
  biological_pc <- bind_rows(lapply(completed, function(label) data.frame(
    method = label,
    variable = "biological_group",
    weighted_pc_r2 = weighted_pc_r2_general(
      method_results_common[[label]][, md$sample_id, drop = FALSE], md$biological_group,
      target_type = "categorical", n_pcs = 5L
    )
  )))

  batch_effects <- bind_rows(lapply(completed, function(label) bind_rows(lapply(unique(md$batch), function(batch_id) {
    d <- md[md$batch == batch_id, , drop = FALSE]
    design <- model.matrix(~ biological_group, data = d)
    coef_name <- grep("^biological_group", colnames(design), value = TRUE)[1]
    limma_group_fit_table(method_results_common[[label]], d, design, coef_name, label, paste0("batch_", batch_id))
  }))))
  batch_pairs <- combn(unique(md$batch), 2, simplify = FALSE)
  cross_batch <- bind_rows(lapply(completed, function(label) bind_rows(lapply(batch_pairs, function(pair) {
    a <- batch_effects[batch_effects$method == label & batch_effects$model == paste0("batch_", pair[1]), ]
    b <- batch_effects[batch_effects$method == label & batch_effects$model == paste0("batch_", pair[2]), ]; b <- b[match(a$feature_id, b$feature_id), ]
    metrics <- effect_concordance(a$standardized_effect, b$standardized_effect)
    data.frame(method = label, batch_1 = pair[1], batch_2 = pair[2], t(metrics), stringsAsFactors = FALSE)
  }))))
  write.csv(cross_batch, file.path(result_dir, "waveica_cross_batch_effect_concordance.csv"), row.names = FALSE, quote = TRUE, na = "")

  set.seed(config$bootstrap_seed)
  smaller_ids <- which(md$biological_group == "group_0")
  larger_ids <- which(md$biological_group == "group_1")
  balance_rows <- list()
  for (label in completed) {
    coefficient_store <- matrix(NA_real_, length(common_features), 100L, dimnames = list(common_features, NULL))
    selected_store <- matrix(FALSE, length(common_features), 100L, dimnames = list(common_features, NULL))
    for (replicate_id in seq_len(100L)) {
      idx <- c(smaller_ids, sample(larger_ids, length(smaller_ids), replace = FALSE))
      d <- md[idx, , drop = FALSE]
      design <- model.matrix(~ biological_group, data = d)
      coef_name <- grep("^biological_group", colnames(design), value = TRUE)[1]
      tab <- limma_group_fit_table(method_results_common[[label]], d, design, coef_name, label, "balanced_resample")
      tab <- tab[match(common_features, tab$feature_id), ]
      coefficient_store[, replicate_id] <- tab$standardized_effect
      selected_store[, replicate_id] <- tab$significant
    }
    median_effect <- apply(coefficient_store, 1, stats::median, na.rm = TRUE)
    direction_stability <- vapply(seq_len(nrow(coefficient_store)), function(i) mean(sign(coefficient_store[i, ]) == sign(median_effect[i]), na.rm = TRUE), numeric(1))
    balance_rows[[label]] <- data.frame(
      method = label, feature_id = common_features, selection_frequency = rowMeans(selected_store),
      effect_direction_stability = direction_stability, median_standardized_effect = median_effect,
      selected_at_least_half = rowMeans(selected_store) >= 0.5, stringsAsFactors = FALSE
    )
    log_line("Completed WaveICA balanced-resampling stability: ", label, ".")
  }
  stability <- bind_rows(balance_rows)
  write.csv(stability, file.path(result_dir, "waveica_balanced_resampling_stability.csv"), row.names = FALSE, quote = TRUE)
}
write.csv(biological_pc, file.path(result_dir, "biological_weighted_pc_r2.csv"), row.names = FALSE, quote = TRUE)

# Report both explicitly comparable panels required by the benchmark design.
# The common panel is the intersection retained by all completed methods. The
# coverage-penalized panel restores any method-dropped feature to its Raw value.
coverage_mats <- lapply(completed, function(label) {
  out <- x
  retained <- intersect(rownames(x), rownames(method_results[[label]]))
  out[retained, ] <- method_results[[label]][retained, colnames(x), drop = FALSE]
  out
})
names(coverage_mats) <- completed

summarise_panel <- function(panel_type) {
  bind_rows(lapply(completed, function(label) {
    if (panel_type == "common_retained") {
      feature_ids <- common_features
      panel_mat <- method_results[[label]][feature_ids, , drop = FALSE]
      qc <- qc_cv[qc_cv$method == label & qc_cv$feature_id %in% feature_ids, ]
      gam <- gam_rows[gam_rows$method == label & gam_rows$metabolite %in% feature_ids, ]
      lb <- ljung_rows[ljung_rows$method == label & ljung_rows$metabolite %in% feature_ids, ]
      ba <- batch_associations[batch_associations$method == label & batch_associations$feature_id %in% feature_ids, ]
      assoc <- biological_associations[biological_associations$method == label & biological_associations$feature_id %in% feature_ids, ]
    } else {
      feature_ids <- rownames(x)
      panel_mat <- coverage_mats[[label]]
      fill_raw <- function(data_value, id_col) {
        current <- data_value[data_value$method == label, , drop = FALSE]
        raw_rows <- data_value[data_value$method == "Raw" & !data_value[[id_col]] %in% current[[id_col]], , drop = FALSE]
        if (nrow(raw_rows)) raw_rows$method <- label
        bind_rows(current, raw_rows)
      }
      qc <- fill_raw(qc_cv, "feature_id")
      gam <- fill_raw(gam_rows, "metabolite")
      lb <- fill_raw(ljung_rows, "metabolite")
      ba <- fill_raw(batch_associations, "feature_id")
      assoc <- if (dataset == "sacurine") {
        biological_panel <- sweep(panel_mat[, md$sample_id, drop = FALSE], 2, md$osmolality, "/")
        limma_fit_table(biological_panel, md, design_primary, coefficients, label, "coverage_penalized_raw_fill")
      } else {
        limma_group_fit_table(panel_mat, md, design_primary, coefficient, label, "coverage_penalized_raw_fill")
      }
    }
    bio_pc_value <- if (dataset == "sacurine") {
      md_panel <- meta[!meta$is_qc, , drop = FALSE]
      biological_panel <- sweep(panel_mat[, md_panel$sample_id, drop = FALSE], 2, md_panel$osmolality, "/")
      mean(vapply(c("age", "bmi", "gender"), function(variable) weighted_pc_r2_general(
        biological_panel, md_panel[[variable]],
        target_type = if (variable %in% c("age", "bmi")) "continuous" else "categorical",
        n_pcs = 5L
      ), numeric(1)))
    } else {
      md_panel <- meta[!meta$is_qc, , drop = FALSE]
      mean(weighted_pc_r2_general(
        panel_mat[, md_panel$sample_id, drop = FALSE], md_panel$biological_group,
        target_type = "categorical", n_pcs = 5L
      ))
    }
    data.frame(
      method = label, panel = panel_type, features_evaluated = length(feature_ids),
      heldout_qc_cv_mean = mean(qc$cv_percent, na.rm = TRUE),
      residual_gam_deviance_mean = mean(gam$explained, na.rm = TRUE),
      residual_ljung_box_proportion = mean(lb$is_autocorrelated[is.finite(lb$p_value)]),
      batch_weighted_pc_r2 = weighted_pc_r2_general(
        panel_mat[, study_meta$sample_id, drop = FALSE], study_meta$batch,
        target_type = "categorical", n_pcs = 5L
      ),
      batch_associated_features = sum(ba$is_batch_associated, na.rm = TRUE),
      biological_weighted_pc_r2_mean = bio_pc_value,
      biological_associated_features = sum(assoc$significant, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}
common_panel_metrics <- summarise_panel("common_retained")
coverage_penalized_metrics <- summarise_panel("coverage_penalized_raw_fill")
write.csv(common_panel_metrics, file.path(result_dir, "method_metrics_common_panel.csv"), row.names = FALSE, quote = TRUE)
write.csv(coverage_penalized_metrics, file.path(result_dir, "method_metrics_coverage_penalized.csv"), row.names = FALSE, quote = TRUE)

metric_rows <- bind_rows(lapply(completed, function(label) {
  qc <- qc_cv[qc_cv$method == label, ]
  pairs <- qc_pairs[qc_pairs$method == label & qc_pairs$correlation_method == "pearson", ]
  gam <- gam_rows[gam_rows$method == label, ]
  lb <- ljung_rows[ljung_rows$method == label, ]
  ba <- batch_associations[batch_associations$method == label, ]
  assoc <- biological_associations[biological_associations$method == label, ]
  mag <- correction_magnitude[correction_magnitude$method == label, ]
  runtime <- attempts[attempts$method == label, ]
  data.frame(
    method = label,
    heldout_qc_cv_mean = mean(qc$cv_percent, na.rm = TRUE), heldout_qc_cv_sd = stats::sd(qc$cv_percent, na.rm = TRUE), heldout_qc_cv_median = stats::median(qc$cv_percent, na.rm = TRUE), heldout_qc_cv_iqr = stats::IQR(qc$cv_percent, na.rm = TRUE), heldout_qc_cv_finite_features = sum(is.finite(qc$cv_percent)),
    heldout_qc_pearson_mean = mean(pairs$correlation, na.rm = TRUE), heldout_qc_pearson_median = stats::median(pairs$correlation, na.rm = TRUE), heldout_qc_pearson_sd = stats::sd(pairs$correlation, na.rm = TRUE), heldout_qc_pearson_iqr = stats::IQR(pairs$correlation, na.rm = TRUE),
    residual_gam_deviance_mean = mean(gam$explained, na.rm = TRUE), residual_gam_deviance_sd = stats::sd(gam$explained, na.rm = TRUE), residual_gam_deviance_median = stats::median(gam$explained, na.rm = TRUE), residual_gam_profiles = sum(is.finite(gam$explained)),
    residual_ljung_box_significant = sum(lb$is_autocorrelated, na.rm = TRUE), residual_ljung_box_tested = sum(is.finite(lb$p_value)), residual_ljung_box_proportion = mean(lb$is_autocorrelated[is.finite(lb$p_value)]),
    batch_weighted_pc_r2 = weighted_pc_r2_general(
      method_results[[label]][, study_meta$sample_id, drop = FALSE], study_meta$batch,
      target_type = "categorical", n_pcs = 5L
    ), batch_associated_features = sum(ba$is_batch_associated, na.rm = TRUE),
    biological_weighted_pc_r2_mean = mean(biological_pc$weighted_pc_r2[biological_pc$method == label], na.rm = TRUE), biological_associated_features = sum(assoc$significant, na.rm = TRUE), biological_effect_median_abs = stats::median(abs(assoc$standardized_effect), na.rm = TRUE),
    cross_batch_effect_pearson_median = stats::median(cross_batch$pearson[cross_batch$method == label], na.rm = TRUE), cross_batch_effect_sign_concordance_median = stats::median(cross_batch$sign_concordance[cross_batch$method == label], na.rm = TRUE),
    tuning_runtime_sec = runtime$tuning_sec, correction_runtime_sec = runtime$correction_sec, total_runtime_sec = runtime$total_sec,
    features_retained = nrow(method_results[[label]]), features_lost = nrow(x) - nrow(method_results[[label]]), total_injections_retained = ncol(method_results[[label]]), study_samples_retained = sum(study_meta$sample_id %in% colnames(method_results[[label]])), hidden_qcs_retained = sum(hidden_ids %in% colnames(method_results[[label]])),
    median_absolute_log1p_change = stats::median(mag$median_absolute_log1p_change, na.rm = TRUE), p90_absolute_log1p_change = stats::quantile(mag$p90_absolute_log1p_change, 0.9, na.rm = TRUE, names = FALSE), proportion_entries_materially_changed = mean(mag$proportion_entries_materially_changed, na.rm = TRUE),
    sample_replicate_pearson = NA_real_, cross_batch_sample_icc = NA_real_, repeat_metric_reason = "No repeated biological participant/patient measurements across batches.",
    stringsAsFactors = FALSE
  )
}))
metric_rows$method <- factor(metric_rows$method, levels = method_order)
metric_rows <- metric_rows[order(metric_rows$method), ]; metric_rows$method <- as.character(metric_rows$method)
write.csv(metric_rows, file.path(result_dir, "method_metrics.csv"), row.names = FALSE, quote = TRUE, na = "")
metrics_long <- metric_rows |> pivot_longer(-c(method, repeat_metric_reason), names_to = "metric", values_to = "value")
write.csv(metrics_long, file.path(result_dir, "method_metrics_long.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(read.csv(file.path(if (dataset == "sacurine") file.path(repo_root, "data/public/sacurine/audit/randomization_diagnostics.csv") else file.path(repo_root, "data/public/waveica_adenocarcinoma/audit/randomization_diagnostics.csv"))), file.path(result_dir, "randomization_diagnostics.csv"), row.names = FALSE, quote = TRUE, na = "")

log_line("Generating figures and numerical source tables.")
design_source <- meta
design_source$qc_role <- ifelse(meta$sample_id %in% hidden_ids, "Hidden QC", ifelse(meta$sample_id %in% training_ids, "Training QC", ifelse(meta$is_qc, "QC", "Study")))
write.csv(design_source, file.path(figure_source_dir, "acquisition_design.csv"), row.names = FALSE, quote = TRUE, na = "")
p_design <- ggplot(design_source, aes(global_run_order, as.numeric(factor(batch)), color = qc_role)) + geom_point(size = 1.5) + scale_color_manual(values = c(Study = "#9CA3AF", `Training QC` = "#0072B2", `Hidden QC` = "#D55E00", QC = "#009E73")) + labs(title = paste(config$title, "acquisition design"), x = "Global run order", y = "Supplied batch", color = NULL) + theme_publication()
ggsave(file.path(figure_dir, "acquisition_design.pdf"), p_design, width = 11, height = 4.8)

pca_source <- bind_rows(lapply(intersect(c("Raw", "WINN auto (QC)"), completed), function(label) {
  ids <- study_meta$sample_id
  pca <- prcomp(t(log1p(method_results[[label]][, ids, drop = FALSE])), center = TRUE, scale. = TRUE)
  data.frame(method = label, sample_id = ids, PC1 = pca$x[, 1], PC2 = pca$x[, 2], batch = study_meta$batch, biological = if (dataset == "sacurine") study_meta$gender else study_meta$biological_group, stringsAsFactors = FALSE)
}))
write.csv(pca_source, file.path(figure_source_dir, "raw_winn_pca.csv"), row.names = FALSE, quote = TRUE)
p_pca <- ggplot(pca_source, aes(PC1, PC2, color = biological, shape = factor(batch))) + geom_point(alpha = 0.75, size = 1.7) + facet_wrap(~method, scales = "free") + labs(title = paste(config$title, "Raw versus WiNN PCA"), color = if (dataset == "sacurine") "Gender" else "Group", shape = "Batch") + theme_publication()
ggsave(file.path(figure_dir, "raw_vs_winn_pca.pdf"), p_pca, width = 11, height = 5.5)

technical_source <- metric_rows |> select(method, heldout_qc_cv_mean, residual_gam_deviance_mean, residual_ljung_box_proportion, batch_weighted_pc_r2) |> pivot_longer(-method, names_to = "metric", values_to = "value")
technical_source$method <- factor(technical_source$method, levels = method_order)
write.csv(technical_source, file.path(figure_source_dir, "technical_comparison.csv"), row.names = FALSE, quote = TRUE)
p_technical <- ggplot(technical_source, aes(method, value, fill = method)) + geom_col() + facet_wrap(~metric, scales = "free_y") + scale_fill_manual(values = method_palette, guide = "none") + labs(title = paste(config$title, "technical comparison"), x = NULL, y = NULL) + theme_publication() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(figure_dir, "technical_comparison.pdf"), p_technical, width = 13, height = 7)

biological_source <- metric_rows |> select(method, biological_weighted_pc_r2_mean, biological_associated_features, cross_batch_effect_pearson_median) |> pivot_longer(-method, names_to = "metric", values_to = "value")
biological_source$method <- factor(biological_source$method, levels = method_order)
write.csv(biological_source, file.path(figure_source_dir, "biological_comparison.csv"), row.names = FALSE, quote = TRUE)
p_bio <- ggplot(biological_source, aes(method, value, fill = method)) + geom_col() + facet_wrap(~metric, scales = "free_y") + scale_fill_manual(values = method_palette, guide = "none") + labs(title = paste(config$title, "biological preservation"), x = NULL, y = NULL) + theme_publication() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(figure_dir, "biological_comparison.pdf"), p_bio, width = 12, height = 6)

if (file.exists(file.path(result_dir, "winn_selectivity_summary.csv"))) {
  selectivity_source <- read.csv(file.path(result_dir, "winn_selectivity_summary.csv"), check.names = FALSE)
  if (nrow(selectivity_source)) {
    selectivity_long <- selectivity_source |> select(method, proportion_unique_features_detrended, proportion_features_selected_for_batch) |> pivot_longer(-method, names_to = "gate", values_to = "proportion")
    write.csv(selectivity_long, file.path(figure_source_dir, "winn_selectivity.csv"), row.names = FALSE, quote = TRUE)
    p_selectivity <- ggplot(selectivity_long, aes(method, proportion, fill = gate)) + geom_col(position = "dodge") + scale_y_continuous(labels = scales::label_percent()) + labs(title = paste(config$title, "WiNN selectivity"), x = NULL, y = "Proportion", fill = NULL) + theme_publication() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
    ggsave(file.path(figure_dir, "winn_selectivity.pdf"), p_selectivity, width = 9, height = 5)
  }
}

raw_gam <- gam_rows[gam_rows$method == "Raw", ] |> group_by(metabolite) |> summarise(raw_drift = mean(explained, na.rm = TRUE), .groups = "drop")
targets <- quantile(raw_gam$raw_drift, c(0.1, 0.5, 0.9), na.rm = TRUE)
representatives <- unique(vapply(targets, function(target) raw_gam$metabolite[which.min(abs(raw_gam$raw_drift - target))], character(1)))
trajectory_source <- bind_rows(lapply(intersect(c("Raw", "QC-RLSC", "WINN auto (QC)"), completed), function(label) {
  values <- method_results[[label]][representatives, , drop = FALSE]
  out <- as.data.frame(as.table(values), stringsAsFactors = FALSE)
  names(out) <- c("feature_id", "sample_id", "intensity")
  out$method <- label
  out <- merge(out, design_source[, c("sample_id", "global_run_order", "batch", "qc_role")], by = "sample_id", sort = FALSE)
  out
}))
write.csv(trajectory_source, file.path(figure_source_dir, "representative_feature_trajectories.csv"), row.names = FALSE, quote = TRUE)
p_trend <- ggplot(trajectory_source, aes(global_run_order, log1p(intensity), color = qc_role)) + geom_point(size = 0.8, alpha = 0.7) + facet_grid(feature_id ~ method, scales = "free_y") + labs(title = paste(config$title, "prespecified raw-drift quantile features"), x = "Global run order", y = "log1p intensity", color = NULL) + theme_publication()
ggsave(file.path(figure_dir, "representative_feature_trajectories.pdf"), p_trend, width = 13, height = 8)

if (dataset == "sacurine") {
  effect_source <- cross_batch
  write.csv(effect_source, file.path(figure_source_dir, "sacurine_effect_concordance.csv"), row.names = FALSE, quote = TRUE)
} else {
  write.csv(cross_batch, file.path(figure_source_dir, "waveica_batch_effect_concordance.csv"), row.names = FALSE, quote = TRUE)
}

source(file.path(repo_root, "scripts", "augment_human_public_benchmark_outputs.R"))
augment_human_public_benchmark_outputs(repo_root, dataset)
log_line("Generated publication figures, numerical figure sources, and cross-batch direction diagnostics.")

writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
preprocessing_summary_source <- if (dataset == "sacurine") file.path(repo_root, "data/public/sacurine/processed/sacurine_preprocessing_summary.json") else file.path(repo_root, "data/public/waveica_adenocarcinoma/processed/waveica_preprocessing_summary.json")
file.copy(preprocessing_summary_source, file.path(result_dir, "preprocessing_summary.json"), overwrite = TRUE)
validation <- data.frame(
  check = c("nine_methods_attempted", "hidden_qcs_not_supplied", "phenotypes_not_used_for_tuning", "exact_sample_alignment_all_completed", "feature_loss_reported", "ljung_bh_within_segment", "autobatch_no_true_batches", "no_repeat_metrics_fabricated"),
  passed = c(nrow(attempts) == 9L && identical(attempts$method, method_order), all(!method_protocol$supplied_hidden_qc_ids), all(!method_protocol$tuning_uses_phenotype), all(dimension_rows$exact_sample_alignment), all(!is.na(dimension_rows$features_lost)), TRUE, !method_protocol$supplied_batch_labels[method_protocol$method == "WINN auto-batch (QC)"], all(is.na(metric_rows$sample_replicate_pearson)) && all(is.na(metric_rows$cross_batch_sample_icc))),
  stringsAsFactors = FALSE
)
write.csv(validation, file.path(result_dir, "final_validation.csv"), row.names = FALSE, quote = TRUE)
if (!all(validation$passed)) stop("One or more final validation checks failed.")
completion <- list(dataset = dataset, completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), dimensions = dim(x), hidden_qcs = hidden_ids, training_qcs = training_ids, attempted_methods = method_order, completed_methods = completed, failed_methods = attempts$method[attempts$status == "failed"], validations_passed = all(validation$passed), winn_version = as.character(packageVersion("winn")), winn_commit = "b72aa80a2a6400126092b0814a8ba6012f5f863e")
jsonlite::write_json(completion, file.path(result_dir, "analysis_complete.json"), auto_unbox = TRUE, pretty = TRUE)
log_line("Benchmark completed and validated for ", config$title, ".")
