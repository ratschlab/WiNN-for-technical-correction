#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript or R --file.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args

required_packages <- c(
  "dplyr", "ggplot2", "jsonlite", "lmtest", "malbacR", "mgcv",
  "qcrlscR", "scales", "statTarget", "sva", "tibble", "tidyr",
  "TIGERr", "winn"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing benchmark package(s): ", paste(missing_packages, collapse = ", "),
    ". Use the R 4.5 benchmark runtime documented in README.md."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(tidyr)
})

source(file.path(repo_root, "scripts", "weighted_pc_r2.R"), local = TRUE)
source(file.path(repo_root, "scripts", "benchmark_helpers.R"), local = TRUE)
source(file.path(repo_root, "scripts", "run_order_drift_helpers.R"), local = TRUE)
source(file.path(repo_root, "scripts", "plot_theme.R"), local = TRUE)

processed_dir <- file.path(repo_root, "data", "public", "batchcorr_set1", "processed")
result_dir <- file.path(repo_root, "results", "batchcorr_set1")
figure_dir <- file.path(result_dir, "figures")
figure_data_dir <- file.path(figure_dir, "source_data")
tuning_dir <- file.path(result_dir, "tuning_candidates")
method_output_dir <- file.path(result_dir, "method_outputs")
if (force && dir.exists(result_dir)) {
  unlink(result_dir, recursive = TRUE, force = TRUE)
}
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tuning_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(method_output_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(result_dir, "analysis_log.txt")
if (force || !file.exists(log_path)) writeLines(character(), log_path)
log_line <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "] ", paste0(..., collapse = ""))
  cat(line, "\n", file = log_path, append = TRUE)
  message(line)
}

read_feature_matrix <- function(path) {
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"feature_id" %in% names(d)) stop("No feature_id column in ", path)
  feature_id <- d$feature_id
  d$feature_id <- NULL
  m <- as.matrix(d)
  storage.mode(m) <- "double"
  rownames(m) <- feature_id
  m
}

metadata_path <- file.path(processed_dir, "BatchCorr_set1_metadata.csv")
matrix_path <- file.path(processed_dir, "BatchCorr_set1_imputed_data.csv")
if (!file.exists(metadata_path) || !file.exists(matrix_path)) {
  stop("Preprocessed Set 1 files are missing. Run scripts/preprocess_batchcorr_set1.R first.")
}

meta <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
x <- read_feature_matrix(matrix_path)
logical_cols <- c("is_qc", "is_reference", "is_study_sample")
meta[logical_cols] <- lapply(meta[logical_cols], as.logical)
meta$run_order <- as.integer(meta$run_order)
meta$within_batch_order <- as.integer(meta$within_batch_order)
meta$batch_order_index <- as.integer(meta$batch_order_index)
meta <- meta[order(meta$run_order), , drop = FALSE]
x <- x[, meta$sample_id, drop = FALSE]

if (anyNA(x) || any(!is.finite(x)) || any(x < 0)) stop("Analysis matrix is not finite and nonnegative.")
if (!identical(colnames(x), meta$sample_id)) stop("Analysis matrix does not align with metadata.")
if (anyDuplicated(colnames(x)) || anyDuplicated(rownames(x))) stop("Duplicate matrix identifiers found.")
if (anyNA(meta$batch) || anyNA(meta$run_order)) stop("Missing batch or run order.")
log_line("Loaded ", nrow(x), " features x ", ncol(x), " injections.")

method_order <- winn_method_order()
method_palette <- winn_method_palette()

# Fixed fresh hidden-QC protocol for the amended benchmark: exactly one
# interior QC per analytical batch. Seed 20260809 selects no IDs used by the
# superseded implementation audit.
qc_meta <- meta[meta$is_qc, , drop = FALSE]
batch_order <- unique(meta$batch[order(meta$run_order)])
holdout_seed <- 20260809L
set.seed(holdout_seed)
holdout_ids <- unlist(lapply(batch_order, function(batch_id) {
  d <- qc_meta[qc_meta$batch == batch_id, , drop = FALSE]
  d <- d[order(d$run_order), , drop = FALSE]
  eligible <- d$sample_id[seq_len(nrow(d)) %in% seq.int(2L, nrow(d) - 1L)]
  if (length(eligible) < 1L) stop("No interior QC is available in batch ", batch_id)
  sample(eligible, 1L)
}), use.names = FALSE)
holdout_ids <- meta$sample_id[meta$sample_id %in% holdout_ids]
train_qc_ids <- meta$sample_id[meta$is_qc & !(meta$sample_id %in% holdout_ids)]

holdout_by_batch <- table(meta$batch[match(holdout_ids, meta$sample_id)])
train_by_batch <- table(meta$batch[match(train_qc_ids, meta$sample_id)])
if (!all(batch_order %in% names(holdout_by_batch)) || any(holdout_by_batch[batch_order] != 1L)) {
  stop("Holdout protocol did not select exactly one QC per batch.")
}
if (any(train_by_batch[batch_order] < 3L)) stop("Fewer than three training QCs remain in a batch.")
for (batch_id in batch_order) {
  q <- qc_meta[qc_meta$batch == batch_id, , drop = FALSE]
  q <- q[order(q$run_order), , drop = FALSE]
  if (q$sample_id[1] %in% holdout_ids || q$sample_id[nrow(q)] %in% holdout_ids) {
    stop("An endpoint QC was selected for holdout in batch ", batch_id)
  }
}

holdout_table <- meta[match(holdout_ids, meta$sample_id), c("sample_id", "batch", "run_order", "within_batch_order")]
holdout_table$seed <- holdout_seed
holdout_table <- holdout_table[, c("seed", "sample_id", "batch", "run_order", "within_batch_order")]
write.csv(holdout_table, file.path(result_dir, "qc_holdout_ids.csv"), row.names = FALSE, quote = TRUE)
training_table <- meta[match(train_qc_ids, meta$sample_id), c("sample_id", "batch", "run_order", "within_batch_order")]
write.csv(training_table, file.path(result_dir, "qc_training_ids.csv"), row.names = FALSE, quote = TRUE)

protocol_lines <- c(
  "# BatchCorr Set 1 hidden-QC protocol",
  "",
  paste0("Ten QCs are held out: one reproducibly sampled interior QC from each of the ten analytical batches (seed ", holdout_seed, ")."),
  "The first and last QC in every batch remain training controls, and every batch retains at least four training QCs.",
  "",
  "Before any method call, each held-out QC is relabelled as an ordinary Sample in the method-facing metadata.",
  "Only the 41 training QC IDs are supplied to QC-aware methods and WiNN auto-selection.",
  "QC-RLSC uses a Set 1-specific sample-ID-preserving implementation: subtractive per-batch LOESS on log1p intensities, degree 1, span 1, followed by batch shifting on the log scale.",
  "QC-RLSC parameters are fixed from the 4-5-training-QC-per-batch design and are not selected using held-out values.",
  "TIGER and SERRF use run_tiger_all_corrected/run_serrf_all_corrected with the sanitized metadata and training controls only.",
  "The held-out IDs are used only after correction for evaluation; no tuning score uses held-out QC values or labels.",
  "ComBat receives batch labels but no QC identities. WiNN auto-batch receives no supplied batch labels. WiNN fixed receives supplied batches and no QC identities."
)
writeLines(protocol_lines, file.path(result_dir, "qc_holdout_protocol.md"))

meta_hidden <- meta
meta_hidden$class[meta_hidden$sample_id %in% holdout_ids] <- "Sample"
meta_hidden$sample_type[meta_hidden$sample_id %in% holdout_ids] <- "Study"
if (any(meta_hidden$class[match(holdout_ids, meta_hidden$sample_id)] == "QC")) {
  stop("Held-out QC labels remain visible in method-facing metadata.")
}

method_protocol <- data.frame(
  method = method_order,
  supplied_batch_labels = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  supplied_training_qc_ids = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  supplied_holdout_qc_ids = FALSE,
  method_facing_holdout_class = "Sample",
  tuning_uses_holdout_values = FALSE,
  stringsAsFactors = FALSE
)
write.csv(method_protocol, file.path(result_dir, "hidden_qc_method_protocol.csv"), row.names = FALSE, quote = TRUE)

# qcrlscR 0.1.3's qc.rlsc.wrap(intra = TRUE) concatenates batches in factor
# level order and does not forward fixed span arguments to qc.rlsc(). Both
# behaviours can silently invalidate a sample-indexed benchmark. Set 1 is also
# deposited on a log scale and retains only 4-5 training QCs per batch after
# holdout. This dataset-specific implementation therefore performs stable
# degree-1, full-span subtractive correction on the exact log1p analysis scale,
# preserves input row IDs explicitly, and shifts batches on that same scale.
run_qc_rlsc_batchcorr_set1 <- function(
  x_mat,
  control_ids,
  meta_df,
  span = 1,
  degree = 1
) {
  md <- align_meta_to_mat(x_mat, meta_df) |>
    arrange(run_order)
  x_ord <- x_mat[, md$sample_id, drop = FALSE]
  dat_log <- as.data.frame(t(log1p(x_ord)), check.names = FALSE)
  qc_labels <- ifelse(md$sample_id %in% control_ids, "QC", "Sample")
  batch_levels <- unique(as.character(md$batch))
  batch_factor <- factor(md$batch, levels = batch_levels)
  training_counts <- table(batch_factor[md$sample_id %in% control_ids])
  if (length(training_counts) != length(batch_levels) || any(training_counts < 4L)) {
    stop("Corrected QC-RLSC requires at least four training QCs in every batch.")
  }

  corrected_log <- matrix(
    NA_real_, nrow = nrow(dat_log), ncol = ncol(dat_log),
    dimnames = dimnames(dat_log)
  )
  for (batch_id in batch_levels) {
    idx <- which(md$batch == batch_id)
    batch_corrected <- qcrlscR::qc.rlsc(
      x = dat_log[idx, , drop = FALSE],
      y = qc_labels[idx],
      method = "subtract",
      opti = FALSE,
      span = span,
      degree = degree
    )
    batch_corrected <- as.matrix(batch_corrected)
    if (!identical(rownames(batch_corrected), rownames(dat_log)[idx])) {
      stop("QC-RLSC changed sample identifiers inside batch ", batch_id, ".")
    }
    corrected_log[idx, ] <- batch_corrected
  }
  if (anyNA(corrected_log) || any(!is.finite(corrected_log))) {
    stop("Corrected QC-RLSC produced non-finite within-batch values.")
  }

  corrected_log <- as.matrix(qcrlscR::batch.shift(
    as.data.frame(corrected_log, check.names = FALSE),
    batch_factor,
    overall_average = TRUE
  ))
  if (!identical(rownames(corrected_log), md$sample_id)) {
    stop("QC-RLSC batch shifting changed sample order or identifiers.")
  }
  corrected <- t(expm1(corrected_log))
  rownames(corrected) <- rownames(x_ord)
  colnames(corrected) <- md$sample_id
  corrected <- clip_nonnegative(corrected)
  corrected[, meta_df$sample_id, drop = FALSE]
}

sample_replicate_metrics <- function(mat, md) {
  ids <- intersect(md$sample_id[md$is_study_sample], colnames(mat))
  md <- md[match(ids, md$sample_id), , drop = FALSE]
  z <- log1p(mat[, ids, drop = FALSE])
  groups <- split(seq_along(ids), md$accession_id)
  groups <- groups[lengths(groups) >= 2L]
  bind_rows(lapply(names(groups), function(accession) {
    idx <- groups[[accession]]
    pair_index <- utils::combn(idx, 2L, simplify = FALSE)
    pair_values <- bind_rows(lapply(pair_index, function(p) {
      data.frame(
        pearson = suppressWarnings(stats::cor(z[, p[1]], z[, p[2]], method = "pearson")),
        icc_a1 = calc_icc_a1_matrix(cbind(z[, p[1]], z[, p[2]])),
        stringsAsFactors = FALSE
      )
    }))
    data.frame(
      accession_id = accession,
      n_replicates = length(idx),
      n_pairs = nrow(pair_values),
      pearson = median(pair_values$pearson, na.rm = TRUE),
      icc_a1 = median(pair_values$icc_a1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

feature_repeatability <- function(mat, md) {
  ids <- intersect(md$sample_id[md$is_study_sample], colnames(mat))
  md <- md[match(ids, md$sample_id), , drop = FALSE]
  z <- log1p(mat[, ids, drop = FALSE])
  groups <- split(seq_along(ids), md$accession_id)
  groups <- groups[lengths(groups) >= 2L]
  out <- vapply(seq_len(nrow(z)), function(feature_index) {
    y <- z[feature_index, ]
    means <- vapply(groups, function(idx) mean(y[idx]), numeric(1))
    between <- stats::var(means)
    within_num <- sum(vapply(groups, function(idx) {
      if (length(idx) < 2L) return(0)
      stats::var(y[idx]) * (length(idx) - 1L)
    }, numeric(1)))
    within_df <- sum(lengths(groups) - 1L)
    within <- if (within_df > 0L) within_num / within_df else NA_real_
    if (!is.finite(between) || !is.finite(within) || between + within <= 0) NA_real_ else between / (between + within)
  }, numeric(1))
  data.frame(feature_id = rownames(z), repeatability = out, stringsAsFactors = FALSE)
}

genotype_associations <- function(mat, md, method_label) {
  ids <- intersect(md$sample_id[md$is_study_sample], colnames(mat))
  md <- md[match(ids, md$sample_id), , drop = FALSE]
  g <- factor(md$accession_id)
  z <- log1p(mat[, ids, drop = FALSE])
  group_index <- split(seq_along(g), g)
  n <- length(g)
  k <- nlevels(g)
  rows <- lapply(seq_len(nrow(z)), function(i) {
    y <- z[i, ]
    grand <- mean(y)
    group_means <- vapply(group_index, function(idx) mean(y[idx]), numeric(1))
    ss_between <- sum(lengths(group_index) * (group_means - grand)^2)
    ss_total <- sum((y - grand)^2)
    ss_within <- max(0, ss_total - ss_between)
    f_stat <- if (k > 1L && n > k && ss_within > 0) {
      (ss_between / (k - 1L)) / (ss_within / (n - k))
    } else {
      NA_real_
    }
    p_value <- if (is.finite(f_stat)) stats::pf(f_stat, k - 1L, n - k, lower.tail = FALSE) else NA_real_
    data.frame(
      method = method_label,
      feature_id = rownames(z)[i],
      n_injections = n,
      n_accessions = k,
      f_statistic = f_stat,
      p_value = p_value,
      eta_squared = if (ss_total > 0) ss_between / ss_total else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- bind_rows(rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out$is_associated <- is.finite(out$p_adj) & out$p_adj < 0.05
  out
}

evaluate_tuning_matrix <- function(mat, raw_subset) {
  study_ids <- meta$sample_id[meta$is_study_sample & meta$sample_id %in% colnames(mat)]
  if (length(study_ids) < 100L) stop("Insufficient study-sample coverage during tuning.")
  m <- mat[, study_ids, drop = FALSE]
  md <- meta[match(study_ids, meta$sample_id), , drop = FALSE]
  reps <- sample_replicate_metrics(m, md)
  feature_rep <- feature_repeatability(m, md)
  associations <- genotype_associations(m, md, "candidate")
  data.frame(
    output_features = nrow(m),
    output_study_samples = ncol(m),
    accession_associated_features = sum(associations$is_associated, na.rm = TRUE),
    accession_pc_r2 = weighted_pc_r2(m, md$accession_id, target_type = "categorical", n_pcs = 5L),
    batch_pc_r2 = weighted_pc_r2(m, md$batch, target_type = "categorical", n_pcs = 5L),
    replicate_pearson_median = median(reps$pearson, na.rm = TRUE),
    replicate_icc_median = median(reps$icc_a1, na.rm = TRUE),
    feature_repeatability_median = median(feature_rep$repeatability, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

rank_metric <- function(x, higher_better = TRUE) {
  x <- as.numeric(x)
  if (higher_better) {
    x[!is.finite(x)] <- -Inf
    rank(x, ties.method = "average")
  } else {
    x[!is.finite(x)] <- Inf
    rank(-x, ties.method = "average")
  }
}

score_candidates <- function(d) {
  if (!nrow(d)) return(d)
  d$score <-
    rank_metric(d$accession_associated_features, TRUE) +
    rank_metric(d$accession_pc_r2, TRUE) +
    rank_metric(d$replicate_pearson_median, TRUE) +
    rank_metric(d$replicate_icc_median, TRUE) +
    rank_metric(d$feature_repeatability_median, TRUE) +
    rank_metric(d$batch_pc_r2, FALSE)
  d[order(-d$score, d$candidate_id), , drop = FALSE]
}

capture_call <- function(fn) {
  warnings <- character()
  messages <- character()
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(
      fn(),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    value = value,
    elapsed_sec = elapsed,
    warnings = unique(warnings),
    messages = unique(messages),
    error = if (inherits(value, "error")) conditionMessage(value) else ""
  )
}

set.seed(20260722)
tune_feature_ids <- sample(rownames(x), min(120L, nrow(x)), replace = FALSE)
x_tune <- x[tune_feature_ids, , drop = FALSE]
write.csv(data.frame(feature_id = tune_feature_ids), file.path(tuning_dir, "tuning_feature_ids.csv"), row.names = FALSE, quote = TRUE)

tuning_grids <- list(
  "ComBat" = data.frame(par_prior = c(TRUE, FALSE)),
  "QC-RFSC" = expand.grid(ntree = c(200L, 500L), coCV = c(20, 30), Frule = 0.8, KEEP.OUT.ATTRS = FALSE),
  "TIGER" = data.frame(use_injection = c(TRUE, FALSE)),
  "SERRF" = data.frame(jitter_eps = c(0, 1e-6, 1e-5))
)

run_candidate_method <- function(method_label, p, x_mat) {
  switch(
    method_label,
    "ComBat" = run_combat(x_mat, meta_df = meta_hidden, par_prior = as.logical(p$par_prior)),
    "QC-RFSC" = run_qc_rfsc_with_controls(
      x_mat,
      control_ids = train_qc_ids,
      meta_df = meta_hidden,
      ntree = as.integer(p$ntree),
      coCV = as.numeric(p$coCV),
      Frule = as.numeric(p$Frule)
    ),
    "TIGER" = run_tiger_all_corrected(
      x_mat,
      control_ids = train_qc_ids,
      meta_df = meta_hidden,
      use_injection = as.logical(p$use_injection),
      mtry_percent = 0.4,
      nodesize_percent = 0.4,
      ntree = 5,
      parallel_cores = 1
    ),
    "SERRF" = run_serrf_all_corrected(
      x_mat,
      control_ids = train_qc_ids,
      meta_df = meta_hidden,
      jitter_eps = as.numeric(p$jitter_eps)
    ),
    stop("Unknown tuning method: ", method_label)
  )
}

best_params <- list()
tuning_seconds <- setNames(rep(0, length(tuning_grids)), names(tuning_grids))
tuning_diagnostics <- list()

for (method_label in names(tuning_grids)) {
  grid <- tuning_grids[[method_label]]
  candidate_rows <- vector("list", nrow(grid))
  log_line("Tuning ", method_label, " across ", nrow(grid), " candidates.")
  method_tuning_start <- proc.time()[["elapsed"]]
  for (i in seq_len(nrow(grid))) {
    p <- as.list(grid[i, , drop = FALSE])
    captured <- capture_call(function() run_candidate_method(method_label, p, x_tune))
    tuning_diagnostics[[length(tuning_diagnostics) + 1L]] <- data.frame(
      method = method_label,
      candidate_id = i,
      warnings = paste(captured$warnings, collapse = " | "),
      messages = paste(captured$messages, collapse = " | "),
      error = captured$error,
      stringsAsFactors = FALSE
    )
    if (!inherits(captured$value, "error") && is.matrix(captured$value)) {
      eval <- tryCatch(evaluate_tuning_matrix(clip_nonnegative(captured$value), x_tune), error = function(e) e)
      if (!inherits(eval, "error")) {
        candidate_rows[[i]] <- cbind(
          data.frame(candidate_id = i, grid[i, , drop = FALSE], stringsAsFactors = FALSE),
          runtime_sec = captured$elapsed_sec,
          status = "completed",
          eval
        )
      } else {
        captured$error <- conditionMessage(eval)
      }
    }
    if (is.null(candidate_rows[[i]])) {
      candidate_rows[[i]] <- cbind(
        data.frame(candidate_id = i, grid[i, , drop = FALSE], stringsAsFactors = FALSE),
        runtime_sec = captured$elapsed_sec,
        status = "failed",
        output_features = NA_integer_,
        output_study_samples = NA_integer_,
        accession_associated_features = NA_integer_,
        accession_pc_r2 = NA_real_,
        batch_pc_r2 = NA_real_,
        replicate_pearson_median = NA_real_,
        replicate_icc_median = NA_real_,
        feature_repeatability_median = NA_real_
      )
    }
    log_line(method_label, " candidate ", i, "/", nrow(grid), ": ", candidate_rows[[i]]$status)
  }
  tuning_seconds[[method_label]] <- proc.time()[["elapsed"]] - method_tuning_start
  candidates <- score_candidates(bind_rows(candidate_rows))
  write.csv(candidates, file.path(tuning_dir, paste0(gsub("[^A-Za-z0-9]+", "_", tolower(method_label)), ".csv")), row.names = FALSE, quote = TRUE, na = "")
  valid <- candidates[candidates$status == "completed" & is.finite(candidates$score), , drop = FALSE]
  if (!nrow(valid)) {
    best_params[[method_label]] <- grid[1, , drop = FALSE]
    best_params[[method_label]]$selection_status <- "fallback_after_all_candidates_failed"
  } else {
    selected_id <- valid$candidate_id[1]
    best_params[[method_label]] <- grid[selected_id, , drop = FALSE]
    best_params[[method_label]]$selection_status <- "selected_by_ranked_biological_preservation_and_batch_metrics"
  }
}
write.csv(bind_rows(tuning_diagnostics), file.path(tuning_dir, "candidate_diagnostics.csv"), row.names = FALSE, quote = TRUE)

selected_param_rows <- bind_rows(lapply(names(best_params), function(label) {
  data.frame(method = label, best_params[[label]], check.names = FALSE, stringsAsFactors = FALSE)
}))
write.csv(selected_param_rows, file.path(tuning_dir, "selected_competitor_parameters.csv"), row.names = FALSE, quote = TRUE, na = "")

method_messages <- list()
method_attempts <- list()
method_results <- list()
cache_version <- "batchcorr_set1_freshqc10_qcrlsc_corrected_v2"

run_full_method <- function(label, fn, tuning_sec = 0) {
  log_line("Starting full method: ", label)
  captured <- capture_call(function() run_method_cached(
    label = label,
    fn = fn,
    plot_dir = result_dir,
    tuning_sec = tuning_sec,
    force = force,
    cache_version = cache_version
  ))
  method_messages[[length(method_messages) + 1L]] <<- data.frame(
    method = label,
    warnings = paste(captured$warnings, collapse = " | "),
    messages = paste(captured$messages, collapse = " | "),
    error = captured$error,
    stringsAsFactors = FALSE
  )
  if (inherits(captured$value, "error")) {
    method_attempts[[length(method_attempts) + 1L]] <<- data.frame(
      method = label,
      status = "failed",
      tuning_sec = tuning_sec,
      process_sec = captured$elapsed_sec,
      total_sec = tuning_sec + captured$elapsed_sec,
      error = captured$error,
      stringsAsFactors = FALSE
    )
    log_line("Method failed: ", label, " — ", captured$error)
    return(NULL)
  }
  validated <- tryCatch({
    result <- captured$value
    mat <- result$data
    if (!is.matrix(mat)) stop(label, " did not return a matrix.")
    if (is.null(rownames(mat)) && nrow(mat) == nrow(x)) rownames(mat) <- rownames(x)
    if (is.null(colnames(mat)) && ncol(mat) == ncol(x)) colnames(mat) <- colnames(x)
    if (anyDuplicated(rownames(mat)) || anyDuplicated(colnames(mat))) stop(label, " returned duplicate identifiers.")
    if (label == "TIGER") {
      missing_ids <- setdiff(colnames(x), colnames(mat))
      if (length(missing_ids)) mat <- cbind(mat, x[, missing_ids, drop = FALSE])
    }
    common_features <- rownames(x)[rownames(x) %in% rownames(mat)]
    common_samples <- meta$sample_id[meta$sample_id %in% colnames(mat)]
    if (!length(common_features) || !length(common_samples)) stop(label, " returned no alignable feature/sample identifiers.")
    mat <- mat[common_features, common_samples, drop = FALSE]
    if (anyNA(mat) || any(!is.finite(mat))) stop(label, " returned non-finite values.")
    result$data <- clip_nonnegative(mat)
    result
  }, error = function(e) e)
  if (inherits(validated, "error")) {
    validation_error <- conditionMessage(validated)
    method_messages[[length(method_messages)]]$error <<- validation_error
    method_attempts[[length(method_attempts) + 1L]] <<- data.frame(
      method = label,
      status = "failed",
      tuning_sec = tuning_sec,
      process_sec = captured$elapsed_sec,
      total_sec = tuning_sec + captured$elapsed_sec,
      error = validation_error,
      stringsAsFactors = FALSE
    )
    log_line("Method output validation failed: ", label, " — ", validation_error)
    return(NULL)
  }
  result <- validated
  mat <- result$data
  saveRDS(mat, file.path(method_output_dir, paste0(gsub("[^A-Za-z0-9]+", "_", tolower(label)), ".rds")))
  method_attempts[[length(method_attempts) + 1L]] <<- data.frame(
    method = label,
    status = "completed",
    tuning_sec = result$tuning_sec,
    process_sec = result$process_sec,
    total_sec = result$elapsed_sec,
    error = "",
    stringsAsFactors = FALSE
  )
  log_line("Completed full method: ", label, " (", nrow(mat), " x ", ncol(mat), ").")
  result
}

method_results[["Raw"]] <- run_full_method("Raw", function() x)
method_results[["ComBat"]] <- run_full_method(
  "ComBat",
  function() run_combat(x, meta_df = meta_hidden, par_prior = as.logical(best_params[["ComBat"]]$par_prior[1])),
  tuning_sec = tuning_seconds[["ComBat"]]
)
method_results[["QC-RLSC"]] <- run_full_method(
  "QC-RLSC",
  function() run_qc_rlsc_batchcorr_set1(x, train_qc_ids, meta_hidden, span = 1, degree = 1),
  tuning_sec = 0
)
method_results[["QC-RFSC"]] <- run_full_method(
  "QC-RFSC",
  function() run_qc_rfsc_with_controls(
    x,
    train_qc_ids,
    meta_hidden,
    ntree = as.integer(best_params[["QC-RFSC"]]$ntree[1]),
    coCV = as.numeric(best_params[["QC-RFSC"]]$coCV[1]),
    Frule = as.numeric(best_params[["QC-RFSC"]]$Frule[1])
  ),
  tuning_sec = tuning_seconds[["QC-RFSC"]]
)
method_results[["TIGER"]] <- run_full_method(
  "TIGER",
  function() run_tiger_all_corrected(
    x,
    train_qc_ids,
    meta_hidden,
    use_injection = as.logical(best_params[["TIGER"]]$use_injection[1]),
    mtry_percent = 0.4,
    nodesize_percent = 0.4,
    ntree = 5,
    parallel_cores = 1
  ),
  tuning_sec = tuning_seconds[["TIGER"]]
)
method_results[["SERRF"]] <- run_full_method(
  "SERRF",
  function() run_serrf_all_corrected(x, train_qc_ids, meta_hidden, jitter_eps = as.numeric(best_params[["SERRF"]]$jitter_eps[1])),
  tuning_sec = tuning_seconds[["SERRF"]]
)
method_results[["WINN auto (QC)"]] <- run_full_method(
  "WINN auto (QC)",
  function() run_winn_mode1_auto_qc(x, train_qc_ids, meta_hidden, params = NULL)
)
method_results[["WINN auto-batch (QC)"]] <- run_full_method(
  "WINN auto-batch (QC)",
  function() run_winn_mode2_auto_batch_qc(x, train_qc_ids, meta_hidden, params = NULL)
)
method_results[["WINN default (no QC)"]] <- run_full_method(
  "WINN default (no QC)",
  function() run_winn_mode3_default_no_qc(x, meta_hidden)
)

attempts_df <- bind_rows(method_attempts)
attempts_df$method <- factor(attempts_df$method, levels = method_order)
attempts_df <- attempts_df[order(attempts_df$method), , drop = FALSE]
attempts_df$method <- as.character(attempts_df$method)
write.csv(attempts_df, file.path(result_dir, "method_attempts.csv"), row.names = FALSE, quote = TRUE)
write.csv(attempts_df[, c("method", "tuning_sec", "process_sec", "total_sec")], file.path(result_dir, "runtime.csv"), row.names = FALSE, quote = TRUE)
messages_df <- bind_rows(method_messages)
write.csv(messages_df, file.path(result_dir, "method_warnings_messages_errors.csv"), row.names = FALSE, quote = TRUE)
write.csv(attempts_df[attempts_df$status == "failed", c("method", "error"), drop = FALSE], file.path(result_dir, "method_failures.csv"), row.names = FALSE, quote = TRUE)

completed_methods <- method_order[vapply(method_results[method_order], function(z) !is.null(z), logical(1))]
method_mats <- lapply(method_results[completed_methods], `[[`, "data")

# Selected-parameter ledger, including the internally selected WiNN messages.
parameter_rows <- list(
  data.frame(method = "Raw", parameters = "none", stringsAsFactors = FALSE),
  data.frame(method = "ComBat", parameters = paste0("par_prior=", best_params[["ComBat"]]$par_prior[1]), stringsAsFactors = FALSE),
  data.frame(method = "QC-RLSC", parameters = "dataset_specific_fixed; analysis_scale=log1p; method=subtract; intra_batch=TRUE; span=1; degree=1; batch_shift=log_scale; controls=41 training QCs; sample_ids_preserved=TRUE", stringsAsFactors = FALSE),
  data.frame(method = "QC-RFSC", parameters = paste0("ntree=", best_params[["QC-RFSC"]]$ntree[1], "; coCV=", best_params[["QC-RFSC"]]$coCV[1], "; Frule=", best_params[["QC-RFSC"]]$Frule[1], "; controls=41 training QCs"), stringsAsFactors = FALSE),
  data.frame(method = "TIGER", parameters = paste0("use_injection=", best_params[["TIGER"]]$use_injection[1], "; mtry_percent=0.4; nodesize_percent=0.4; ntree=5; controls=41 training QCs"), stringsAsFactors = FALSE),
  data.frame(method = "SERRF", parameters = paste0("jitter_eps=", best_params[["SERRF"]]$jitter_eps[1], "; controls=41 training QCs"), stringsAsFactors = FALSE),
  data.frame(method = "WINN auto (QC)", parameters = paste(messages_df$messages[messages_df$method == "WINN auto (QC)"], collapse = " | "), stringsAsFactors = FALSE),
  data.frame(method = "WINN auto-batch (QC)", parameters = paste(messages_df$messages[messages_df$method == "WINN auto-batch (QC)"], collapse = " | "), stringsAsFactors = FALSE),
  data.frame(method = "WINN default (no QC)", parameters = "parameters=fixed; supplied batches; no QCs; test=Ljung-Box; fdr_threshold=0.05; median_adjustment=shrink; spline_method=conservative; remove_batch_effects=anova; scale_by_batch=FALSE", stringsAsFactors = FALSE)
)
parameters_df <- bind_rows(parameter_rows)
write.csv(parameters_df, file.path(result_dir, "method_parameters.csv"), row.names = FALSE, quote = TRUE)

log_line("Computing held-out QC metrics.")
qc_cv_rows <- bind_rows(lapply(completed_methods, function(label) {
  mat <- method_mats[[label]]
  qc_ids <- holdout_ids[holdout_ids %in% colnames(mat)]
  if (length(qc_ids) < 2L) return(tibble())
  q <- mat[, qc_ids, drop = FALSE]
  means <- rowMeans(q)
  cvs <- apply(q, 1, stats::sd) / means * 100
  tibble(method = label, feature_id = rownames(q), n_holdout_qc = length(qc_ids), mean = means, sd = apply(q, 1, stats::sd), cv_percent = cvs)
}))
write.csv(qc_cv_rows, file.path(result_dir, "qc_cv_by_feature.csv"), row.names = FALSE, quote = TRUE)

qc_pairwise <- bind_rows(lapply(completed_methods, function(label) {
  mat <- method_mats[[label]]
  qc_ids <- holdout_ids[holdout_ids %in% colnames(mat)]
  if (length(qc_ids) < 2L) return(tibble())
  cmat <- stats::cor(log1p(mat[, qc_ids, drop = FALSE]), method = "pearson")
  vals <- cmat[upper.tri(cmat)]
  tibble(method = label, median_pairwise_correlation = median(vals, na.rm = TRUE), mean_pairwise_correlation = mean(vals, na.rm = TRUE))
}))
write.csv(qc_pairwise, file.path(result_dir, "heldout_qc_pairwise_correlation.csv"), row.names = FALSE, quote = TRUE)

study_meta <- meta[meta$is_study_sample, , drop = FALSE]
log_line("Computing residual run-order GAM profiles.")
gam_rows <- bind_rows(lapply(completed_methods, function(label) {
  compute_metabolite_segment_gam(
    method_mats[[label]],
    study_meta,
    method_label = label,
    sample_id_col = "sample_id",
    order_col = "run_order",
    batch_col = "batch",
    transform_fun = log1p,
    min_obs = 6L,
    k_max = 6L
  )
}))
write.csv(gam_rows, file.path(result_dir, "run_order_gam_by_feature_batch.csv"), row.names = FALSE, quote = TRUE)

log_line("Computing residual Ljung-Box profiles.")
ljung_rows <- bind_rows(lapply(completed_methods, function(label) {
  compute_metabolite_segment_ljung_box(
    method_mats[[label]],
    study_meta,
    method_label = label,
    sample_id_col = "sample_id",
    order_col = "run_order",
    batch_col = "batch",
    transform_fun = log1p,
    min_obs = 8L,
    max_lag = 10L
  )
}))
ljung_rows <- annotate_metabolite_segment_ljung_box(ljung_rows, alpha = 0.05, adjust_method = "BH")
write.csv(ljung_rows, file.path(result_dir, "ljung_box_by_feature_batch.csv"), row.names = FALSE, quote = TRUE)

log_line("Computing accession associations and replicate agreement.")
association_rows <- bind_rows(lapply(completed_methods, function(label) genotype_associations(method_mats[[label]], study_meta, label)))
write.csv(association_rows, file.path(result_dir, "accession_associations.csv"), row.names = FALSE, quote = TRUE)

replicate_rows <- bind_rows(lapply(completed_methods, function(label) {
  d <- sample_replicate_metrics(method_mats[[label]], study_meta)
  d$method <- label
  d[, c("method", setdiff(names(d), "method"))]
}))
write.csv(replicate_rows, file.path(result_dir, "replicate_agreement_by_accession.csv"), row.names = FALSE, quote = TRUE)

feature_repeat_rows <- bind_rows(lapply(completed_methods, function(label) {
  d <- feature_repeatability(method_mats[[label]], study_meta)
  d$method <- label
  d[, c("method", "feature_id", "repeatability")]
}))
write.csv(feature_repeat_rows, file.path(result_dir, "feature_replicate_repeatability.csv"), row.names = FALSE, quote = TRUE)

raw_associated <- association_rows$feature_id[association_rows$method == "Raw" & association_rows$is_associated]
association_recovery <- bind_rows(lapply(completed_methods, function(label) {
  d <- association_rows[association_rows$method == label, , drop = FALSE]
  ids <- d$feature_id[d$is_associated]
  tibble(
    method = label,
    n_associated = length(ids),
    change_from_raw = length(ids) - length(raw_associated),
    overlap_with_raw = length(intersect(ids, raw_associated)),
    lost_from_raw = length(setdiff(raw_associated, ids)),
    gained_vs_raw = length(setdiff(ids, raw_associated)),
    jaccard_with_raw = if (length(union(ids, raw_associated))) length(intersect(ids, raw_associated)) / length(union(ids, raw_associated)) else NA_real_
  )
}))
write.csv(association_recovery, file.path(result_dir, "accession_association_recovery.csv"), row.names = FALSE, quote = TRUE)

gam_summary <- gam_rows |>
  group_by(method) |>
  summarise(n_profiles = sum(is.finite(explained)), run_order_gam_mean = mean(explained, na.rm = TRUE), run_order_gam_sd = sd(explained, na.rm = TRUE), .groups = "drop")
ljung_summary <- ljung_rows |>
  group_by(method) |>
  summarise(ljung_profiles_tested = sum(is.finite(p_value)), ljung_significant = sum(is_autocorrelated, na.rm = TRUE), ljung_significant_proportion = ljung_significant / ljung_profiles_tested, .groups = "drop")
qc_summary <- qc_cv_rows |>
  group_by(method) |>
  summarise(heldout_qc_cv_mean = mean(cv_percent, na.rm = TRUE), heldout_qc_cv_sd = sd(cv_percent, na.rm = TRUE), heldout_qc_cv_median = median(cv_percent, na.rm = TRUE), .groups = "drop")
rep_summary <- replicate_rows |>
  group_by(method) |>
  summarise(replicate_pearson_mean = mean(pearson, na.rm = TRUE), replicate_pearson_sd = sd(pearson, na.rm = TRUE), replicate_pearson_median = median(pearson, na.rm = TRUE), replicate_icc_median = median(icc_a1, na.rm = TRUE), .groups = "drop")
feature_rep_summary <- feature_repeat_rows |>
  group_by(method) |>
  summarise(feature_repeatability_mean = mean(repeatability, na.rm = TRUE), feature_repeatability_sd = sd(repeatability, na.rm = TRUE), feature_repeatability_median = median(repeatability, na.rm = TRUE), .groups = "drop")
assoc_summary <- association_rows |>
  group_by(method) |>
  summarise(accession_associated_features = sum(is_associated, na.rm = TRUE), accession_eta_squared_median = median(eta_squared, na.rm = TRUE), accession_eta_squared_mean = mean(eta_squared, na.rm = TRUE), .groups = "drop")

metrics_base <- bind_rows(lapply(method_order, function(label) {
  attempt <- attempts_df[attempts_df$method == label, , drop = FALSE]
  if (label %in% completed_methods) {
    mat <- method_mats[[label]]
    study_ids <- study_meta$sample_id[study_meta$sample_id %in% colnames(mat)]
    md <- study_meta[match(study_ids, study_meta$sample_id), , drop = FALSE]
    tibble(
      method = label,
      status = "completed",
      retained_features = nrow(mat),
      output_samples = ncol(mat),
      study_sample_coverage = length(study_ids) / sum(meta$is_study_sample),
      heldout_qc_coverage = sum(holdout_ids %in% colnames(mat)) / length(holdout_ids),
      batch_pc_r2 = weighted_pc_r2(mat[, study_ids, drop = FALSE], md$batch, target_type = "categorical", n_pcs = 5L),
      accession_pc_r2 = weighted_pc_r2(mat[, study_ids, drop = FALSE], md$accession_id, target_type = "categorical", n_pcs = 5L),
      tuning_sec = attempt$tuning_sec,
      process_sec = attempt$process_sec,
      runtime_sec = attempt$total_sec,
      error = ""
    )
  } else {
    tibble(
      method = label,
      status = "failed",
      retained_features = NA_integer_,
      output_samples = NA_integer_,
      study_sample_coverage = NA_real_,
      heldout_qc_coverage = NA_real_,
      batch_pc_r2 = NA_real_,
      accession_pc_r2 = NA_real_,
      tuning_sec = attempt$tuning_sec,
      process_sec = attempt$process_sec,
      runtime_sec = attempt$total_sec,
      error = attempt$error
    )
  }
}))

method_metrics <- metrics_base |>
  left_join(qc_summary, by = "method") |>
  left_join(qc_pairwise, by = "method") |>
  left_join(gam_summary, by = "method") |>
  left_join(ljung_summary, by = "method") |>
  left_join(assoc_summary, by = "method") |>
  left_join(rep_summary, by = "method") |>
  left_join(feature_rep_summary, by = "method") |>
  left_join(association_recovery, by = "method")
method_metrics$method <- factor(method_metrics$method, levels = method_order)
method_metrics <- method_metrics[order(method_metrics$method), , drop = FALSE]
method_metrics$method <- as.character(method_metrics$method)
write.csv(method_metrics, file.path(result_dir, "method_metrics.csv"), row.names = FALSE, quote = TRUE, na = "")

method_metrics_long <- method_metrics |>
  select(-status, -error) |>
  pivot_longer(cols = -method, names_to = "metric", values_to = "value")
write.csv(method_metrics_long, file.path(result_dir, "method_metrics_long.csv"), row.names = FALSE, quote = TRUE, na = "")

# Figure 1: dataset design.
design_data <- meta |>
  transmute(sample_id, batch, batch_order_index, run_order, within_batch_order, sample_type = ifelse(is_qc, "QC", "Study"), accession_id)
write.csv(design_data, file.path(figure_data_dir, "dataset_design.csv"), row.names = FALSE, quote = TRUE, na = "")
p_design <- ggplot(design_data, aes(run_order, batch_order_index, color = sample_type, shape = sample_type)) +
  geom_point(alpha = 0.75, size = 1.5) +
  scale_y_continuous(breaks = seq_along(batch_order), labels = batch_order) +
  scale_color_manual(values = c(Study = "#9CA3AF", QC = "#D55E00")) +
  labs(title = "BatchCorr Set 1 acquisition design", subtitle = "True SeqNr order; pooled references shown in orange", x = "Global acquisition sequence number", y = "Analytical batch", color = NULL, shape = NULL) +
  theme_publication()
ggsave(file.path(figure_dir, "dataset_design.pdf"), p_design, width = 12, height = 5.8)

# Figure 2: randomization and biological replication diagnostics.
batch_counts <- meta |>
  mutate(sample_type_plot = ifelse(is_qc, "QC", "Study")) |>
  count(batch, batch_order_index, sample_type_plot, name = "n")
accession_balance <- read.csv(file.path(repo_root, "data", "public", "batchcorr_set1", "audit", "accession_batch_balance.csv"), stringsAsFactors = FALSE)
write.csv(batch_counts, file.path(figure_data_dir, "batch_composition.csv"), row.names = FALSE, quote = TRUE)
write.csv(accession_balance, file.path(figure_data_dir, "accession_batch_balance.csv"), row.names = FALSE, quote = TRUE)
p_batch <- ggplot(batch_counts, aes(factor(batch, levels = batch_order), n, fill = sample_type_plot)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c(Study = "#56B4E9", QC = "#D55E00")) +
  labs(title = "Batch composition", subtitle = "All batches retain 5–6 references", x = "Batch in acquisition order", y = "Injections", fill = NULL) +
  theme_publication()
p_reps <- ggplot(accession_balance, aes(n_biological_replicates, fill = cross_batch_replicates)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(`TRUE` = "#009E73", `FALSE` = "#9CA3AF")) +
  labs(title = "Accession replicate structure", subtitle = "346 accessions span more than one batch", x = "Retained biological injections per accession", y = "Accessions", fill = "Cross-batch") +
  theme_publication()
ggsave(file.path(figure_dir, "randomization_batch_composition.pdf"), p_batch, width = 9, height = 5.5)
ggsave(file.path(figure_dir, "randomization_accession_replicates.pdf"), p_reps, width = 8, height = 5.5)

# Figure 3: Raw versus WiNN auto PCA with replicate connections.
pca_methods <- intersect(c("Raw", "WINN auto (QC)"), completed_methods)
pca_data <- bind_rows(lapply(pca_methods, function(label) {
  mat <- method_mats[[label]]
  ids <- study_meta$sample_id[study_meta$sample_id %in% colnames(mat)]
  z <- t(log1p(mat[, ids, drop = FALSE]))
  keep <- apply(z, 2, stats::sd) > 0
  pc <- stats::prcomp(z[, keep, drop = FALSE], center = TRUE, scale. = TRUE)
  variance <- pc$sdev^2 / sum(pc$sdev^2)
  md <- study_meta[match(ids, study_meta$sample_id), , drop = FALSE]
  data.frame(
    method = label,
    sample_id = ids,
    accession_id = md$accession_id,
    batch = md$batch,
    PC1 = pc$x[, 1],
    PC2 = pc$x[, 2],
    PC1_variance = variance[1],
    PC2_variance = variance[2],
    stringsAsFactors = FALSE
  )
}))
write.csv(pca_data, file.path(figure_data_dir, "raw_vs_winn_pca.csv"), row.names = FALSE, quote = TRUE)
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = factor(batch, levels = batch_order), group = interaction(method, accession_id))) +
  geom_line(color = "#D1D5DB", alpha = 0.25, linewidth = 0.25) +
  geom_point(alpha = 0.68, size = 1.2) +
  facet_wrap(~method, scales = "free") +
  labs(title = "Raw versus WiNN auto (QC)", subtitle = "Lines connect genuine biological replicates; color denotes supplied analytical batch", color = "Batch", x = "PC1", y = "PC2") +
  theme_publication()
ggsave(file.path(figure_dir, "raw_vs_winn_pca.pdf"), p_pca, width = 12, height = 5.8)

# Figure 4: method comparison.
comparison_source <- method_metrics |>
  filter(status == "completed") |>
  select(method, batch_pc_r2, accession_pc_r2, heldout_qc_cv_mean, run_order_gam_mean) |>
  pivot_longer(-method, names_to = "metric", values_to = "value") |>
  mutate(
    method = factor(method, levels = method_order),
    metric = factor(metric, levels = c("batch_pc_r2", "accession_pc_r2", "heldout_qc_cv_mean", "run_order_gam_mean"), labels = c("Batch variance (weighted-PC R²)", "Accession variance (weighted-PC R²)", "Held-out QC CV (%)", "Residual run-order GAM deviance"))
  )
write.csv(comparison_source, file.path(figure_data_dir, "method_comparison.csv"), row.names = FALSE, quote = TRUE)
p_compare <- ggplot(comparison_source, aes(method, value, fill = method)) +
  geom_col() +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = method_palette, guide = "none") +
  labs(title = "Technical correction and biological preservation", subtitle = "Lower is better for batch variance, QC CV, and run-order deviance; accession variance is a preservation metric", x = NULL, y = NULL) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 32, hjust = 1))
ggsave(file.path(figure_dir, "method_comparison.pdf"), p_compare, width = 13, height = 8)

# Figure 5: association recovery relative to Raw.
association_plot_data <- association_recovery |>
  mutate(method = factor(method, levels = method_order))
write.csv(association_plot_data, file.path(figure_data_dir, "accession_association_change.csv"), row.names = FALSE, quote = TRUE)
p_assoc <- ggplot(association_plot_data, aes(method, change_from_raw, fill = method)) +
  geom_col() +
  geom_hline(yintercept = 0, color = "#111827", linewidth = 0.4) +
  scale_fill_manual(values = method_palette, guide = "none") +
  labs(title = "Change in accession-associated features relative to Raw", subtitle = "Counts use one-way accession F-tests with BH-adjusted p < 0.05; a larger count alone is not evidence of better correction", x = NULL, y = "Associated-feature count change") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 32, hjust = 1))
ggsave(file.path(figure_dir, "accession_association_change.pdf"), p_assoc, width = 11, height = 5.8)

# Figure 6: reproducibly selected representative QC trajectories.
raw_qc_cv <- qc_cv_rows[qc_cv_rows$method == "Raw" & is.finite(qc_cv_rows$cv_percent), , drop = FALSE]
targets <- stats::quantile(raw_qc_cv$cv_percent, probs = c(0.05, 0.25, 0.45, 0.65, 0.85, 0.95), na.rm = TRUE)
selected_features <- unique(vapply(targets, function(target) raw_qc_cv$feature_id[which.min(abs(raw_qc_cv$cv_percent - target))], character(1)))
trend_methods <- intersect(c("Raw", "WINN auto (QC)"), completed_methods)
qc_trend_data <- bind_rows(lapply(trend_methods, function(label) {
  mat <- method_mats[[label]]
  qc_ids <- qc_meta$sample_id[qc_meta$sample_id %in% colnames(mat)]
  grid <- expand.grid(feature_id = selected_features, sample_id = qc_ids, stringsAsFactors = FALSE)
  grid$method <- label
  grid$log_intensity <- mapply(function(f, s) log1p(mat[f, s]), grid$feature_id, grid$sample_id)
  md <- meta[match(grid$sample_id, meta$sample_id), , drop = FALSE]
  grid$batch <- md$batch
  grid$run_order <- md$run_order
  grid$qc_role <- ifelse(grid$sample_id %in% holdout_ids, "Held out", "Training")
  grid
}))
write.csv(qc_trend_data, file.path(figure_data_dir, "representative_qc_trends.csv"), row.names = FALSE, quote = TRUE)
p_trend <- ggplot(qc_trend_data, aes(run_order, log_intensity, color = qc_role, shape = qc_role)) +
  geom_line(aes(group = interaction(method, feature_id, batch)), color = "#D1D5DB", linewidth = 0.3) +
  geom_point(size = 1.4) +
  facet_grid(feature_id ~ method, scales = "free_y") +
  scale_color_manual(values = c(Training = "#0072B2", `Held out` = "#D55E00")) +
  labs(title = "Representative pooled-QC trajectories", subtitle = "Features selected at fixed quantiles of Raw QC CV; held-out QCs were never controls", x = "Global acquisition sequence number", y = "Deposited/imputed log analysis scale", color = NULL, shape = NULL) +
  theme_publication()
ggsave(file.path(figure_dir, "representative_qc_trends.pdf"), p_trend, width = 12, height = 12)

writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
completion <- list(
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  suitability = "SUITABLE WITH LIMITATIONS",
  input_dimensions = c(features = nrow(x), injections = ncol(x)),
  holdout_qcs = length(holdout_ids),
  attempted_methods = method_order,
  completed_methods = completed_methods,
  failed_methods = setdiff(method_order, completed_methods),
  heldout_labels_supplied_to_methods = FALSE
)
jsonlite::write_json(completion, file.path(result_dir, "analysis_complete.json"), pretty = TRUE, auto_unbox = TRUE)
log_line("Analysis complete. Attempted all nine methods; completed ", length(completed_methods), ".")
