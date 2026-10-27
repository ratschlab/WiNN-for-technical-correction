read_feature_matrix_csv <- function(path, id_column = "feature_id") {
  value <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!id_column %in% names(value)) stop("Missing ", id_column, " in ", path)
  ids <- as.character(value[[id_column]])
  value[[id_column]] <- NULL
  out <- as.matrix(value)
  storage.mode(out) <- "double"
  rownames(out) <- ids
  out
}

read_tab_feature_matrix <- function(path) {
  value <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- as.character(value[[1]])
  value[[1]] <- NULL
  out <- as.matrix(value)
  storage.mode(out) <- "double"
  rownames(out) <- ids
  out
}

load_winn_ablation_dataset <- function(repo_root, dataset) {
  dataset <- match.arg(dataset, c(
    "simulation", "mtbls79", "clinical_fiams", "batchcorr_set1",
    "sacurine", "waveica"
  ))

  if (dataset == "simulation") {
    data_dir <- file.path(repo_root, "data", "simulated", "canonical", "SIM01")
    x <- readRDS(file.path(data_dir, "raw_intensity.rds"))
    truth <- readRDS(file.path(data_dir, "clean_ground_truth.rds"))
    meta <- read.csv(file.path(data_dir, "sample_metadata.csv"), check.names = FALSE)
    meta$sample_id <- as.character(meta$sample_id)
    meta$is_qc <- as.logical(meta$is_qc)
    meta$is_study <- !meta$is_qc
    meta$biological_label <- ifelse(meta$is_qc, "QC", "study")
    meta <- meta[match(colnames(x), meta$sample_id), , drop = FALSE]
    x <- x[, meta$sample_id, drop = FALSE]
    truth <- truth[rownames(x), meta$sample_id, drop = FALSE]
    hidden_ids <- as.character(meta$sample_id[as.logical(meta$is_hidden_reference)])
    truth_annotations <- read.csv(file.path(data_dir, "feature_truth.csv"), check.names = FALSE)
    extra <- list(truth = truth, truth_annotations = truth_annotations)
  } else if (dataset == "mtbls79") {
    processed <- file.path(repo_root, "data", "public", "processed")
    meta <- read.csv(
      file.path(processed, "MTBLS79_metadata.csv"), row.names = 1,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    data_df <- read.csv(
      file.path(processed, "MTBLS79_imputed_data.csv"), row.names = 1,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    meta$sample_id <- rownames(meta)
    meta$batch <- as.integer(meta$Batch)
    meta$run_order <- seq_len(nrow(meta))
    meta$within_batch_order <- ave(meta$run_order, meta$batch, FUN = seq_along)
    meta$class <- as.character(meta$Class)
    meta$sample_rep <- as.integer(meta$Sample_Rep)
    meta$is_qc <- meta$class == "QC"
    meta$is_study <- !meta$is_qc
    meta$biological_label <- meta$class
    x <- t(as.matrix(data_df))
    storage.mode(x) <- "double"
    colnames(x) <- meta$sample_id
    holdout <- read.csv(
      file.path(repo_root, "config", "holdouts", "mtbls79.csv"),
      stringsAsFactors = FALSE
    )
    hidden_ids <- as.character(holdout$sample_id)
    extra <- list()
  } else if (dataset == "clinical_fiams") {
    private_path <- file.path(
      repo_root, "data", "private", "clinical_fiams", "prepared_injection_level.rds"
    )
    if (!file.exists(private_path)) {
      stop(
        "The clinical input is restricted. Place the authorized prepared bundle at ",
        private_path, "."
      )
    }
    private <- readRDS(private_path)
    x <- as.matrix(private$x)
    meta <- private$meta[match(colnames(x), private$meta$sample_id), , drop = FALSE]
    meta$within_batch_order <- ave(meta$run_order, meta$batch, FUN = seq_along)
    meta$is_qc <- meta$role != "clinical"
    meta$is_study <- meta$role == "clinical"
    meta$biological_label <- meta$sample
    hidden_ids <- as.character(meta$sample_id[meta$role == "NIST1950"])
    extra <- list(training_reference_ids = meta$sample_id[meta$role == "control"])
  } else if (dataset == "batchcorr_set1") {
    processed <- file.path(repo_root, "data", "public", "batchcorr_set1", "processed")
    meta <- read.csv(
      file.path(processed, "BatchCorr_set1_metadata.csv"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    x <- read_feature_matrix_csv(file.path(processed, "BatchCorr_set1_imputed_data.csv"))
    meta$is_qc <- as.logical(meta$is_qc)
    meta$is_study <- as.logical(meta$is_study_sample)
    meta$run_order <- as.integer(meta$run_order)
    meta$within_batch_order <- as.integer(meta$within_batch_order)
    meta$biological_label <- as.character(meta$accession_id)
    meta <- meta[order(meta$run_order), , drop = FALSE]
    x <- x[, meta$sample_id, drop = FALSE]
    hidden_ids <- read.csv(
      file.path(repo_root, "config", "holdouts", "batchcorr_set1.csv"),
      stringsAsFactors = FALSE
    )$sample_id
    extra <- list()
  } else if (dataset == "sacurine") {
    processed <- file.path(repo_root, "data", "public", "sacurine", "processed")
    x <- readRDS(file.path(processed, "sacurine_intensity_matrix.rds"))
    meta <- read.csv(
      file.path(processed, "sacurine_metadata.csv"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    meta$is_qc <- as.logical(meta$is_qc)
    meta$is_study <- !meta$is_qc
    meta$run_order <- as.integer(meta$run_order)
    meta$within_batch_order <- as.integer(meta$within_batch_order)
    meta$biological_label <- as.character(meta$gender)
    meta <- meta[match(colnames(x), meta$sample_id), , drop = FALSE]
    hidden_ids <- read.csv(
      file.path(repo_root, "config", "holdouts", "sacurine.csv"),
      stringsAsFactors = FALSE
    )$sample_id
    extra <- list()
  } else {
    processed <- file.path(repo_root, "data", "public", "waveica_adenocarcinoma", "processed")
    x <- readRDS(file.path(processed, "waveica_intensity_matrix.rds"))
    meta <- read.csv(
      file.path(processed, "waveica_metadata.csv"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    meta$is_qc <- as.logical(meta$is_qc)
    meta$is_study <- !meta$is_qc
    meta$run_order <- as.integer(meta$run_order)
    meta$within_batch_order <- as.integer(meta$within_batch_order)
    meta$biological_label <- as.character(meta$biological_group)
    meta <- meta[match(colnames(x), meta$sample_id), , drop = FALSE]
    hidden_ids <- read.csv(
      file.path(repo_root, "config", "holdouts", "waveica.csv"),
      stringsAsFactors = FALSE
    )$sample_id
    extra <- list()
  }

  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (!identical(colnames(x), meta$sample_id)) stop(dataset, ": matrix/metadata sample order mismatch")
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) stop(dataset, ": duplicate identifiers")
  if (any(!is.finite(x)) || any(x < 0)) stop(dataset, ": input is not finite and nonnegative")
  if (anyNA(meta$batch) || anyNA(meta$run_order)) stop(dataset, ": missing batch or run order")
  if (!all(hidden_ids %in% colnames(x))) stop(dataset, ": hidden reference IDs are not all present")

  c(list(
    dataset = dataset,
    x = x,
    meta = meta,
    hidden_ids = as.character(hidden_ids)
  ), extra)
}

ablation_variant_table <- function() {
  data.frame(
    variant_id = c(
      "C0_RAW", "C1_OUTLIER", "C2_SELECTIVE_DRIFT", "C3_SELECTIVE_BATCH",
      "C4_FULL_FIXED", "G_SS", "G_AS", "G_SA", "G_AA"
    ),
    variant_label = c(
      "Raw", "Outlier", "+ selective drift", "+ selective batch", "+ PQN",
      "selective/selective", "all/selective", "selective/all", "all/all"
    ),
    analysis_type = c(rep("cumulative", 5L), rep("gate_factorial", 4L)),
    outlier_shrinkage = c(FALSE, rep(TRUE, 8L)),
    drift_gate = c("none", "none", "selective", "selective", "selective", "selective", "all", "selective", "all"),
    batch_gate = c("none", "none", "none", "selective", "selective", "selective", "selective", "all", "all"),
    pqn_mode = c("none", "none", "none", "none", "shrink", rep("shrink", 4L)),
    stringsAsFactors = FALSE
  )
}

capture_ablation_call <- function(fn) {
  warnings <- character()
  messages <- character()
  trace <- ""
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
    error = function(e) {
      trace <<- paste(capture.output(traceback(2)), collapse = "\n")
      structure(list(message = conditionMessage(e)), class = "ablation_error")
    }
  )
  list(
    value = value,
    elapsed_sec = proc.time()[["elapsed"]] - started,
    warnings = warnings,
    messages = messages,
    traceback = trace,
    error = if (inherits(value, "ablation_error")) value$message else ""
  )
}

derive_forced_batch_variant <- function(source_result, batch, control_samples = NULL) {
  batch_started <- proc.time()[["elapsed"]]
  batch_result <- winn::anova_batch_correction(
    source_result$intermediates$post_drift,
    batch = batch,
    fdr_threshold = 0.05,
    gate = "all",
    return_diagnostics = TRUE
  )
  batch_runtime <- proc.time()[["elapsed"]] - batch_started
  pqn_started <- proc.time()[["elapsed"]]
  pqn_result <- winn::normalize_by_dilution_factor(
    batch_result$data,
    processing = "shrink",
    control_samples = control_samples,
    return_diagnostics = TRUE
  )
  pqn_runtime <- proc.time()[["elapsed"]] - pqn_started
  final <- pqn_result$data
  dimnames(final) <- dimnames(source_result$intermediates$raw)
  stage_runtime <- source_result$stage_runtime_sec
  stage_runtime[["batch"]] <- batch_runtime
  stage_runtime[["pqn"]] <- pqn_runtime
  intermediates <- source_result$intermediates
  intermediates$post_batch <- batch_result$data
  intermediates$final <- final
  diagnostics <- source_result$diagnostics
  diagnostics$batch <- batch_result$diagnostics
  diagnostics$pqn <- pqn_result$diagnostics
  diagnostics$matrix_state <- do.call(rbind, lapply(names(intermediates), function(stage) {
    value <- intermediates[[stage]]
    data.frame(
      stage = stage, n_features = nrow(value), n_samples = ncol(value),
      na_count = sum(is.na(value)), nan_count = sum(is.nan(value)),
      infinite_count = sum(is.infinite(value)), minimum = min(value), maximum = max(value),
      stringsAsFactors = FALSE
    )
  }))
  list(
    data = final,
    intermediates = intermediates,
    diagnostics = diagnostics,
    stage_runtime_sec = stage_runtime,
    total_runtime_sec = sum(stage_runtime),
    configuration = within(source_result$configuration, batch_gate <- "all")
  )
}

empty_raw_result <- function(x) {
  stage_runtime <- stats::setNames(rep(0, 4L), c("outlier", "drift", "batch", "pqn"))
  list(
    data = x,
    intermediates = list(raw = x, post_outlier = x, post_drift = x, post_batch = x, final = x),
    diagnostics = list(outlier = list(
      entries_changed = 0L, proportion_entries_changed = 0,
      features_changed = 0L, samples_changed = 0L,
      median_absolute_log1p_change = 0, p90_absolute_log1p_change = 0
    ), drift = data.frame(), batch = data.frame(), pqn = data.frame(), matrix_state = data.frame()),
    stage_runtime_sec = stage_runtime,
    total_runtime_sec = 0,
    configuration = list(
      parameters = "fixed", supplied_batch = TRUE,
      control_samples_supplied = FALSE, use_outlier_shrinkage = FALSE,
      drift_gate = "none", batch_gate = "none", pqn_mode = "none",
      fdr_threshold = 0.05, test = "Ljung-Box", lag = "adaptive",
      spline_method = "conservative", remove_batch_effects = "anova",
      scale_by_batch = FALSE
    )
  )
}

subset_cumulative_result <- function(base_result, stage) {
  stage <- match.arg(stage, c("outlier", "drift", "batch", "final"))
  matrix_name <- switch(stage,
    outlier = "post_outlier", drift = "post_drift", batch = "post_batch", final = "final"
  )
  stage_index <- match(stage, c("outlier", "drift", "batch", "final"))
  keep_runtime <- c("outlier", "drift", "batch", "pqn")[seq_len(stage_index)]
  runtime <- base_result$stage_runtime_sec
  runtime[!names(runtime) %in% keep_runtime] <- 0
  diagnostics <- base_result$diagnostics
  if (stage_index < 2L) diagnostics$drift <- data.frame()
  if (stage_index < 3L) diagnostics$batch <- data.frame()
  if (stage_index < 4L) diagnostics$pqn <- data.frame()
  configuration <- base_result$configuration
  if (stage_index < 2L) configuration$drift_gate <- "none"
  if (stage_index < 3L) configuration$batch_gate <- "none"
  if (stage_index < 4L) configuration$pqn_mode <- "none"
  out <- base_result
  out$data <- base_result$intermediates[[matrix_name]]
  out$intermediates$final <- out$data
  out$stage_runtime_sec <- runtime
  out$total_runtime_sec <- sum(runtime)
  out$diagnostics <- diagnostics
  out$configuration <- configuration
  out
}

safe_cor_ablation <- function(a, b, method = "pearson", min_n = 3L) {
  keep <- is.finite(a) & is.finite(b)
  if (sum(keep) < min_n || stats::sd(a[keep]) <= 0 || stats::sd(b[keep]) <= 0) return(NA_real_)
  suppressWarnings(stats::cor(a[keep], b[keep], method = method))
}

icc_a1_ablation <- function(ratings) {
  ratings <- as.matrix(ratings)
  ratings <- ratings[complete.cases(ratings), , drop = FALSE]
  n <- nrow(ratings); k <- ncol(ratings)
  if (n < 2L || k < 2L) return(NA_real_)
  grand <- mean(ratings); row_mean <- rowMeans(ratings); col_mean <- colMeans(ratings)
  ss_row <- k * sum((row_mean - grand)^2)
  ss_col <- n * sum((col_mean - grand)^2)
  ss_error <- sum((ratings - grand)^2) - ss_row - ss_col
  ms_row <- ss_row / (n - 1L); ms_col <- ss_col / (k - 1L)
  ms_error <- ss_error / ((n - 1L) * (k - 1L))
  denominator <- ms_row + (k - 1L) * ms_error + k * (ms_col - ms_error) / n
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  (ms_row - ms_error) / denominator
}

weighted_pc_r2_ablation <- function(mat, target, target_type = NULL, n_pcs = 5L) {
  if (!exists("weighted_pc_r2_explicit", mode = "function", inherits = TRUE)) {
    stop("Source scripts/weighted_pc_r2.R before using weighted_pc_r2_ablation().")
  }
  weighted_pc_r2_explicit(
    mat = mat,
    target = target,
    target_type = target_type,
    n_pcs = n_pcs
  )
}

qc_cv_ablation <- function(mat, ids) {
  ids <- intersect(ids, colnames(mat))
  if (length(ids) < 2L) return(stats::setNames(rep(NA_real_, nrow(mat)), rownames(mat)))
  apply(mat[, ids, drop = FALSE], 1, function(v) {
    mean_value <- mean(v, na.rm = TRUE)
    if (!is.finite(mean_value) || abs(mean_value) <= .Machine$double.eps) return(NA_real_)
    100 * stats::sd(v, na.rm = TRUE) / abs(mean_value)
  })
}

hidden_profile_correlations <- function(mat, ids) {
  ids <- intersect(ids, colnames(mat))
  if (length(ids) < 2L) return(data.frame(method = character(), value = numeric()))
  z <- log1p(mat[, ids, drop = FALSE])
  bind_one <- function(method) {
    cm <- suppressWarnings(stats::cor(z, use = "pairwise.complete.obs", method = method))
    data.frame(correlation_method = method, value = cm[upper.tri(cm)], stringsAsFactors = FALSE)
  }
  rbind(bind_one("pearson"), bind_one("spearman"))
}

correction_magnitude_summary <- function(current, previous, raw, tolerance = 1e-10) {
  current_log <- log1p(current)
  previous_log <- log1p(previous)
  raw_log <- log1p(raw)
  delta_previous <- abs(current_log - previous_log)
  delta_raw <- abs(current_log - raw_log)
  sample_cor <- vapply(seq_len(ncol(current)), function(i) {
    safe_cor_ablation(current_log[, i], raw_log[, i])
  }, numeric(1))
  feature_cor <- vapply(seq_len(nrow(current)), function(i) {
    safe_cor_ablation(current_log[i, ], raw_log[i, ])
  }, numeric(1))
  data.frame(
    median_abs_log1p_change_preceding = stats::median(delta_previous, na.rm = TRUE),
    p90_abs_log1p_change_preceding = as.numeric(stats::quantile(delta_previous, 0.9, na.rm = TRUE, names = FALSE)),
    median_abs_log1p_change_raw = stats::median(delta_raw, na.rm = TRUE),
    p90_abs_log1p_change_raw = as.numeric(stats::quantile(delta_raw, 0.9, na.rm = TRUE, names = FALSE)),
    proportion_entries_changed = mean(delta_raw > tolerance, na.rm = TRUE),
    median_sample_profile_pearson_raw = stats::median(sample_cor, na.rm = TRUE),
    median_feature_profile_pearson_raw = stats::median(feature_cor, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

effect_concordance_ablation <- function(a, b) {
  keep <- is.finite(a) & is.finite(b)
  if (sum(keep) < 3L) return(c(pearson = NA, spearman = NA, direction = NA))
  c(
    pearson = stats::cor(a[keep], b[keep], method = "pearson"),
    spearman = stats::cor(a[keep], b[keep], method = "spearman"),
    direction = mean(sign(a[keep]) == sign(b[keep]))
  )
}

mtbls79_replicate_metrics <- function(mat, meta) {
  ids <- meta$sample_id[meta$is_study]
  md <- meta[match(ids, meta$sample_id), , drop = FALSE]
  z <- log1p(mat[, ids, drop = FALSE])
  groups <- split(seq_along(ids), md$sample_rep)
  groups <- groups[lengths(groups) >= 2L]
  rows <- lapply(names(groups), function(group_id) {
    index <- groups[[group_id]]
    cm <- suppressWarnings(stats::cor(z[, index, drop = FALSE], use = "pairwise.complete.obs"))
    data.frame(
      sample_rep = group_id,
      pearson = stats::median(cm[lower.tri(cm)], na.rm = TRUE),
      icc_a1 = icc_a1_ablation(z[, index, drop = FALSE]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mtbls79_feature_cross_batch_icc <- function(mat, meta) {
  ids <- meta$sample_id[meta$is_study]
  md <- meta[match(ids, meta$sample_id), , drop = FALSE]
  z <- log1p(mat[, ids, drop = FALSE])
  batches <- sort(unique(md$batch))
  pairs <- utils::combn(batches, 2L, simplify = FALSE)
  pair_index <- lapply(pairs, function(pair) {
    i1 <- which(md$batch == pair[1]); i2 <- which(md$batch == pair[2])
    shared <- intersect(md$sample_rep[i1], md$sample_rep[i2])
    if (length(shared) < 3L) return(NULL)
    list(
      i1 = i1[match(shared, md$sample_rep[i1])],
      i2 = i2[match(shared, md$sample_rep[i2])]
    )
  })
  pair_index <- Filter(Negate(is.null), pair_index)
  vapply(seq_len(nrow(z)), function(feature_index) {
    values <- vapply(pair_index, function(pair) {
      icc_a1_ablation(cbind(z[feature_index, pair$i1], z[feature_index, pair$i2]))
    }, numeric(1))
    stats::median(values, na.rm = TRUE)
  }, numeric(1))
}

mtbls79_associations <- function(mat, meta, variant_id) {
  ids <- meta$sample_id[meta$is_study]
  md <- meta[match(ids, meta$sample_id), , drop = FALSE]
  z <- log1p(mat[, ids, drop = FALSE])
  groups <- split(seq_along(ids), md$sample_rep)
  collapsed <- vapply(groups, function(index) rowMeans(z[, index, drop = FALSE]), numeric(nrow(z)))
  classes <- vapply(groups, function(index) names(which.max(table(md$class[index]))), character(1))
  p_value <- apply(collapsed, 1, function(v) {
    tryCatch(stats::t.test(v[classes == "C"], v[classes == "S"])$p.value, error = function(e) NA_real_)
  })
  effect <- rowMeans(collapsed[, classes == "C", drop = FALSE]) -
    rowMeans(collapsed[, classes == "S", drop = FALSE])
  p_adj <- stats::p.adjust(p_value, method = "BH")
  data.frame(
    variant_id = variant_id, feature_id = rownames(mat), p_value = p_value,
    p_adj = p_adj, effect = effect,
    significant = is.finite(p_adj) & p_adj < 0.05 & abs(effect) >= 0.2,
    stringsAsFactors = FALSE
  )
}

batchcorr_replicate_metrics <- function(mat, meta) {
  ids <- meta$sample_id[meta$is_study]
  md <- meta[match(ids, meta$sample_id), , drop = FALSE]
  z <- log1p(mat[, ids, drop = FALSE])
  groups <- split(seq_along(ids), md$accession_id)
  groups <- groups[lengths(groups) >= 2L]
  do.call(rbind, lapply(names(groups), function(accession) {
    index <- groups[[accession]]
    pairs <- utils::combn(index, 2L, simplify = FALSE)
    pair_values <- do.call(rbind, lapply(pairs, function(pair) {
      data.frame(
        pearson = safe_cor_ablation(z[, pair[1]], z[, pair[2]]),
        icc_a1 = icc_a1_ablation(cbind(z[, pair[1]], z[, pair[2]]))
      )
    }))
    data.frame(
      accession_id = accession, n_replicates = length(index),
      pearson = stats::median(pair_values$pearson, na.rm = TRUE),
      icc_a1 = stats::median(pair_values$icc_a1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

batchcorr_feature_repeatability <- function(mat, meta) {
  ids <- meta$sample_id[meta$is_study]
  md <- meta[match(ids, meta$sample_id), , drop = FALSE]
  z <- log1p(mat[, ids, drop = FALSE])
  groups <- split(seq_along(ids), md$accession_id)
  groups <- groups[lengths(groups) >= 2L]
  vapply(seq_len(nrow(z)), function(feature_index) {
    y <- z[feature_index, ]
    means <- vapply(groups, function(index) mean(y[index]), numeric(1))
    between <- stats::var(means)
    within_numerator <- sum(vapply(groups, function(index) {
      stats::var(y[index]) * (length(index) - 1L)
    }, numeric(1)))
    within_df <- sum(lengths(groups) - 1L)
    within <- if (within_df > 0L) within_numerator / within_df else NA_real_
    if (!is.finite(between) || !is.finite(within) || between + within <= 0) NA_real_ else between / (between + within)
  }, numeric(1))
}

batchcorr_associations <- function(mat, meta, variant_id) {
  ids <- meta$sample_id[meta$is_study]
  md <- meta[match(ids, meta$sample_id), , drop = FALSE]
  group <- factor(md$accession_id)
  z <- log1p(mat[, ids, drop = FALSE])
  group_index <- split(seq_along(group), group)
  n <- length(group); k <- nlevels(group)
  rows <- lapply(seq_len(nrow(z)), function(i) {
    y <- z[i, ]; grand <- mean(y)
    means <- vapply(group_index, function(index) mean(y[index]), numeric(1))
    ss_between <- sum(lengths(group_index) * (means - grand)^2)
    ss_total <- sum((y - grand)^2); ss_within <- max(0, ss_total - ss_between)
    f_stat <- if (k > 1L && n > k && ss_within > 0) {
      (ss_between / (k - 1L)) / (ss_within / (n - k))
    } else NA_real_
    data.frame(
      variant_id = variant_id, feature_id = rownames(z)[i],
      f_statistic = f_stat,
      p_value = if (is.finite(f_stat)) stats::pf(f_stat, k - 1L, n - k, lower.tail = FALSE) else NA_real_,
      eta_squared = if (ss_total > 0) ss_between / ss_total else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out$significant <- is.finite(out$p_adj) & out$p_adj < 0.05
  out
}

limma_fit_ablation <- function(mat, meta, design, coefficients, variant_id, model) {
  z <- log1p(mat[, meta$sample_id, drop = FALSE])
  fit <- limma::eBayes(limma::lmFit(z, design), robust = TRUE)
  do.call(rbind, lapply(names(coefficients), function(variable) {
    index <- match(coefficients[[variable]], colnames(design))
    if (is.na(index)) return(NULL)
    coefficient <- fit$coefficients[, index]
    t_value <- fit$t[, index]
    outcome_sd <- apply(z, 1, stats::sd)
    p_value <- fit$p.value[, index]
    p_adj <- stats::p.adjust(p_value, method = "BH")
    data.frame(
      variant_id = variant_id, model = model, variable = variable,
      feature_id = rownames(z), coefficient = coefficient,
      standardized_effect = coefficient / outcome_sd,
      t_statistic = t_value, p_value = p_value, p_adj = p_adj,
      partial_r_squared = t_value^2 / (t_value^2 + fit$df.total),
      significant = is.finite(p_adj) & p_adj < 0.05,
      stringsAsFactors = FALSE
    )
  }))
}

sacurine_biology <- function(mat, meta, variant_id) {
  md <- meta[meta$is_study, , drop = FALSE]
  md$age_z <- as.numeric(scale(md$age))
  md$bmi_z <- as.numeric(scale(md$bmi))
  md$gender <- factor(md$gender)
  md$sampling_factor <- factor(md$sampling)
  biological <- sweep(mat[, md$sample_id, drop = FALSE], 2, md$osmolality, "/")
  design <- stats::model.matrix(~ age_z + bmi_z + gender + sampling_factor, data = md)
  gender_name <- grep("^gender", colnames(design), value = TRUE)[1]
  coefficients <- c(age = "age_z", bmi = "bmi_z", gender = gender_name)
  primary <- limma_fit_ablation(
    biological, md, design, coefficients, variant_id,
    "primary_osmolality_normalized"
  )
  by_batch <- do.call(rbind, lapply(unique(md$batch), function(batch_id) {
    d <- md[md$batch == batch_id, , drop = FALSE]
    batch_design <- stats::model.matrix(~ age_z + bmi_z + gender + sampling_factor, data = d)
    batch_gender <- grep("^gender", colnames(batch_design), value = TRUE)[1]
    limma_fit_ablation(
      biological, d, batch_design,
      c(age = "age_z", bmi = "bmi_z", gender = batch_gender),
      variant_id, paste0("batch_", batch_id)
    )
  }))
  pairs <- utils::combn(unique(md$batch), 2L, simplify = FALSE)
  concordance <- do.call(rbind, lapply(c("age", "bmi", "gender"), function(variable) {
    do.call(rbind, lapply(pairs, function(pair) {
      a <- by_batch[by_batch$variable == variable & by_batch$model == paste0("batch_", pair[1]), ]
      b <- by_batch[by_batch$variable == variable & by_batch$model == paste0("batch_", pair[2]), ]
      b <- b[match(a$feature_id, b$feature_id), , drop = FALSE]
      value <- effect_concordance_ablation(a$standardized_effect, b$standardized_effect)
      data.frame(
        variant_id = variant_id, variable = variable,
        batch_1 = pair[1], batch_2 = pair[2],
        pearson = value[["pearson"]], spearman = value[["spearman"]],
        direction_concordance = value[["direction"]], stringsAsFactors = FALSE
      )
    }))
  }))
  pc <- data.frame(
    variable = c("age", "bmi", "gender"),
    value = c(
      weighted_pc_r2_ablation(biological, md$age, target_type = "continuous"),
      weighted_pc_r2_ablation(biological, md$bmi, target_type = "continuous"),
      weighted_pc_r2_ablation(biological, md$gender, target_type = "categorical")
    )
  )
  list(primary = primary, by_batch = by_batch, concordance = concordance, pc = pc)
}

waveica_biology <- function(mat, meta, variant_id) {
  md <- meta[meta$is_study, , drop = FALSE]
  md$biological_group <- factor(md$biological_group, levels = c("group_0", "group_1"))
  fit_group <- function(model_md, model) {
    design <- stats::model.matrix(~ biological_group, data = model_md)
    coefficient <- grep("^biological_group", colnames(design), value = TRUE)[1]
    tab <- limma_fit_ablation(mat, model_md, design, c(group = coefficient), variant_id, model)
    z <- log1p(mat[, model_md$sample_id, drop = FALSE])
    group <- model_md$biological_group
    n0 <- sum(group == "group_0"); n1 <- sum(group == "group_1")
    sd0 <- apply(z[, group == "group_0", drop = FALSE], 1, stats::sd)
    sd1 <- apply(z[, group == "group_1", drop = FALSE], 1, stats::sd)
    pooled <- sqrt(((n0 - 1L) * sd0^2 + (n1 - 1L) * sd1^2) / (n0 + n1 - 2L))
    tab$standardized_effect <- tab$coefficient / pooled[match(tab$feature_id, rownames(z))]
    tab
  }
  primary <- fit_group(md, "primary")
  by_batch <- do.call(rbind, lapply(unique(md$batch), function(batch_id) {
    fit_group(md[md$batch == batch_id, , drop = FALSE], paste0("batch_", batch_id))
  }))
  pairs <- utils::combn(unique(md$batch), 2L, simplify = FALSE)
  concordance <- do.call(rbind, lapply(pairs, function(pair) {
    a <- by_batch[by_batch$model == paste0("batch_", pair[1]), ]
    b <- by_batch[by_batch$model == paste0("batch_", pair[2]), ]
    b <- b[match(a$feature_id, b$feature_id), , drop = FALSE]
    value <- effect_concordance_ablation(a$standardized_effect, b$standardized_effect)
    data.frame(
      variant_id = variant_id, batch_1 = pair[1], batch_2 = pair[2],
      pearson = value[["pearson"]], spearman = value[["spearman"]],
      direction_concordance = value[["direction"]], stringsAsFactors = FALSE
    )
  }))
  wide <- reshape(
    by_batch[, c("feature_id", "model", "standardized_effect")],
    idvar = "feature_id", timevar = "model", direction = "wide"
  )
  effect_columns <- grep("^standardized_effect", names(wide), value = TRUE)
  consistent <- apply(wide[, effect_columns, drop = FALSE], 1, function(v) {
    all(is.finite(v)) && (all(v > 0) || all(v < 0))
  })
  list(
    primary = primary, by_batch = by_batch, concordance = concordance,
    group_pc_r2 = weighted_pc_r2_ablation(
      mat[, md$sample_id, drop = FALSE], md$biological_group,
      target_type = "categorical"
    ),
    direction_consistency = mean(consistent[primary$significant], na.rm = TRUE)
  )
}

metric_row <- function(metric, value, direction, denominator, units = "", notes = "") {
  data.frame(
    metric = metric, value = as.numeric(value), metric_direction = direction,
    denominator = as.character(denominator), units = units, notes = notes,
    stringsAsFactors = FALSE
  )
}
