select_hidden_qcs <- function(metadata, per_batch, seed, min_train_per_batch) {
  set.seed(seed)
  selected <- character()
  audit_rows <- list()
  for (batch_id in unique(as.character(metadata$batch))) {
    q <- metadata[metadata$is_qc & metadata$batch == batch_id, , drop = FALSE]
    q <- q[order(q$within_batch_order), , drop = FALSE]
    if (nrow(q) < min_train_per_batch + 1L) stop("Insufficient QCs in batch ", batch_id)
    protected <- c(q$sample_id[1], q$sample_id[nrow(q)])
    candidates <- setdiff(q$sample_id, protected)
    feasible <- min(as.integer(per_batch), length(candidates), nrow(q) - min_train_per_batch)
    if (feasible < 1L) stop("No feasible interior QC holdout in batch ", batch_id)
    chosen <- sample(candidates, feasible, replace = FALSE)
    selected <- c(selected, chosen)
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      batch = batch_id,
      total_qc = nrow(q),
      target_hidden = per_batch,
      selected_hidden = feasible,
      training_qc = nrow(q) - feasible,
      protected_first_qc = protected[1],
      protected_last_qc = protected[2],
      stringsAsFactors = FALSE
    )
  }
  selected <- metadata$sample_id[metadata$sample_id %in% selected]
  training <- metadata$sample_id[metadata$is_qc & !metadata$sample_id %in% selected]
  list(hidden = selected, training = training, audit = dplyr::bind_rows(audit_rows))
}

capture_call_full <- function(fn) {
  warnings <- character()
  messages <- character()
  start <- proc.time()[["elapsed"]]
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
  list(
    value = value,
    elapsed_sec = proc.time()[["elapsed"]] - start,
    warnings = unique(warnings),
    messages = unique(messages),
    error = if (inherits(value, "error")) conditionMessage(value) else ""
  )
}

training_qc_score <- function(mat, training_ids) {
  ids <- training_ids[training_ids %in% colnames(mat)]
  if (length(ids) < 2L) return(data.frame(mean_qc_correlation = NA_real_, mean_qc_cv = NA_real_, qc_score = NA_real_))
  q <- mat[, ids, drop = FALSE]
  cvs <- apply(q, 1, function(v) {
    mu <- mean(v, na.rm = TRUE)
    if (!is.finite(mu) || abs(mu) <= .Machine$double.eps) NA_real_ else stats::sd(v, na.rm = TRUE) / abs(mu)
  })
  correlations <- suppressWarnings(stats::cor(q, use = "pairwise.complete.obs"))
  correlation_values <- correlations[lower.tri(correlations)]
  mean_correlation <- mean(correlation_values, na.rm = TRUE)
  mean_cv <- mean(cvs, na.rm = TRUE)
  data.frame(
    mean_qc_correlation = mean_correlation,
    mean_qc_cv = mean_cv,
    qc_score = mean_correlation - 0.5 * mean_cv,
    stringsAsFactors = FALSE
  )
}

run_qc_rlsc_id_safe <- function(
    x_mat, control_ids, meta_df, span = 0.75, degree = 1L,
    shift_batches = TRUE, min_qcs_per_batch = 4L) {
  min_qcs_per_batch <- as.integer(min_qcs_per_batch)
  if (!is.finite(min_qcs_per_batch) || min_qcs_per_batch < 2L) {
    stop("min_qcs_per_batch must be an integer of at least two.")
  }
  md <- align_meta_to_mat(x_mat, meta_df)
  md <- md[order(md$run_order), , drop = FALSE]
  x_ord <- x_mat[, md$sample_id, drop = FALSE]
  dat_log <- as.data.frame(t(log1p(x_ord)), check.names = FALSE)
  qc_labels <- ifelse(md$sample_id %in% control_ids, "QC", "Sample")
  batch_levels <- unique(as.character(md$batch))
  batch_factor <- factor(md$batch, levels = batch_levels)
  counts <- table(batch_factor[md$sample_id %in% control_ids])
  if (length(counts) != length(batch_levels) || any(counts < min_qcs_per_batch)) {
    stop(
      "QC-RLSC requires at least ", min_qcs_per_batch,
      " training QCs in every batch."
    )
  }

  corrected_log <- matrix(NA_real_, nrow = nrow(dat_log), ncol = ncol(dat_log), dimnames = dimnames(dat_log))
  for (batch_id in batch_levels) {
    idx <- which(md$batch == batch_id)
    batch_corrected <- as.matrix(qcrlscR::qc.rlsc(
      x = dat_log[idx, , drop = FALSE],
      y = qc_labels[idx],
      method = "subtract",
      opti = FALSE,
      span = span,
      degree = degree
    ))
    if (!identical(rownames(batch_corrected), md$sample_id[idx])) stop("QC-RLSC changed sample identifiers in batch ", batch_id)
    corrected_log[idx, ] <- batch_corrected
  }
  if (any(!is.finite(corrected_log))) stop("QC-RLSC returned nonfinite within-batch values.")
  if (isTRUE(shift_batches)) {
    corrected_log <- as.matrix(qcrlscR::batch.shift(as.data.frame(corrected_log, check.names = FALSE), batch_factor, overall_average = TRUE))
  }
  if (!identical(rownames(corrected_log), md$sample_id)) stop("QC-RLSC changed sample identifiers during batch shifting.")
  corrected <- t(expm1(corrected_log))
  rownames(corrected) <- rownames(x_ord)
  colnames(corrected) <- md$sample_id
  corrected <- pmax(corrected, 0)
  corrected[, colnames(x_mat), drop = FALSE]
}

audit_qcrlsc_wrapper <- function(x_mat, control_ids, meta_df, audit_dir) {
  subset_features <- rownames(x_mat)[seq_len(min(20L, nrow(x_mat)))]
  xx <- x_mat[subset_features, , drop = FALSE]
  md <- align_meta_to_mat(xx, meta_df)
  md <- md[order(md$run_order), , drop = FALSE]
  dat <- as.data.frame(t(xx[, md$sample_id, drop = FALSE]), check.names = FALSE)
  labels <- ifelse(md$sample_id %in% control_ids, "QC", "Sample")
  run_one <- function(span) {
    qcrlscR::qc.rlsc.wrap(
      dat = dat, cls.qc = labels,
      cls.bl = factor(md$batch, levels = unique(md$batch)),
      method = "divide", intra = TRUE, opti = FALSE,
      log10 = FALSE, outl = FALSE, shift = TRUE, span = span
    )
  }
  a <- capture_call_full(function() run_one(0.30))
  b <- capture_call_full(function() run_one(0.90))
  a_ok <- !inherits(a$value, "error")
  b_ok <- !inherits(b$value, "error")
  sample_order_ok <- a_ok && identical(rownames(as.matrix(a$value)), md$sample_id)
  span_outputs_identical <- a_ok && b_ok && isTRUE(all.equal(as.matrix(a$value), as.matrix(b$value), tolerance = 1e-12, check.attributes = FALSE))

  direct_batch <- unique(md$batch)[1]
  idx <- which(md$batch == direct_batch)
  direct_input <- as.data.frame(t(log1p(xx[, md$sample_id[idx], drop = FALSE])), check.names = FALSE)
  direct_labels <- labels[idx]
  direct <- capture_call_full(function() qcrlscR::qc.rlsc(direct_input, direct_labels, method = "subtract", opti = FALSE, span = 0.75, degree = 1))
  direct_id_safe <- !inherits(direct$value, "error") && identical(rownames(as.matrix(direct$value)), md$sample_id[idx])

  audit <- data.frame(
    check = c("canonical_wrapper_completed_span_0.30", "canonical_wrapper_completed_span_0.90", "canonical_wrapper_preserved_sample_order", "canonical_wrapper_respects_requested_span", "direct_single_batch_qc_rlsc_preserved_ids"),
    passed = c(a_ok, b_ok, sample_order_ok, !span_outputs_identical, direct_id_safe),
    observed = c(a$error, b$error, as.character(sample_order_ok), as.character(span_outputs_identical), as.character(direct_id_safe)),
    interpretation = c(
      "smoke test", "smoke test", "required for matrix alignment",
      "failed because outputs_identical=TRUE proves qc.rlsc.wrap ignored the supplied span",
      "validates the ID-safe per-batch primitive used by the benchmark"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(audit, file.path(audit_dir, "qc_rlsc_wrapper_audit.csv"), row.names = FALSE, quote = TRUE)
  list(audit = audit, custom_required = !sample_order_ok || span_outputs_identical || !direct_id_safe)
}

parse_winn_auto_messages <- function(messages) {
  joined <- paste(messages, collapse = "\n")
  extract <- function(pattern, default = NA_character_) {
    match_value <- regexec(pattern, joined, perl = TRUE)
    found <- regmatches(joined, match_value)[[1]]
    if (length(found) >= 2L) trimws(found[2]) else default
  }
  list(
    test = extract("Autocorrelation test: ([^ ]+)", "Ljung-Box"),
    acorr_fdr = as.numeric(extract("Autocorrelation test: [^\\n]+\\(FDR: ([0-9.]+)\\)", "0.05")),
    spline_method = extract("Spline method: ([^\\n]+)", "conservative"),
    anova_fdr = as.numeric(extract("Batch correction FDR: ([0-9.]+)", "0.05")),
    normalization = extract("Normalization: ([^\\n]+)", "shrink"),
    scale_by_batch = identical(tolower(extract("Scale by batch: ([^\\n]+)", "false")), "true"),
    pelt_penalty = suppressWarnings(as.numeric(extract("fkPELT penalty: ([0-9.eE+-]+)", NA_character_))),
    mean_cv = suppressWarnings(as.numeric(extract("Quality metrics - CV: ([0-9.eE+-]+)", NA_character_))),
    mean_correlation = suppressWarnings(as.numeric(extract("Correlation: ([0-9.eE+-]+)", NA_character_))),
    messages = joined
  )
}

run_winn_auto_tuning_subset <- function(x_tune, control_ids, meta_df, auto_batch) {
  captured <- capture_call_full(function() run_winn_with_controls(
    x_tune, control_ids = control_ids, meta_df = meta_df,
    auto_batch = auto_batch, parameters = "auto"
  ))
  if (inherits(captured$value, "error")) stop(captured$error)
  params <- parse_winn_auto_messages(captured$messages)
  if (!is.finite(params$acorr_fdr) || !is.finite(params$anova_fdr)) stop("Could not parse WiNN automatic parameters.")
  list(params = params, captured = captured)
}

winn_drift_test_table <- function(data, batch, run_order, test, fdr_threshold) {
  model_df <- get(".autocorrelation_model_df", envir = asNamespace("winn"))(test)
  resolve_lag <- get(".resolve_autocorrelation_lag", envir = asNamespace("winn"))
  rows <- list()
  for (batch_id in unique(batch)) {
    idx <- which(batch == batch_id)
    segment <- data[, idx, drop = FALSE]
    segment_order <- run_order[idx]
    lag_value <- resolve_lag(NULL, length(idx), model_df)
    p_values <- apply(segment, 1, function(v) {
      z <- log1p(v)
      if (test == "Ljung-Box") {
        if (is.na(lag_value)) return(NA_real_)
        tryCatch(stats::Box.test(z, lag = lag_value, type = "Ljung-Box", fitdf = model_df)$p.value, error = function(e) NA_real_)
      } else {
        tryCatch(lmtest::dwtest(z ~ segment_order)$p.value, error = function(e) NA_real_)
      }
    })
    adjusted <- stats::p.adjust(p_values, method = "fdr")
    rows[[length(rows) + 1L]] <- data.frame(
      feature_id = rownames(data), batch = as.character(batch_id), test = test,
      lag = if (test == "Ljung-Box") lag_value else NA_integer_,
      p_value = p_values, p_adj = adjusted,
      selected_for_drift = is.finite(adjusted) & adjusted < fdr_threshold,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

winn_batch_test_table <- function(data, batch, fdr_threshold) {
  batch <- factor(batch)
  z <- log1p(data)
  if (nlevels(batch) < 2L) {
    return(data.frame(feature_id = rownames(data), p_value = NA_real_, p_adj = NA_real_, selected_for_batch = FALSE))
  }
  p_values <- vapply(seq_len(nrow(z)), function(i) {
    fit <- stats::aov(z[i, ] ~ batch)
    summary(fit)[[1]]$`Pr(>F)`[1]
  }, numeric(1))
  adjusted <- stats::p.adjust(p_values, method = "fdr")
  data.frame(feature_id = rownames(data), p_value = p_values, p_adj = adjusted, selected_for_batch = is.finite(adjusted) & adjusted < fdr_threshold, stringsAsFactors = FALSE)
}

run_winn_instrumented <- function(x_mat, meta_df, control_ids, auto_batch, params, mode_label) {
  md <- align_meta_to_mat(x_mat, meta_df)
  control_idx <- if (is.null(control_ids)) NULL else which(md$sample_id %in% control_ids)
  adjust_outliers <- get("adjust_outliers_mad", envir = asNamespace("winn"))
  detect_batch <- get(".auto_detect_batch", envir = asNamespace("winn"))
  norm_data <- adjust_outliers(as.matrix(x_mat))
  batch_vec <- if (auto_batch) detect_batch(norm_data, pelt_penalty = params$pelt_penalty) else md$batch
  drift_table <- winn_drift_test_table(norm_data, batch_vec, md$run_order, params$test, params$acorr_fdr)
  drift_corrected <- winn::autocorrelation_correct(norm_data, run_order = md$run_order, batch = batch_vec, lag = NULL, test = params$test, detrend = "mean", fdr_threshold = params$acorr_fdr, spline_method = params$spline_method)
  batch_table <- winn_batch_test_table(drift_corrected, batch_vec, params$anova_fdr)
  batch_corrected <- winn::anova_batch_correction(drift_corrected, batch_vec, fdr_threshold = params$anova_fdr)
  normalized <- if (identical(params$normalization, "none")) batch_corrected else winn::normalize_by_dilution_factor(batch_corrected, processing = params$normalization, control_samples = control_idx)
  final <- if (isTRUE(params$scale_by_batch)) winn::scale_by_batch(normalized, batch_vec) else normalized
  rownames(final) <- rownames(x_mat)
  colnames(final) <- colnames(x_mat)

  drift_features <- unique(drift_table$feature_id[drift_table$selected_for_drift])
  batch_features <- batch_table$feature_id[batch_table$selected_for_batch]
  pqn_changed <- colSums(abs(log1p(pmax(final, 0)) - log1p(pmax(batch_corrected, 0))) > 1e-10) > 0
  by_feature <- merge(
    aggregate(selected_for_drift ~ feature_id, drift_table, any),
    batch_table[, c("feature_id", "selected_for_batch")], by = "feature_id", all = TRUE
  )
  by_feature$method <- mode_label
  by_feature$drift_and_batch <- by_feature$selected_for_drift & by_feature$selected_for_batch
  summary <- data.frame(
    method = mode_label,
    profiles_tested_for_drift = nrow(drift_table),
    profiles_flagged_for_drift = sum(drift_table$selected_for_drift),
    proportion_profiles_flagged_for_drift = mean(drift_table$selected_for_drift),
    unique_features_detrended = length(drift_features),
    proportion_unique_features_detrended = length(drift_features) / nrow(x_mat),
    features_tested_for_batch = nrow(batch_table),
    features_selected_for_batch = length(batch_features),
    proportion_features_selected_for_batch = length(batch_features) / nrow(x_mat),
    drift_batch_overlap_features = length(intersect(drift_features, batch_features)),
    samples_modified_by_pqn = sum(pqn_changed),
    proportion_samples_modified_by_pqn = mean(pqn_changed),
    inferred_segments = length(unique(batch_vec)),
    inferred_changepoints = paste(which(diff(as.integer(factor(batch_vec, levels = unique(batch_vec)))) != 0), collapse = ";"),
    pelt_penalty = if (auto_batch) attr(batch_vec, "pelt_penalty") else NA_real_,
    test = params$test, acorr_fdr = params$acorr_fdr, anova_fdr = params$anova_fdr,
    normalization = params$normalization, spline_method = params$spline_method,
    scale_by_batch = params$scale_by_batch,
    smoother_fallback_count = NA_integer_,
    smoother_fallback_note = "not exposed by the installed WiNN API; stage output is validated directly",
    diagnostic_matrix_equal = NA,
    stringsAsFactors = FALSE
  )
  list(data = pmax(final, 0), drift = drift_table, batch = batch_table, by_feature = by_feature, summary = summary, detected_batch = batch_vec)
}

weighted_pc_r2_general <- function(mat, target, target_type = NULL, n_pcs = 5L) {
  if (!exists("weighted_pc_r2_explicit", mode = "function", inherits = TRUE)) {
    stop("Source scripts/weighted_pc_r2.R before using weighted_pc_r2_general().")
  }
  weighted_pc_r2_explicit(
    mat = mat,
    target = target,
    target_type = target_type,
    n_pcs = n_pcs
  )
}

feature_batch_associations <- function(mat, metadata, method_label) {
  ids <- metadata$sample_id[!metadata$is_qc & metadata$sample_id %in% colnames(mat)]
  md <- metadata[match(ids, metadata$sample_id), , drop = FALSE]
  batch <- factor(md$batch)
  z <- log1p(mat[, ids, drop = FALSE])
  rows <- lapply(seq_len(nrow(z)), function(i) {
    fit <- stats::aov(z[i, ] ~ batch)
    tab <- summary(fit)[[1]]
    ss <- tab$`Sum Sq`
    data.frame(method = method_label, feature_id = rownames(z)[i], f_statistic = tab$`F value`[1], p_value = tab$`Pr(>F)`[1], eta_squared = ss[1] / sum(ss), stringsAsFactors = FALSE)
  })
  out <- dplyr::bind_rows(rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out$is_batch_associated <- is.finite(out$p_adj) & out$p_adj < 0.05
  out
}

qc_metrics_by_feature <- function(mat, hidden_ids, method_label) {
  ids <- hidden_ids[hidden_ids %in% colnames(mat)]
  q <- mat[, ids, drop = FALSE]
  means <- rowMeans(q)
  sds <- apply(q, 1, stats::sd)
  data.frame(method = method_label, feature_id = rownames(q), n_hidden_qc = length(ids), mean = means, sd = sds, cv_percent = 100 * sds / abs(means), stringsAsFactors = FALSE)
}

qc_pairwise_metrics <- function(mat, hidden_ids, method_label) {
  ids <- hidden_ids[hidden_ids %in% colnames(mat)]
  out <- list()
  for (cor_method in c("pearson", "spearman")) {
    cm <- suppressWarnings(stats::cor(log1p(mat[, ids, drop = FALSE]), method = cor_method))
    pairs <- which(upper.tri(cm), arr.ind = TRUE)
    out[[cor_method]] <- data.frame(method = method_label, correlation_method = cor_method, sample_id_1 = colnames(cm)[pairs[, 1]], sample_id_2 = colnames(cm)[pairs[, 2]], correlation = cm[pairs], stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(out)
}

correction_magnitude_table <- function(mat, raw, method_label) {
  features <- intersect(rownames(raw), rownames(mat))
  samples <- intersect(colnames(raw), colnames(mat))
  delta <- abs(log1p(mat[features, samples, drop = FALSE]) - log1p(raw[features, samples, drop = FALSE]))
  data.frame(
    method = method_label,
    feature_id = features,
    median_absolute_log1p_change = apply(delta, 1, stats::median),
    p90_absolute_log1p_change = apply(delta, 1, stats::quantile, probs = 0.9, names = FALSE),
    proportion_entries_materially_changed = rowMeans(delta > 1e-6),
    stringsAsFactors = FALSE
  )
}

effect_concordance <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L) return(c(pearson = NA, spearman = NA, sign_concordance = NA, median_abs_difference = NA, top_quartile_sign_concordance = NA))
  top <- abs((a + b) / 2) >= stats::quantile(abs((a + b) / 2), 0.75, na.rm = TRUE)
  c(
    pearson = stats::cor(a[ok], b[ok], method = "pearson"),
    spearman = stats::cor(a[ok], b[ok], method = "spearman"),
    sign_concordance = mean(sign(a[ok]) == sign(b[ok])),
    median_abs_difference = stats::median(abs(a[ok] - b[ok])),
    top_quartile_sign_concordance = mean(sign(a[ok & top]) == sign(b[ok & top]))
  )
}
