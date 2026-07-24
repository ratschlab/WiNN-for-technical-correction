# Full evaluation engine for one canonical repeated-simulation unit.

canonical_write_csv_gz <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  utils::write.csv(value, connection, row.names = FALSE, quote = TRUE, na = "")
  invisible(path)
}

canonical_metric_row <- function(
  seed_id, method, metric, value, direction, denominator, units,
  status = "calculated", notes = ""
) {
  data.frame(
    seed_id = seed_id,
    method = method,
    metric = metric,
    value = as.numeric(value),
    metric_direction = direction,
    denominator = as.numeric(denominator),
    units = units,
    status = status,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

canonical_gate_stats <- function(predicted, actual) {
  predicted <- as.logical(predicted)
  actual <- as.logical(actual)
  tp <- sum(predicted & actual)
  tn <- sum(!predicted & !actual)
  fp <- sum(predicted & !actual)
  fn <- sum(!predicted & actual)
  data.frame(
    true_positive = tp,
    true_negative = tn,
    false_positive = fp,
    false_negative = fn,
    sensitivity = if (tp + fn > 0L) tp / (tp + fn) else NA_real_,
    specificity = if (tn + fp > 0L) tn / (tn + fp) else NA_real_,
    precision = if (tp + fp > 0L) tp / (tp + fp) else NA_real_,
    false_positive_rate = if (fp + tn > 0L) fp / (fp + tn) else NA_real_,
    stringsAsFactors = FALSE
  )
}

compute_canonical_seed_metrics <- function(
  seed_id,
  method_matrices,
  raw,
  truth,
  meta,
  hidden_ids,
  feature_truth,
  runtime,
  ablation_diagnostics,
  ablation_configurations,
  metric_dictionary,
  run_dir,
  log_line = function(...) invisible(NULL)
) {
  metric_dir <- file.path(run_dir, "metrics")
  dir.create(metric_dir, recursive = TRUE, showWarnings = FALSE)
  method_order <- names(method_matrices)
  study_ids <- meta$sample_id[!meta$is_qc]
  study_meta <- meta[match(study_ids, meta$sample_id), , drop = FALSE]
  truth_log <- log1p(truth[, study_ids, drop = FALSE])
  tolerance <- 1e-10

  truth_annotation <- feature_truth[match(rownames(raw), feature_truth$metabolite), , drop = FALSE]
  if (anyNA(truth_annotation$metabolite)) stop("Feature truth is misaligned.", call. = FALSE)
  drift_truth <- as.logical(truth_annotation$drift_effect_applied_any_plate)
  batch_truth <- as.logical(truth_annotation$batch_effect_applied)
  unaffected_truth <- !(drift_truth | batch_truth)
  affected_truth <- drift_truth | batch_truth

  metric_rows <- list()
  qc_rows <- list()
  profile_rows <- list()
  gam_rows <- list()
  ljung_rows <- list()
  feature_recovery_rows <- list()
  sample_recovery_rows <- list()
  gate_rows <- list()

  # Exact aliases reuse expensive evaluation calculations but retain separate
  # method rows.  Equality is revalidated by the seed runner before this point.
  alias_source <- c(
    C0_RAW = "Raw",
    C4_FULL_FIXED = "WINN default (no QC)",
    G_SS = "WINN default (no QC)"
  )
  calculation_cache <- list()

  for (method in method_order) {
    log_line("Computing full metrics for ", method, ".")
    matrix <- method_matrices[[method]]
    alias <- if (method %in% names(alias_source)) alias_source[[method]] else NA_character_
    use_alias <- !is.na(alias) && alias %in% names(calculation_cache) &&
      max(abs(matrix - method_matrices[[alias]]), na.rm = TRUE) <= 1e-8

    if (use_alias) {
      calculated <- calculation_cache[[alias]]
    } else {
      qc_cv <- qc_cv_ablation(matrix, hidden_ids)
      profile_correlations <- hidden_profile_correlations(matrix, hidden_ids)

      gam <- compute_metabolite_segment_gam(
        matrix[, study_ids, drop = FALSE], study_meta, method,
        sample_id_col = "sample_id", order_col = "run_order",
        batch_col = "batch", transform_fun = log1p,
        min_obs = 6L, k_max = 6L
      )
      ljung <- compute_metabolite_segment_ljung_box(
        matrix[, study_ids, drop = FALSE], study_meta, method,
        sample_id_col = "sample_id", order_col = "run_order",
        batch_col = "batch", transform_fun = log1p,
        min_obs = 8L, max_lag = 10L
      ) |>
        dplyr::group_by(method, batch, segment_id) |>
        dplyr::mutate(
          p_adj = stats::p.adjust(p_value, method = "BH"),
          is_autocorrelated = is.finite(p_adj) & p_adj < 0.05
        ) |>
        dplyr::ungroup()

      candidate_log <- log1p(matrix[, study_ids, drop = FALSE])
      if (any(!is.finite(candidate_log))) stop(method, " has nonfinite log1p values.", call. = FALSE)
      feature_icc <- vapply(seq_len(nrow(candidate_log)), function(index) {
        icc_a1_ablation(cbind(truth_log[index, ], candidate_log[index, ]))
      }, numeric(1))
      feature_cor <- vapply(seq_len(nrow(candidate_log)), function(index) {
        safe_cor_ablation(truth_log[index, ], candidate_log[index, ], method = "pearson")
      }, numeric(1))
      sample_cor <- vapply(seq_len(ncol(candidate_log)), function(index) {
        safe_cor_ablation(truth_log[, index], candidate_log[, index], method = "pearson")
      }, numeric(1))
      batch_r2 <- weighted_pc_r2_explicit(
        matrix[, study_ids, drop = FALSE], study_meta$batch,
        target_type = "categorical"
      )
      calculated <- list(
        qc_cv = qc_cv,
        profile_correlations = profile_correlations,
        gam = gam,
        ljung = ljung,
        feature_icc = feature_icc,
        feature_cor = feature_cor,
        sample_cor = sample_cor,
        rmse = sqrt(mean((candidate_log - truth_log)^2)),
        batch_r2 = batch_r2
      )
    }
    calculation_cache[[method]] <- calculated

    qc_table <- data.frame(
      seed_id = seed_id, method = method,
      feature_id = names(calculated$qc_cv),
      cv_percent = as.numeric(calculated$qc_cv),
      stringsAsFactors = FALSE
    )
    qc_rows[[method]] <- qc_table
    profile <- calculated$profile_correlations
    profile$seed_id <- seed_id
    profile$method <- method
    profile_rows[[method]] <- profile[, c("seed_id", "method", "correlation_method", "value")]
    gam <- calculated$gam
    gam$method <- method
    gam$seed_id <- seed_id
    gam_rows[[method]] <- gam[, c("seed_id", names(gam)[names(gam) != "seed_id"])]
    ljung <- calculated$ljung
    ljung$method <- method
    ljung$seed_id <- seed_id
    ljung_rows[[method]] <- ljung[, c("seed_id", names(ljung)[names(ljung) != "seed_id"])]
    feature_recovery_rows[[method]] <- data.frame(
      seed_id = seed_id, method = method, feature_id = rownames(raw),
      icc_a1 = calculated$feature_icc,
      pearson = calculated$feature_cor,
      stringsAsFactors = FALSE
    )
    sample_recovery_rows[[method]] <- data.frame(
      seed_id = seed_id, method = method, sample_id = study_ids,
      pearson = calculated$sample_cor,
      stringsAsFactors = FALSE
    )

    n_ljung <- sum(is.finite(ljung$p_value))
    metric_rows[[paste0(method, "::universal")]] <- dplyr::bind_rows(
      canonical_metric_row(seed_id, method, "heldout_qc_cv_mean", mean(calculated$qc_cv, na.rm = TRUE), "lower", sum(is.finite(calculated$qc_cv)), "percent"),
      canonical_metric_row(seed_id, method, "heldout_qc_cv_median", stats::median(calculated$qc_cv, na.rm = TRUE), "lower", sum(is.finite(calculated$qc_cv)), "percent"),
      canonical_metric_row(seed_id, method, "heldout_qc_cv_sd", stats::sd(calculated$qc_cv, na.rm = TRUE), "lower", sum(is.finite(calculated$qc_cv)), "percent"),
      canonical_metric_row(seed_id, method, "hidden_qc_profile_pearson_median", stats::median(profile$value[profile$correlation_method == "pearson"], na.rm = TRUE), "higher", sum(profile$correlation_method == "pearson"), "correlation"),
      canonical_metric_row(seed_id, method, "hidden_qc_profile_spearman_median", stats::median(profile$value[profile$correlation_method == "spearman"], na.rm = TRUE), "higher", sum(profile$correlation_method == "spearman"), "correlation"),
      canonical_metric_row(seed_id, method, "residual_gam_deviance_mean", mean(gam$explained, na.rm = TRUE), "lower", sum(is.finite(gam$explained)), "proportion"),
      canonical_metric_row(seed_id, method, "residual_gam_deviance_median", stats::median(gam$explained, na.rm = TRUE), "lower", sum(is.finite(gam$explained)), "proportion"),
      canonical_metric_row(seed_id, method, "residual_ljung_box_significant", sum(ljung$is_autocorrelated, na.rm = TRUE), "lower", n_ljung, "profiles"),
      canonical_metric_row(seed_id, method, "residual_ljung_box_proportion", if (n_ljung) mean(ljung$is_autocorrelated[is.finite(ljung$p_value)]) else NA_real_, "lower", n_ljung, "proportion", notes = "BH adjusted within supplied batch segment."),
      canonical_metric_row(seed_id, method, "batch_weighted_pc_r2", calculated$batch_r2, "lower", length(study_ids), "weighted R-squared"),
      canonical_metric_row(seed_id, method, "truth_sample_profile_pearson_mean", mean(calculated$sample_cor, na.rm = TRUE), "higher", sum(is.finite(calculated$sample_cor)), "correlation"),
      canonical_metric_row(seed_id, method, "truth_metabolite_profile_icc_mean", mean(calculated$feature_icc, na.rm = TRUE), "higher", sum(is.finite(calculated$feature_icc)), "ICC(A,1)"),
      canonical_metric_row(seed_id, method, "truth_log1p_rmse", calculated$rmse, "lower", length(truth_log), "log1p intensity"),
      canonical_metric_row(seed_id, method, "features_retained", nrow(matrix), "higher", nrow(raw), "features"),
      canonical_metric_row(seed_id, method, "samples_retained", ncol(matrix), "higher", ncol(raw), "samples"),
      canonical_metric_row(seed_id, method, "feature_coverage", nrow(matrix) / nrow(raw), "higher", nrow(raw), "proportion"),
      canonical_metric_row(seed_id, method, "sample_coverage", ncol(matrix) / ncol(raw), "higher", ncol(raw), "proportion"),
      canonical_metric_row(seed_id, method, "runtime_sec", runtime$runtime_sec[match(method, runtime$method)], "lower", 1, "seconds")
    )

    diagnostic_method <- if (method == "WINN default (no QC)") "C4_FULL_FIXED" else method
    diagnostics <- ablation_diagnostics[[diagnostic_method]]
    configuration <- ablation_configurations[[diagnostic_method]]
    drift_supported <- is.list(diagnostics) && is.data.frame(diagnostics$drift) &&
      nrow(diagnostics$drift) > 0L && !is.null(configuration$drift_gate) &&
      configuration$drift_gate != "none"
    batch_supported <- is.list(diagnostics) && is.data.frame(diagnostics$batch) &&
      nrow(diagnostics$batch) > 0L && !is.null(configuration$batch_gate) &&
      configuration$batch_gate != "none"

    drift_stats <- canonical_gate_stats(rep(FALSE, nrow(raw)), drift_truth)
    batch_stats <- canonical_gate_stats(rep(FALSE, nrow(raw)), batch_truth)
    if (drift_supported) {
      drift_predicted <- rownames(raw) %in% unique(
        diagnostics$drift$feature_id[diagnostics$drift$actual_correction]
      )
      drift_stats <- canonical_gate_stats(drift_predicted, drift_truth)
    }
    if (batch_supported) {
      batch_predicted <- rownames(raw) %in%
        diagnostics$batch$feature_id[diagnostics$batch$actual_correction]
      batch_stats <- canonical_gate_stats(batch_predicted, batch_truth)
    }
    changed_features <- rowSums(abs(log1p(matrix) - log1p(raw)) > tolerance) > 0L
    gate_rows[[method]] <- data.frame(
      seed_id = seed_id, method = method,
      drift_supported = drift_supported, batch_supported = batch_supported,
      drift_stats, batch_stats,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    # Disambiguate the two sets of confusion-count column names.
    names(gate_rows[[method]])[5:12] <- paste0("drift_", names(gate_rows[[method]])[5:12])
    names(gate_rows[[method]])[13:20] <- paste0("batch_", names(gate_rows[[method]])[13:20])

    unsupported <- function(supported) if (supported) "calculated_known_feature_truth" else "not_supported_no_gate_diagnostics"
    gate_metric_rows <- dplyr::bind_rows(
      canonical_metric_row(seed_id, method, "drift_gate_sensitivity", if (drift_supported) drift_stats$sensitivity else NA_real_, "higher", sum(drift_truth), "proportion", unsupported(drift_supported)),
      canonical_metric_row(seed_id, method, "drift_gate_specificity", if (drift_supported) drift_stats$specificity else NA_real_, "higher", sum(!drift_truth), "proportion", unsupported(drift_supported)),
      canonical_metric_row(seed_id, method, "drift_gate_precision", if (drift_supported) drift_stats$precision else NA_real_, "higher", if (drift_supported) drift_stats$true_positive + drift_stats$false_positive else NA_real_, "proportion", unsupported(drift_supported)),
      canonical_metric_row(seed_id, method, "drift_gate_false_positive_rate", if (drift_supported) drift_stats$false_positive_rate else NA_real_, "lower", sum(!drift_truth), "proportion", unsupported(drift_supported)),
      canonical_metric_row(seed_id, method, "batch_gate_sensitivity", if (batch_supported) batch_stats$sensitivity else NA_real_, "higher", sum(batch_truth), "proportion", unsupported(batch_supported)),
      canonical_metric_row(seed_id, method, "batch_gate_specificity", if (batch_supported) batch_stats$specificity else NA_real_, "higher", sum(!batch_truth), "proportion", unsupported(batch_supported)),
      canonical_metric_row(seed_id, method, "batch_gate_precision", if (batch_supported) batch_stats$precision else NA_real_, "higher", if (batch_supported) batch_stats$true_positive + batch_stats$false_positive else NA_real_, "proportion", unsupported(batch_supported)),
      canonical_metric_row(seed_id, method, "batch_gate_false_positive_rate", if (batch_supported) batch_stats$false_positive_rate else NA_real_, "lower", sum(!batch_truth), "proportion", unsupported(batch_supported)),
      canonical_metric_row(seed_id, method, "unaffected_features_modified", mean(changed_features[unaffected_truth]), "lower", sum(unaffected_truth), "proportion", notes = "Final-matrix intervention; dilution truth is not separately annotated."),
      canonical_metric_row(seed_id, method, "affected_features_left_unchanged", mean(!changed_features[affected_truth]), "lower", sum(affected_truth), "proportion")
    )
    metric_rows[[paste0(method, "::gate")]] <- gate_metric_rows
  }

  metrics <- dplyr::bind_rows(metric_rows)
  metrics$method <- factor(metrics$method, levels = method_order)
  metrics <- metrics[order(metrics$method, match(metrics$metric, metric_dictionary$metric)), , drop = FALSE]
  metrics$method <- as.character(metrics$method)
  if (nrow(metrics) != length(method_order) * nrow(metric_dictionary)) {
    stop("Method metric table is incomplete.", call. = FALSE)
  }
  if (!identical(sort(unique(metrics$metric)), sort(metric_dictionary$metric))) {
    stop("Method metric names differ from the frozen dictionary.", call. = FALSE)
  }

  utils::write.csv(metrics, file.path(run_dir, "method_metrics.csv"), row.names = FALSE, quote = TRUE, na = "")
  canonical_write_csv_gz(dplyr::bind_rows(qc_rows), file.path(metric_dir, "heldout_qc_cv_by_feature.csv.gz"))
  canonical_write_csv_gz(dplyr::bind_rows(profile_rows), file.path(metric_dir, "hidden_profile_correlations.csv.gz"))
  canonical_write_csv_gz(dplyr::bind_rows(gam_rows), file.path(metric_dir, "run_order_gam_by_feature_batch.csv.gz"))
  canonical_write_csv_gz(dplyr::bind_rows(ljung_rows), file.path(metric_dir, "ljung_box_by_feature_batch.csv.gz"))
  canonical_write_csv_gz(dplyr::bind_rows(feature_recovery_rows), file.path(metric_dir, "truth_recovery_by_feature.csv.gz"))
  canonical_write_csv_gz(dplyr::bind_rows(sample_recovery_rows), file.path(metric_dir, "truth_recovery_by_sample.csv.gz"))
  gate_table <- dplyr::bind_rows(gate_rows)
  utils::write.csv(gate_table, file.path(metric_dir, "gate_truth_confusion.csv"), row.names = FALSE, quote = TRUE, na = "")

  list(
    method_metrics = metrics,
    gate_truth = gate_table,
    detailed_row_counts = data.frame(
      table = c(
        "heldout_qc_cv_by_feature", "hidden_profile_correlations",
        "run_order_gam_by_feature_batch", "ljung_box_by_feature_batch",
        "truth_recovery_by_feature", "truth_recovery_by_sample"
      ),
      rows = c(
        nrow(dplyr::bind_rows(qc_rows)), nrow(dplyr::bind_rows(profile_rows)),
        nrow(dplyr::bind_rows(gam_rows)), nrow(dplyr::bind_rows(ljung_rows)),
        nrow(dplyr::bind_rows(feature_recovery_rows)),
        nrow(dplyr::bind_rows(sample_recovery_rows))
      ),
      stringsAsFactors = FALSE
    )
  )
}
