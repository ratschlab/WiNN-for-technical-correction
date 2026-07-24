# Metric helpers for the partial-confounding pilot.  These functions are
# evaluation-only and must never be called from a correction-method closure.

.partial_metric_safe_cor <- function(x, y, method = "pearson") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || stats::sd(x[keep]) <= 0 || stats::sd(y[keep]) <= 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[keep], y[keep], method = method))
}

#' Estimate the phenotype coefficient for every feature
#'
#' The fitted model is prespecified as phenotype + plate + scaled within-plate
#' order.  Rank deficiency returns an explicit non-estimable result rather than
#' silently dropping or forcing a coefficient.
fit_partial_confounding_effects <- function(
  mat,
  allocation,
  feature_truth,
  alpha = 0.05
) {
  required_allocation <- c(
    "sample_id", "plate", "sample_type", "phenotype_numeric",
    "within_plate_position_scaled"
  )
  missing_allocation <- setdiff(required_allocation, names(allocation))
  if (length(missing_allocation)) {
    stop(
      "Allocation lacks effect-model field(s): ",
      paste(missing_allocation, collapse = ", "), call. = FALSE
    )
  }
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("Effect input must be a numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(mat)) || any(mat <= -1)) {
    stop("Effect input must be finite and strictly greater than -1.",
         call. = FALSE)
  }
  study <- allocation[allocation$sample_type == "sample", , drop = FALSE]
  study <- study[match(colnames(mat)[colnames(mat) %in% study$sample_id],
                       study$sample_id), , drop = FALSE]
  if (anyNA(study$sample_id) || nrow(study) != 450L) {
    stop("Expected all 450 study samples in effect-model order.", call. = FALSE)
  }
  study_ids <- study$sample_id
  phenotype <- as.numeric(study$phenotype_numeric)
  plate <- factor(study$plate, levels = unique(allocation$plate))
  within_plate_position_scaled <- as.numeric(
    study$within_plate_position_scaled
  )
  design <- stats::model.matrix(
    ~ phenotype + plate + within_plate_position_scaled
  )
  design_rank <- qr(design)$rank
  design_columns <- ncol(design)
  singular_values <- svd(design, nu = 0L, nv = 0L)$d
  condition_number <- if (min(singular_values) <= 0) {
    Inf
  } else {
    max(singular_values) / min(singular_values)
  }
  phenotype_column <- match("phenotype", colnames(design))
  estimable <- design_rank == design_columns && !is.na(phenotype_column)

  base <- data.frame(
    feature_id = rownames(mat),
    injected_effect_log = feature_truth$phenotype_effect_log[
      match(rownames(mat), feature_truth$metabolite)
    ],
    responsive = feature_truth$responsive[
      match(rownames(mat), feature_truth$metabolite)
    ],
    recovered_effect_log1p = NA_real_,
    standard_error = NA_real_,
    t_statistic = NA_real_,
    p_value = NA_real_,
    p_adjust_bh = NA_real_,
    significant_bh_0p05 = FALSE,
    estimable = estimable,
    stringsAsFactors = FALSE
  )
  if (anyNA(base$injected_effect_log) || anyNA(base$responsive)) {
    stop("Feature truth does not align to the effect matrix.", call. = FALSE)
  }

  if (estimable) {
    response <- t(log1p(mat[, study_ids, drop = FALSE]))
    xtx_inverse <- tryCatch(
      solve(crossprod(design)),
      error = function(e) NULL
    )
    if (!is.null(xtx_inverse)) {
      coefficients <- xtx_inverse %*% crossprod(design, response)
      residuals <- response - design %*% coefficients
      residual_df <- nrow(design) - design_rank
      residual_variance <- colSums(residuals^2) / residual_df
      effect <- as.numeric(coefficients[phenotype_column, ])
      standard_error <- sqrt(
        pmax(0, residual_variance * xtx_inverse[phenotype_column,
                                               phenotype_column])
      )
      statistic <- effect / standard_error
      p_value <- 2 * stats::pt(abs(statistic), df = residual_df,
                               lower.tail = FALSE)
      p_adjust <- stats::p.adjust(p_value, method = "BH")
      base$recovered_effect_log1p <- effect
      base$standard_error <- standard_error
      base$t_statistic <- statistic
      base$p_value <- p_value
      base$p_adjust_bh <- p_adjust
      base$significant_bh_0p05 <- is.finite(p_adjust) & p_adjust < alpha
    } else {
      base$estimable <- FALSE
      estimable <- FALSE
    }
  }

  design_diagnostics <- data.frame(
    formula = "~ phenotype + plate + within_plate_position_scaled",
    n_study = nrow(design),
    n_columns = design_columns,
    rank = design_rank,
    full_rank = design_rank == design_columns,
    phenotype_column_present = !is.na(phenotype_column),
    estimable = estimable,
    condition_number = condition_number,
    minimum_singular_value = min(singular_values),
    stringsAsFactors = FALSE
  )
  list(estimates = base, design = design_diagnostics)
}

#' Summarize biological-effect recovery relative to injected and clean effects
summarise_partial_effect_recovery <- function(recovered, clean_reference) {
  if (!identical(recovered$feature_id, clean_reference$feature_id)) {
    stop("Recovered and clean-reference effect tables are misaligned.",
         call. = FALSE)
  }
  true <- recovered$injected_effect_log
  clean <- clean_reference$recovered_effect_log1p
  estimate <- recovered$recovered_effect_log1p
  responsive <- recovered$responsive
  null <- !responsive
  discoveries <- recovered$significant_bh_0p05
  estimable <- all(recovered$estimable) && all(clean_reference$estimable)

  slope <- intercept <- NA_real_
  if (estimable && stats::sd(true) > 0) {
    slope_fit <- stats::lm(estimate ~ true)
    intercept <- unname(stats::coef(slope_fit)[[1L]])
    slope <- unname(stats::coef(slope_fit)[[2L]])
  }
  n_discoveries <- sum(discoveries, na.rm = TRUE)
  true_positives <- sum(discoveries & responsive, na.rm = TRUE)
  false_positives <- sum(discoveries & null, na.rm = TRUE)
  ratio_injected <- estimate[responsive] / true[responsive]
  ratio_clean <- estimate[responsive] / clean[responsive]

  metrics <- c(
    coefficient_pearson_vs_injected_all = .partial_metric_safe_cor(
      estimate, true, "pearson"
    ),
    coefficient_spearman_vs_injected_all = .partial_metric_safe_cor(
      estimate, true, "spearman"
    ),
    coefficient_pearson_vs_injected_responsive = .partial_metric_safe_cor(
      estimate[responsive], true[responsive], "pearson"
    ),
    coefficient_pearson_vs_clean_all = .partial_metric_safe_cor(
      estimate, clean, "pearson"
    ),
    effect_rmse_vs_injected_all = sqrt(mean((estimate - true)^2,
                                             na.rm = TRUE)),
    effect_rmse_vs_injected_responsive = sqrt(mean(
      (estimate[responsive] - true[responsive])^2, na.rm = TRUE
    )),
    effect_rmse_vs_clean_all = sqrt(mean((estimate - clean)^2,
                                          na.rm = TRUE)),
    effect_rmse_vs_clean_responsive = sqrt(mean(
      (estimate[responsive] - clean[responsive])^2, na.rm = TRUE
    )),
    recovered_vs_injected_slope = slope,
    recovered_vs_injected_intercept = intercept,
    median_attenuation_ratio_injected_responsive = stats::median(
      ratio_injected, na.rm = TRUE
    ),
    median_attenuation_ratio_clean_responsive = stats::median(
      ratio_clean, na.rm = TRUE
    ),
    sign_concordance_responsive = mean(
      sign(estimate[responsive]) == sign(true[responsive]), na.rm = TRUE
    ),
    n_discoveries_bh_0p05 = n_discoveries,
    true_positive_rate = true_positives / sum(responsive),
    empirical_fdr = if (n_discoveries) {
      false_positives / n_discoveries
    } else {
      0
    },
    null_feature_false_positive_rate = false_positives / sum(null)
  )
  if (!estimable) metrics[] <- NA_real_
  data.frame(
    metric = names(metrics),
    value = unname(metrics),
    status = if (estimable) "calculated" else "not_estimable_rank_deficient",
    stringsAsFactors = FALSE
  )
}

#' Calculate matrix recovery against the clean log1p truth
partial_clean_matrix_recovery <- function(candidate, clean, study_ids) {
  candidate_log <- log1p(candidate[, study_ids, drop = FALSE])
  clean_log <- log1p(clean[, study_ids, drop = FALSE])
  if (any(!is.finite(candidate_log)) || any(!is.finite(clean_log))) {
    stop("Nonfinite log1p values in clean-matrix recovery.", call. = FALSE)
  }
  sample_correlations <- vapply(seq_len(ncol(candidate_log)), function(index) {
    .partial_metric_safe_cor(
      candidate_log[, index], clean_log[, index], "pearson"
    )
  }, numeric(1))
  feature_correlations <- vapply(seq_len(nrow(candidate_log)), function(index) {
    .partial_metric_safe_cor(
      candidate_log[index, ], clean_log[index, ], "pearson"
    )
  }, numeric(1))
  feature_icc <- vapply(seq_len(nrow(candidate_log)), function(index) {
    icc_a1_ablation(cbind(
      clean_log[index, ], candidate_log[index, ]
    ))
  }, numeric(1))
  feature_rmse <- sqrt(rowMeans((candidate_log - clean_log)^2))
  sample_rmse <- sqrt(colMeans((candidate_log - clean_log)^2))
  overall_rmse <- sqrt(mean((candidate_log - clean_log)^2))
  list(
    summary = data.frame(
      metric = c(
        "clean_log1p_rmse_study",
        "clean_sample_profile_pearson_mean",
        "clean_sample_profile_pearson_median",
        "clean_feature_profile_pearson_mean",
        "clean_feature_profile_pearson_median",
        "clean_feature_icc_a1_mean",
        "clean_feature_icc_a1_median"
      ),
      value = c(
        overall_rmse,
        mean(sample_correlations, na.rm = TRUE),
        stats::median(sample_correlations, na.rm = TRUE),
        mean(feature_correlations, na.rm = TRUE),
        stats::median(feature_correlations, na.rm = TRUE),
        mean(feature_icc, na.rm = TRUE),
        stats::median(feature_icc, na.rm = TRUE)
      ),
      status = "calculated",
      stringsAsFactors = FALSE
    ),
    by_feature = data.frame(
      feature_id = rownames(candidate_log),
      pearson = feature_correlations,
      icc_a1 = feature_icc,
      rmse_log1p = feature_rmse,
      stringsAsFactors = FALSE
    ),
    by_sample = data.frame(
      sample_id = colnames(candidate_log),
      pearson = sample_correlations,
      rmse_log1p = sample_rmse,
      stringsAsFactors = FALSE
    )
  )
}

#' Summarize intervention rates exposed by winn::winn_ablation diagnostics
summarise_partial_winn_interventions <- function(diagnostics) {
  drift <- diagnostics$drift
  batch <- diagnostics$batch
  pqn <- diagnostics$pqn
  outlier <- diagnostics$outlier
  data.frame(
    outlier_entries_changed = if (is.list(outlier) &&
                                  !is.null(outlier$entries_changed)) {
      outlier$entries_changed
    } else {
      NA_real_
    },
    outlier_proportion_entries_changed = if (is.list(outlier) &&
                                             !is.null(outlier$proportion_entries_changed)) {
      outlier$proportion_entries_changed
    } else {
      NA_real_
    },
    drift_profiles_tested = if (is.data.frame(drift)) nrow(drift) else 0L,
    drift_profiles_selected = if (is.data.frame(drift) &&
                                  "selective_significance" %in% names(drift)) {
      sum(drift$selective_significance, na.rm = TRUE)
    } else {
      0L
    },
    drift_profiles_corrected = if (is.data.frame(drift) &&
                                   "actual_correction" %in% names(drift)) {
      sum(drift$actual_correction, na.rm = TRUE)
    } else {
      0L
    },
    drift_proportion_corrected = if (is.data.frame(drift) && nrow(drift) &&
                                     "actual_correction" %in% names(drift)) {
      mean(drift$actual_correction, na.rm = TRUE)
    } else {
      NA_real_
    },
    drift_unique_features_corrected = if (is.data.frame(drift) &&
                                         "actual_correction" %in% names(drift)) {
      length(unique(drift$feature_id[drift$actual_correction]))
    } else {
      0L
    },
    batch_features_tested = if (is.data.frame(batch)) nrow(batch) else 0L,
    batch_features_selected = if (is.data.frame(batch) &&
                                  "selective_significance" %in% names(batch)) {
      sum(batch$selective_significance, na.rm = TRUE)
    } else {
      0L
    },
    batch_features_corrected = if (is.data.frame(batch) &&
                                   "actual_correction" %in% names(batch)) {
      sum(batch$actual_correction, na.rm = TRUE)
    } else {
      0L
    },
    batch_proportion_corrected = if (is.data.frame(batch) && nrow(batch) &&
                                     "actual_correction" %in% names(batch)) {
      mean(batch$actual_correction, na.rm = TRUE)
    } else {
      NA_real_
    },
    pqn_samples_tested = if (is.data.frame(pqn)) nrow(pqn) else 0L,
    pqn_samples_altered = if (is.data.frame(pqn) && "altered" %in% names(pqn)) {
      sum(pqn$altered, na.rm = TRUE)
    } else {
      0L
    },
    stringsAsFactors = FALSE
  )
}
