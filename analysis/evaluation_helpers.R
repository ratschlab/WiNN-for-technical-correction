source(
  file.path(source_root, "scripts", "run_order_drift_helpers.R"),
  local = FALSE
)

evaluation_metric_row <- function(metric, value, direction, denominator, units = "", notes = "") {
  data.frame(
    metric = metric, value = as.numeric(value), direction = direction,
    denominator = as.character(denominator), units = units, notes = notes,
    stringsAsFactors = FALSE
  )
}

evaluation_study_panel <- function(matrix, prepared) {
  meta <- prepared$meta
  study_flag <- if ("is_study" %in% names(meta)) {
    as.logical(meta$is_study)
  } else if ("role" %in% names(meta)) {
    meta$role == "clinical"
  } else {
    !meta$sample_id %in% c(prepared$training_ids, prepared$hidden_ids, prepared$external_ids)
  }
  ids <- as.character(meta$sample_id[study_flag])
  ids <- ids[ids %in% colnames(matrix)]
  list(
    ids = ids,
    meta = meta[match(ids, meta$sample_id), , drop = FALSE]
  )
}

evaluation_reference_ids <- function(prepared) {
  if (length(prepared$external_ids)) as.character(prepared$external_ids) else as.character(prepared$hidden_ids)
}

evaluation_feature_repeatability <- function(matrix, meta, group_field) {
  panel <- evaluation_study_panel(matrix, list(
    meta = meta, training_ids = character(), hidden_ids = character(), external_ids = character()
  ))
  z <- log1p(matrix[, panel$ids, drop = FALSE])
  groups <- split(seq_along(panel$ids), panel$meta[[group_field]])
  groups <- groups[lengths(groups) >= 2L]
  if (!length(groups)) return(stats::setNames(rep(NA_real_, nrow(z)), rownames(z)))
  vapply(seq_len(nrow(z)), function(feature_index) {
    values <- z[feature_index, ]
    means <- vapply(groups, function(index) mean(values[index]), numeric(1))
    between <- stats::var(means)
    within_df <- sum(lengths(groups) - 1L)
    within <- if (within_df > 0L) {
      sum(vapply(groups, function(index) stats::var(values[index]) * (length(index) - 1L), numeric(1))) / within_df
    } else NA_real_
    if (!is.finite(between) || !is.finite(within) || between + within <= 0) NA_real_ else between / (between + within)
  }, numeric(1))
}

evaluation_group_replicates <- function(matrix, meta, group_field) {
  study_flag <- if ("is_study" %in% names(meta)) as.logical(meta$is_study) else meta$role == "clinical"
  ids <- as.character(meta$sample_id[study_flag])
  md <- meta[match(ids, meta$sample_id), , drop = FALSE]
  z <- log1p(matrix[, ids, drop = FALSE])
  groups <- split(seq_along(ids), md[[group_field]])
  groups <- groups[lengths(groups) >= 2L]
  if (!length(groups)) return(data.frame(group = character(), pearson = numeric(), icc_a1 = numeric()))
  do.call(rbind, lapply(names(groups), function(group_id) {
    index <- groups[[group_id]]
    correlations <- suppressWarnings(stats::cor(z[, index, drop = FALSE], use = "pairwise.complete.obs"))
    data.frame(
      group = group_id,
      pearson = stats::median(correlations[lower.tri(correlations)], na.rm = TRUE),
      icc_a1 = icc_a1_ablation(z[, index, drop = FALSE]),
      stringsAsFactors = FALSE
    )
  }))
}

evaluation_winn_counts <- function(details, n_features) {
  drift <- details$stage_decisions$drift
  batch <- details$stage_decisions$batch
  corrected <- function(value) {
    if (is.data.frame(value) && "actual_correction" %in% names(value)) {
      sum(as.logical(value$actual_correction), na.rm = TRUE)
    } else NA_real_
  }
  unique_drift <- if (is.data.frame(drift) && all(c("feature_id", "actual_correction") %in% names(drift))) {
    length(unique(drift$feature_id[as.logical(drift$actual_correction)]))
  } else NA_real_
  list(
    drift_profiles_corrected = corrected(drift),
    drift_features_corrected = unique_drift,
    batch_features_corrected = corrected(batch),
    drift_feature_proportion = unique_drift / n_features
  )
}

evaluate_release_matrix <- function(dataset_key, method_id, method_label, matrix, prepared, runtime_sec = NA_real_, details = list()) {
  matrix <- release_validate_output(matrix, prepared, method_id)
  study <- evaluation_study_panel(matrix, prepared)
  reference_ids <- evaluation_reference_ids(prepared)
  qc <- qc_cv_ablation(matrix, reference_ids)

  gam <- compute_metabolite_segment_gam(
    matrix[, study$ids, drop = FALSE], study$meta, method_label,
    sample_id_col = "sample_id", order_col = "run_order", batch_col = "batch",
    transform_fun = log1p, min_obs = 6L, k_max = 6L
  )
  lb <- compute_metabolite_segment_ljung_box(
    matrix[, study$ids, drop = FALSE], study$meta, method_label,
    sample_id_col = "sample_id", order_col = "run_order", batch_col = "batch",
    transform_fun = log1p, min_obs = 8L, max_lag = 10L
  )
  if (nrow(lb)) {
    grouping <- interaction(lb$method, lb$batch, lb$segment_id, drop = TRUE)
    lb$p_adj <- ave(lb$p_value, grouping, FUN = function(value) stats::p.adjust(value, method = "BH"))
    lb$is_autocorrelated <- is.finite(lb$p_adj) & lb$p_adj < 0.05
  }

  batch_r2 <- weighted_pc_r2_explicit(
    matrix[, study$ids, drop = FALSE], factor(study$meta$batch),
    target_type = "categorical"
  )
  counts <- evaluation_winn_counts(details, nrow(matrix))
  metrics <- rbind(
    evaluation_metric_row("heldout_qc_cv_mean", mean(qc, na.rm = TRUE), "lower", sum(is.finite(qc)), "percent"),
    evaluation_metric_row("heldout_qc_cv_sd", stats::sd(qc, na.rm = TRUE), "lower", sum(is.finite(qc)), "percent"),
    evaluation_metric_row("heldout_qc_cv_median", stats::median(qc, na.rm = TRUE), "lower", sum(is.finite(qc)), "percent"),
    evaluation_metric_row("residual_run_order_gam_mean_deviance", mean(gam$explained, na.rm = TRUE), "lower", sum(is.finite(gam$explained)), "proportion"),
    evaluation_metric_row("residual_run_order_gam_sd_deviance", stats::sd(gam$explained, na.rm = TRUE), "lower", sum(is.finite(gam$explained)), "proportion"),
    evaluation_metric_row("residual_ljung_box_proportion_significant", mean(lb$is_autocorrelated, na.rm = TRUE), "lower", sum(is.finite(lb$p_adj)), "proportion"),
    evaluation_metric_row("batch_weighted_pc_r2_categorical", batch_r2, "lower", length(study$ids), "weighted R-squared", "Batch is explicitly categorical."),
    evaluation_metric_row("retained_features", nrow(matrix), "higher", nrow(prepared$x), "features"),
    evaluation_metric_row("retained_samples", ncol(matrix), "higher", ncol(prepared$x), "samples"),
    evaluation_metric_row("runtime_sec", runtime_sec, "lower", 1, "seconds"),
    evaluation_metric_row("drift_profiles_corrected", counts$drift_profiles_corrected, "descriptive", nrow(matrix) * length(unique(prepared$meta$batch)), "profiles"),
    evaluation_metric_row("drift_features_corrected", counts$drift_features_corrected, "descriptive", nrow(matrix), "features"),
    evaluation_metric_row("batch_features_corrected", counts$batch_features_corrected, "descriptive", nrow(matrix), "features")
  )

  associations <- data.frame()
  replicates <- data.frame()
  if (identical(dataset_key, "simulation")) {
    truth <- log1p(prepared$truth[, study$ids, drop = FALSE])
    candidate <- log1p(matrix[, study$ids, drop = FALSE])
    feature_icc <- vapply(seq_len(nrow(candidate)), function(index) {
      icc_a1_ablation(cbind(truth[index, ], candidate[index, ]))
    }, numeric(1))
    sample_cor <- vapply(seq_len(ncol(candidate)), function(index) {
      safe_cor_ablation(truth[, index], candidate[, index])
    }, numeric(1))
    metrics <- rbind(metrics,
      evaluation_metric_row("truth_metabolite_profile_icc_mean", mean(feature_icc, na.rm = TRUE), "higher", sum(is.finite(feature_icc)), "ICC(A,1)"),
      evaluation_metric_row("truth_sample_profile_pearson_mean", mean(sample_cor, na.rm = TRUE), "higher", sum(is.finite(sample_cor)), "correlation"),
      evaluation_metric_row("truth_log1p_rmse", sqrt(mean((candidate - truth)^2)), "lower", length(candidate), "log1p intensity")
    )
  } else if (identical(dataset_key, "mtbls79")) {
    replicates <- mtbls79_replicate_metrics(matrix, prepared$meta)
    feature_icc <- mtbls79_feature_cross_batch_icc(matrix, prepared$meta)
    associations <- mtbls79_associations(matrix, prepared$meta, method_id)
    metrics <- rbind(metrics,
      evaluation_metric_row("biological_weighted_pc_r2", weighted_pc_r2_explicit(matrix[, study$ids], factor(study$meta$class), "categorical"), "higher", length(study$ids), "weighted R-squared"),
      evaluation_metric_row("metabolite_icc_a1_median", stats::median(feature_icc, na.rm = TRUE), "higher", sum(is.finite(feature_icc)), "ICC(A,1)"),
      evaluation_metric_row("sample_replicate_pearson_median", stats::median(replicates$pearson, na.rm = TRUE), "higher", nrow(replicates), "correlation"),
      evaluation_metric_row("biological_associated_features", sum(associations$significant), "context", nrow(associations), "features")
    )
  } else if (identical(dataset_key, "clinical_fiams")) {
    replicates <- evaluation_group_replicates(matrix, prepared$meta, "sample")
    repeatability <- evaluation_feature_repeatability(matrix, prepared$meta, "sample")
    metrics <- rbind(metrics,
      evaluation_metric_row("biological_weighted_pc_r2", weighted_pc_r2_explicit(matrix[, study$ids], factor(study$meta$sample), "categorical"), "higher", length(study$ids), "weighted R-squared"),
      evaluation_metric_row("feature_repeatability_ratio_median", stats::median(repeatability, na.rm = TRUE), "higher", sum(is.finite(repeatability)), "repeatability ratio"),
      evaluation_metric_row("sample_replicate_pearson_median", stats::median(replicates$pearson, na.rm = TRUE), "higher", nrow(replicates), "correlation"),
      evaluation_metric_row("genuine_replicate_icc_median", stats::median(replicates$icc_a1, na.rm = TRUE), "higher", nrow(replicates), "ICC(A,1)")
    )
  } else if (identical(dataset_key, "batchcorr_set1")) {
    replicates <- batchcorr_replicate_metrics(matrix, prepared$meta)
    repeatability <- batchcorr_feature_repeatability(matrix, prepared$meta)
    associations <- batchcorr_associations(matrix, prepared$meta, method_id)
    metrics <- rbind(metrics,
      evaluation_metric_row("biological_weighted_pc_r2", weighted_pc_r2_explicit(matrix[, study$ids], factor(study$meta$accession_id), "categorical"), "higher", length(study$ids), "weighted R-squared"),
      evaluation_metric_row("feature_repeatability_ratio_median", stats::median(repeatability, na.rm = TRUE), "higher", sum(is.finite(repeatability)), "repeatability ratio"),
      evaluation_metric_row("sample_replicate_pearson_median", stats::median(replicates$pearson, na.rm = TRUE), "higher", nrow(replicates), "correlation"),
      evaluation_metric_row("genuine_replicate_icc_median", stats::median(replicates$icc_a1, na.rm = TRUE), "higher", nrow(replicates), "ICC(A,1)"),
      evaluation_metric_row("biological_associated_features", sum(associations$significant), "context", nrow(associations), "features")
    )
  } else if (identical(dataset_key, "sacurine")) {
    biology <- sacurine_biology(matrix, prepared$meta, method_id)
    associations <- biology$primary
    metrics <- rbind(metrics,
      evaluation_metric_row("biological_weighted_pc_r2", mean(biology$pc$value), "higher", 3, "weighted R-squared"),
      evaluation_metric_row("biological_associated_features", sum(biology$primary$significant), "context", nrow(biology$primary), "feature-variable tests"),
      evaluation_metric_row("cross_batch_effect_pearson_median", stats::median(biology$concordance$pearson, na.rm = TRUE), "higher", nrow(biology$concordance), "correlation")
    )
  } else if (identical(dataset_key, "waveica")) {
    biology <- waveica_biology(matrix, prepared$meta, method_id)
    associations <- biology$primary
    metrics <- rbind(metrics,
      evaluation_metric_row("biological_weighted_pc_r2", biology$group_pc_r2, "higher", length(study$ids), "weighted R-squared"),
      evaluation_metric_row("biological_associated_features", sum(biology$primary$significant), "context", nrow(biology$primary), "features"),
      evaluation_metric_row("cross_batch_effect_pearson_median", stats::median(biology$concordance$pearson, na.rm = TRUE), "higher", nrow(biology$concordance), "correlation")
    )
  }

  metrics$dataset_key <- dataset_key
  metrics$method_id <- method_id
  metrics$method <- method_label
  qc_table <- data.frame(
    dataset_key = dataset_key, method_id = method_id, method = method_label,
    feature_id = names(qc), cv_percent = as.numeric(qc), stringsAsFactors = FALSE
  )
  gam$dataset_key <- dataset_key; gam$method_id <- method_id
  lb$dataset_key <- dataset_key; lb$method_id <- method_id
  if (nrow(associations)) {
    associations$dataset_key <- dataset_key; associations$method_id <- method_id
  }
  if (nrow(replicates)) {
    replicates$dataset_key <- dataset_key; replicates$method_id <- method_id
  }
  list(
    metrics = metrics, qc = qc_table, gam = gam, ljung_box = lb,
    associations = associations, replicates = replicates
  )
}
