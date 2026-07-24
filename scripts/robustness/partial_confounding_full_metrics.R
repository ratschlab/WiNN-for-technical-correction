# Full evaluation engine for one partial-confounding scenario.
#
# Correction methods must have finished before this file is called. Phenotype
# allocation and responsive-feature truth are evaluation-only inputs here.

partial_full_write_csv_gz <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- gzfile(path, open = "wt", compression = 9L)
  on.exit(close(connection), add = TRUE)
  utils::write.csv(value, connection, row.names = FALSE, quote = TRUE, na = "")
  invisible(path)
}

partial_full_metric_dictionary <- function() {
  metrics <- c(
    "coefficient_pearson_vs_injected_all",
    "coefficient_spearman_vs_injected_all",
    "coefficient_pearson_vs_injected_responsive",
    "coefficient_pearson_vs_clean_all",
    "effect_rmse_vs_injected_all",
    "effect_rmse_vs_injected_responsive",
    "effect_rmse_vs_clean_all",
    "effect_rmse_vs_clean_responsive",
    "recovered_vs_injected_slope",
    "recovered_vs_injected_intercept",
    "median_attenuation_ratio_injected_responsive",
    "median_attenuation_ratio_clean_responsive",
    "sign_concordance_responsive",
    "n_discoveries_bh_0p05",
    "true_positive_rate",
    "empirical_fdr",
    "null_feature_false_positive_rate",
    "clean_log1p_rmse_study",
    "clean_sample_profile_pearson_mean",
    "clean_sample_profile_pearson_median",
    "clean_feature_profile_pearson_mean",
    "clean_feature_profile_pearson_median",
    "clean_feature_icc_a1_mean",
    "clean_feature_icc_a1_median",
    "heldout_qc_cv_mean",
    "heldout_qc_cv_median",
    "heldout_qc_cv_sd",
    "batch_weighted_pc_r2_categorical",
    "phenotype_weighted_pc_r2_categorical",
    "residual_ljung_box_profiles_tested",
    "residual_ljung_box_profiles_significant",
    "residual_ljung_box_proportion_significant",
    "residual_ljung_box_median_statistic",
    "residual_run_order_gam_mean_deviance",
    "retained_features",
    "retained_samples",
    "runtime_sec"
  )
  higher <- c(
    "coefficient_pearson_vs_injected_all",
    "coefficient_spearman_vs_injected_all",
    "coefficient_pearson_vs_injected_responsive",
    "coefficient_pearson_vs_clean_all",
    "median_attenuation_ratio_injected_responsive",
    "median_attenuation_ratio_clean_responsive",
    "sign_concordance_responsive", "true_positive_rate",
    "clean_sample_profile_pearson_mean",
    "clean_sample_profile_pearson_median",
    "clean_feature_profile_pearson_mean",
    "clean_feature_profile_pearson_median",
    "clean_feature_icc_a1_mean", "clean_feature_icc_a1_median",
    "retained_features", "retained_samples"
  )
  lower <- c(
    "effect_rmse_vs_injected_all", "effect_rmse_vs_injected_responsive",
    "effect_rmse_vs_clean_all", "effect_rmse_vs_clean_responsive",
    "empirical_fdr", "null_feature_false_positive_rate",
    "clean_log1p_rmse_study", "heldout_qc_cv_mean",
    "heldout_qc_cv_median", "heldout_qc_cv_sd",
    "batch_weighted_pc_r2_categorical",
    "residual_ljung_box_profiles_significant",
    "residual_ljung_box_proportion_significant",
    "residual_ljung_box_median_statistic",
    "residual_run_order_gam_mean_deviance", "runtime_sec"
  )
  neutral <- setdiff(metrics, c(higher, lower))
  direction <- ifelse(metrics %in% higher, "higher",
                      ifelse(metrics %in% lower, "lower", "descriptive"))
  units <- rep("value", length(metrics))
  units[grepl("cv_", metrics)] <- "percent"
  units[grepl("correlation|pearson|spearman|concordance", metrics)] <-
    "correlation_or_proportion"
  units[grepl("r2", metrics)] <- "weighted_R2"
  units[grepl("rmse|effect|slope|intercept|attenuation", metrics)] <-
    "log1p_effect_or_ratio"
  units[grepl("proportion|rate|fdr", metrics)] <- "proportion"
  units[grepl("profiles_tested|profiles_significant|discoveries|retained", metrics)] <-
    "count"
  units[metrics == "runtime_sec"] <- "seconds"
  data.frame(
    metric = metrics,
    metric_direction = direction,
    units = units,
    primary_preservation_endpoint = metrics %in% c(
      "median_attenuation_ratio_clean_responsive",
      "effect_rmse_vs_clean_responsive", "true_positive_rate"
    ),
    injected_scale_secondary = grepl("injected", metrics) |
      metrics %in% c("recovered_vs_injected_slope",
                     "recovered_vs_injected_intercept"),
    stringsAsFactors = FALSE
  )
}

partial_full_metric_rows <- function(
  scenario_id, seed_id, method, table, dictionary
) {
  if (!all(c("metric", "value", "status") %in% names(table))) {
    stop("Metric table lacks metric/value/status columns.", call. = FALSE)
  }
  if (!setequal(table$metric, dictionary$metric) ||
      anyDuplicated(table$metric)) {
    stop(method, " produced an incomplete or duplicate metric set.",
         call. = FALSE)
  }
  table <- table[match(dictionary$metric, table$metric), , drop = FALSE]
  data.frame(
    scenario_id = scenario_id,
    seed_id = seed_id,
    method = method,
    metric = table$metric,
    value = as.numeric(table$value),
    metric_direction = dictionary$metric_direction,
    units = dictionary$units,
    status = table$status,
    stringsAsFactors = FALSE
  )
}

compute_partial_confounding_full_metrics <- function(
  scenario_id,
  seed_id,
  method_matrices,
  raw,
  clean,
  allocation,
  metadata_evaluation,
  hidden_ids,
  feature_truth,
  runtime,
  ablation_diagnostics,
  metric_dictionary,
  run_dir,
  log_line = function(...) invisible(NULL)
) {
  metric_dir <- file.path(run_dir, "metrics")
  dir.create(metric_dir, recursive = TRUE, showWarnings = FALSE)
  method_order <- names(method_matrices)
  study_ids <- allocation$sample_id[allocation$sample_type == "sample"]
  study_metadata <- data.frame(
    sample_id = study_ids,
    batch = allocation$plate[match(study_ids, allocation$sample_id)],
    run_order = allocation$run_order[match(study_ids, allocation$sample_id)],
    stringsAsFactors = FALSE
  )
  phenotype <- allocation$phenotype[match(study_ids, allocation$sample_id)]
  if (anyNA(study_metadata) || anyNA(phenotype) || length(study_ids) != 450L) {
    stop("Study metadata or phenotype is incomplete.", call. = FALSE)
  }

  clean_effect <- fit_partial_confounding_effects(
    clean, allocation, feature_truth
  )
  clean_effect_table <- clean_effect$estimates
  clean_effect_table$scenario_id <- scenario_id
  clean_effect_table$seed_id <- seed_id

  alias_source <- c(
    C0_RAW = "Raw",
    C4_FULL_FIXED = "WINN default (no QC)",
    G_SS = "WINN default (no QC)"
  )
  calculation_cache <- list()
  metric_rows <- list()
  effect_rows <- list()
  effect_design_rows <- list()
  qc_rows <- list()
  gam_rows <- list()
  ljung_rows <- list()
  clean_feature_rows <- list()
  clean_sample_rows <- list()
  intervention_rows <- list()

  max_difference <- function(a, b) max(abs(a - b), na.rm = TRUE)
  for (method in method_order) {
    log_line("Computing full metrics for ", method, ".")
    matrix <- method_matrices[[method]]
    alias <- if (method %in% names(alias_source)) {
      alias_source[[method]]
    } else {
      NA_character_
    }
    use_alias <- !is.na(alias) && alias %in% names(calculation_cache) &&
      max_difference(matrix, method_matrices[[alias]]) <= 1e-8

    if (use_alias) {
      calculated <- calculation_cache[[alias]]
    } else {
      effect <- fit_partial_confounding_effects(
        matrix, allocation, feature_truth
      )
      recovery <- partial_clean_matrix_recovery(matrix, clean, study_ids)
      qc_cv <- qc_cv_ablation(matrix, hidden_ids)
      gam <- compute_metabolite_segment_gam(
        matrix[, study_ids, drop = FALSE], study_metadata, method,
        sample_id_col = "sample_id", order_col = "run_order",
        batch_col = "batch", transform_fun = log1p,
        min_obs = 6L, k_max = 6L
      )
      ljung <- compute_metabolite_segment_ljung_box(
        matrix[, study_ids, drop = FALSE], study_metadata, method,
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
      calculated <- list(
        effect = effect,
        recovery = recovery,
        qc_cv = qc_cv,
        gam = gam,
        ljung = ljung,
        batch_r2 = weighted_pc_r2_explicit(
          matrix[, study_ids, drop = FALSE], study_metadata$batch,
          target_type = "categorical"
        ),
        phenotype_r2 = weighted_pc_r2_explicit(
          matrix[, study_ids, drop = FALSE], phenotype,
          target_type = "categorical"
        )
      )
    }
    calculation_cache[[method]] <- calculated

    effect_table <- calculated$effect$estimates
    effect_table$scenario_id <- scenario_id
    effect_table$seed_id <- seed_id
    effect_table$method <- method
    effect_rows[[method]] <- effect_table
    effect_design <- calculated$effect$design
    effect_design$scenario_id <- scenario_id
    effect_design$seed_id <- seed_id
    effect_design$method <- method
    effect_design_rows[[method]] <- effect_design

    qc_table <- data.frame(
      scenario_id = scenario_id, seed_id = seed_id, method = method,
      feature_id = names(calculated$qc_cv),
      cv_percent = as.numeric(calculated$qc_cv),
      stringsAsFactors = FALSE
    )
    qc_rows[[method]] <- qc_table

    gam <- calculated$gam
    gam$method <- method
    gam$scenario_id <- scenario_id
    gam$seed_id <- seed_id
    gam_rows[[method]] <- gam
    ljung <- calculated$ljung
    ljung$method <- method
    ljung$scenario_id <- scenario_id
    ljung$seed_id <- seed_id
    ljung_rows[[method]] <- ljung

    recovery_feature <- calculated$recovery$by_feature
    recovery_feature$scenario_id <- scenario_id
    recovery_feature$seed_id <- seed_id
    recovery_feature$method <- method
    clean_feature_rows[[method]] <- recovery_feature
    recovery_sample <- calculated$recovery$by_sample
    recovery_sample$scenario_id <- scenario_id
    recovery_sample$seed_id <- seed_id
    recovery_sample$method <- method
    clean_sample_rows[[method]] <- recovery_sample

    n_ljung <- sum(is.finite(ljung$p_value))
    n_gam <- sum(is.finite(gam$explained))
    technical <- data.frame(
      metric = c(
        "heldout_qc_cv_mean", "heldout_qc_cv_median",
        "heldout_qc_cv_sd", "batch_weighted_pc_r2_categorical",
        "phenotype_weighted_pc_r2_categorical",
        "residual_ljung_box_profiles_tested",
        "residual_ljung_box_profiles_significant",
        "residual_ljung_box_proportion_significant",
        "residual_ljung_box_median_statistic",
        "residual_run_order_gam_mean_deviance",
        "retained_features", "retained_samples", "runtime_sec"
      ),
      value = c(
        mean(calculated$qc_cv, na.rm = TRUE),
        stats::median(calculated$qc_cv, na.rm = TRUE),
        stats::sd(calculated$qc_cv, na.rm = TRUE),
        calculated$batch_r2, calculated$phenotype_r2,
        n_ljung, sum(ljung$is_autocorrelated, na.rm = TRUE),
        if (n_ljung) mean(ljung$is_autocorrelated[
          is.finite(ljung$p_value)
        ]) else NA_real_,
        stats::median(ljung$lb_stat, na.rm = TRUE),
        if (n_gam) mean(gam$explained, na.rm = TRUE) else NA_real_,
        nrow(matrix), ncol(matrix),
        runtime$runtime_sec[match(method, runtime$method)]
      ),
      status = c(
        rep("calculated", 5L),
        rep("calculated_full_feature_batch_profiles", 5L),
        rep("calculated", 3L)
      ),
      stringsAsFactors = FALSE
    )
    effect_metrics <- summarise_partial_effect_recovery(
      calculated$effect$estimates, clean_effect$estimates
    )
    all_metrics <- dplyr::bind_rows(
      effect_metrics, calculated$recovery$summary, technical
    )
    metric_rows[[method]] <- partial_full_metric_rows(
      scenario_id, seed_id, method, all_metrics, metric_dictionary
    )

    if (method %in% names(ablation_diagnostics)) {
      intervention <- summarise_partial_winn_interventions(
        ablation_diagnostics[[method]]
      )
      intervention$scenario_id <- scenario_id
      intervention$seed_id <- seed_id
      intervention$method <- method
      intervention_rows[[method]] <- intervention
    }
  }

  metrics <- dplyr::bind_rows(metric_rows)
  effects <- dplyr::bind_rows(effect_rows)
  designs <- dplyr::bind_rows(effect_design_rows)
  qc <- dplyr::bind_rows(qc_rows)
  gam <- dplyr::bind_rows(gam_rows)
  ljung <- dplyr::bind_rows(ljung_rows)
  clean_feature <- dplyr::bind_rows(clean_feature_rows)
  clean_sample <- dplyr::bind_rows(clean_sample_rows)
  interventions <- if (length(intervention_rows)) {
    dplyr::bind_rows(intervention_rows)
  } else {
    data.frame(
      scenario_id = character(), seed_id = character(), method = character(),
      stringsAsFactors = FALSE
    )
  }

  expected_metric_rows <- length(method_order) * nrow(metric_dictionary)
  if (nrow(metrics) != expected_metric_rows ||
      !identical(unique(metrics$method), method_order)) {
    stop("Full partial-confounding metric table is incomplete.",
         call. = FALSE)
  }
  utils::write.csv(metrics, file.path(run_dir, "method_metrics.csv"),
                   row.names = FALSE, quote = TRUE, na = "")
  utils::write.csv(clean_effect_table,
                   file.path(metric_dir, "clean_effect_reference.csv"),
                   row.names = FALSE, quote = TRUE, na = "")
  utils::write.csv(designs,
                   file.path(metric_dir, "effect_design_diagnostics.csv"),
                   row.names = FALSE, quote = TRUE, na = "")
  utils::write.csv(interventions,
                   file.path(metric_dir, "winn_intervention_summary.csv"),
                   row.names = FALSE, quote = TRUE, na = "")
  partial_full_write_csv_gz(effects,
    file.path(metric_dir, "phenotype_effect_estimates.csv.gz"))
  partial_full_write_csv_gz(qc,
    file.path(metric_dir, "heldout_qc_cv_by_feature.csv.gz"))
  partial_full_write_csv_gz(gam,
    file.path(metric_dir, "run_order_gam_by_feature_batch.csv.gz"))
  partial_full_write_csv_gz(ljung,
    file.path(metric_dir, "ljung_box_by_feature_batch.csv.gz"))
  partial_full_write_csv_gz(clean_feature,
    file.path(metric_dir, "clean_recovery_by_feature.csv.gz"))
  partial_full_write_csv_gz(clean_sample,
    file.path(metric_dir, "clean_recovery_by_sample.csv.gz"))

  list(
    method_metrics = metrics,
    effect_designs = designs,
    interventions = interventions,
    detailed_row_counts = data.frame(
      table = c(
        "phenotype_effect_estimates", "heldout_qc_cv_by_feature",
        "run_order_gam_by_feature_batch", "ljung_box_by_feature_batch",
        "clean_recovery_by_feature", "clean_recovery_by_sample"
      ),
      rows = c(nrow(effects), nrow(qc), nrow(gam), nrow(ljung),
               nrow(clean_feature), nrow(clean_sample)),
      stringsAsFactors = FALSE
    )
  )
}
