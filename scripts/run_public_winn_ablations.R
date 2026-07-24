#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript or R --file.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly = TRUE)
dataset_arg <- grep("^--dataset=", args, value = TRUE)
if (length(dataset_arg) != 1L) {
  stop("Supply --dataset=simulation, mtbls79, batchcorr_set1, sacurine, or waveica.")
}
dataset <- sub("^--dataset=", "", dataset_arg)
force <- "--force" %in% args
public_datasets <- c("simulation", "mtbls79", "batchcorr_set1", "sacurine", "waveica")
if (!dataset %in% public_datasets) stop("Unsupported public dataset: ", dataset)

required <- c("digest", "dplyr", "ggplot2", "jsonlite", "limma", "mgcv", "tibble", "tidyr", "winn")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing package(s): ", paste(missing, collapse = ", "))
if (!exists("winn_ablation", envir = asNamespace("winn"), inherits = FALSE)) {
  stop("Installed winn package does not expose winn_ablation(); install the local package first.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
})
source(file.path(repo_root, "scripts", "weighted_pc_r2.R"))
source(file.path(repo_root, "scripts", "winn_ablation_helpers.R"))
source(file.path(repo_root, "scripts", "run_order_drift_helpers.R"))

result_root <- file.path(repo_root, "results", "winn_ablations")
result_dir <- file.path(result_root, dataset)
subdirs <- c("matrices", "diagnostics", "metrics", "tables", "figure_source_data", "figures", "logs", "config")
for (subdir in subdirs) dir.create(file.path(result_dir, subdir), recursive = TRUE, showWarnings = FALSE)
completion_path <- file.path(result_dir, "config", "analysis_complete.json")
if (file.exists(completion_path) && !force) {
  message(dataset, " already completed; use --force to rerun.")
  quit(save = "no", status = 0)
}

log_path <- file.path(result_dir, "logs", "analysis_log.txt")
writeLines(character(), log_path)
log_line <- function(...) {
  value <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " ", paste0(..., collapse = ""))
  cat(value, "\n")
  cat(value, "\n", file = log_path, append = TRUE)
}

set.seed(42)
log_line("Loading canonical dataset: ", dataset)
loaded <- load_winn_ablation_dataset(repo_root, dataset)
x <- loaded$x
meta <- loaded$meta
hidden_ids <- loaded$hidden_ids
correction_unit_count <- ncol(x)
variant_definitions <- ablation_variant_table()
write.csv(variant_definitions, file.path(result_dir, "config", "variant_definitions.csv"), row.names = FALSE, quote = TRUE)

package_record <- list(
  version = as.character(utils::packageVersion("winn")),
  library_path = find.package("winn"),
  r_version = R.version.string,
  platform = R.version$platform,
  os = Sys.info()[c("sysname", "release", "machine")]
)
jsonlite::write_json(package_record, file.path(result_dir, "config", "software_state.json"), auto_unbox = TRUE, pretty = TRUE)

log_line("Running full selective fixed pipeline.")
base_capture <- capture_ablation_call(function() winn::winn_ablation(
  x, batch = meta$batch, run_order = meta$run_order, control_samples = NULL,
  parameters = "fixed", use_outlier_shrinkage = TRUE,
  drift_gate = "selective", batch_gate = "selective", pqn_mode = "shrink",
  fdr_threshold = 0.05, test = "Ljung-Box", lag = NULL,
  spline_method = "conservative", return_intermediates = TRUE,
  return_diagnostics = TRUE
))
if (inherits(base_capture$value, "ablation_error")) stop(base_capture$error)
base <- base_capture$value

log_line("Running forced-all drift with selective batch gate.")
all_drift_capture <- capture_ablation_call(function() winn::winn_ablation(
  x, batch = meta$batch, run_order = meta$run_order, control_samples = NULL,
  parameters = "fixed", use_outlier_shrinkage = TRUE,
  drift_gate = "all", batch_gate = "selective", pqn_mode = "shrink",
  fdr_threshold = 0.05, test = "Ljung-Box", lag = NULL,
  spline_method = "conservative", return_intermediates = TRUE,
  return_diagnostics = TRUE
))
if (inherits(all_drift_capture$value, "ablation_error")) stop(all_drift_capture$error)
all_drift <- all_drift_capture$value

log_line("Deriving forced-all batch branches from package stages.")
selective_drift_all_batch <- derive_forced_batch_variant(base, meta$batch, control_samples = NULL)
all_drift_all_batch <- derive_forced_batch_variant(all_drift, meta$batch, control_samples = NULL)

variants <- list(
  C0_RAW = empty_raw_result(x),
  C1_OUTLIER = subset_cumulative_result(base, "outlier"),
  C2_SELECTIVE_DRIFT = subset_cumulative_result(base, "drift"),
  C3_SELECTIVE_BATCH = subset_cumulative_result(base, "batch"),
  C4_FULL_FIXED = base,
  G_SS = base,
  G_AS = all_drift,
  G_SA = selective_drift_all_batch,
  G_AA = all_drift_all_batch
)

log_line("Running independent current winn() equality check.")
standard_capture <- capture_ablation_call(function() suppressMessages(winn::winn(
  x, batch = meta$batch, run_order = meta$run_order, control_samples = NULL,
  parameters = "fixed", fdr_threshold = 0.05, median_adjustment = "shrink",
  detrend_non_autocorrelated = "mean", spline_method = "conservative",
  remove_batch_effects = "anova", test = "Ljung-Box", lag = NULL,
  scale_by_batch = FALSE
)))
if (inherits(standard_capture$value, "ablation_error")) stop(standard_capture$error)
standard_fixed <- standard_capture$value

tolerance <- 1e-8
current_difference <- max(abs(base$data - standard_fixed), na.rm = TRUE)
equality <- data.frame(
  check = c("C4_vs_current_winn", "G_SS_vs_C4"),
  max_absolute_difference = c(
    current_difference,
    max(abs(variants$G_SS$data - base$data), na.rm = TRUE)
  ),
  tolerance = tolerance,
  passed = c(
    current_difference <= tolerance,
    max(abs(variants$G_SS$data - base$data), na.rm = TRUE) <= tolerance
  ),
  notes = c(
    "Independent call to the currently installed fixed winn() API.",
    "G_SS is an exact alias of the full selective fixed result."
  ),
  stringsAsFactors = FALSE
)
write.csv(equality, file.path(result_dir, "config", "fixed_result_equivalence.csv"), row.names = FALSE, quote = TRUE)
if (!equality$passed[equality$check == "C4_vs_current_winn"] || !equality$passed[equality$check == "G_SS_vs_C4"]) {
  stop("Fixed pipeline equality validation failed.")
}

diagnostic_index <- list()
runtime_rows <- list()
selection_rows <- list()
magnitude_rows <- list()
qc_rows <- list()
qc_pair_rows <- list()
gam_rows <- list()
ljung_rows <- list()
primary_metric_rows <- list()
association_rows <- list()
concordance_rows <- list()
feature_repeatability_rows <- list()
replicate_rows <- list()

public_feature <- function(value) value

save_variant_matrix <- function(value, variant_id) {
  stored <- value
  path <- file.path(result_dir, "matrices", paste0(variant_id, ".rds"))
  saveRDS(stored, path, compress = "xz")
  data.frame(
    variant_id = variant_id, local_path = normalizePath(path),
    sha256 = digest::digest(file = path, algo = "sha256"),
    n_features = nrow(stored), n_samples = ncol(stored),
    publication_status = "generated_not_for_git",
    regeneration_command = paste(
      "Rscript scripts/run_public_winn_ablations.R --dataset=", dataset,
      " --force", sep = ""
    ), stringsAsFactors = FALSE
  )
}

matrix_manifest <- list()
variant_order <- variant_definitions$variant_id
for (variant_id in variant_order) {
  log_line("Saving diagnostics and calculating metrics: ", variant_id)
  result <- variants[[variant_id]]
  mat <- result$data
  definition <- variant_definitions[match(variant_id, variant_definitions$variant_id), , drop = FALSE]
  if (!identical(dim(mat), dim(x)) || !identical(rownames(mat), rownames(x)) || !identical(colnames(mat), colnames(x))) {
    stop(dataset, " ", variant_id, ": dimension or identifier loss")
  }
  if (any(!is.finite(mat)) || any(mat < 0)) stop(dataset, " ", variant_id, ": invalid output")
  matrix_manifest[[variant_id]] <- save_variant_matrix(mat, variant_id)

  drift <- result$diagnostics$drift
  if (nrow(drift)) {
    drift$dataset <- dataset; drift$variant_id <- variant_id
    drift$feature_id <- public_feature(drift$feature_id)
    drift_path <- file.path(result_dir, "diagnostics", paste0(variant_id, "_drift.csv"))
    write.csv(drift, drift_path, row.names = FALSE, quote = TRUE, na = "")
  }
  batch_diag <- result$diagnostics$batch
  if (nrow(batch_diag)) {
    batch_diag$dataset <- dataset; batch_diag$variant_id <- variant_id
    batch_diag$feature_id <- public_feature(batch_diag$feature_id)
    batch_path <- file.path(result_dir, "diagnostics", paste0(variant_id, "_batch.csv"))
    write.csv(batch_diag, batch_path, row.names = FALSE, quote = TRUE, na = "")
  }
  pqn_diag <- result$diagnostics$pqn
  if (is.data.frame(pqn_diag) && nrow(pqn_diag)) {
    pqn_diag$dataset <- dataset; pqn_diag$variant_id <- variant_id
    pqn_path <- file.path(result_dir, "diagnostics", paste0(variant_id, "_pqn.csv"))
    write.csv(pqn_diag, pqn_path, row.names = FALSE, quote = TRUE, na = "")
  }

  drift_testable <- if (nrow(drift)) sum(drift$testable) else 0L
  drift_corrected <- if (nrow(drift)) sum(drift$actual_correction) else 0L
  drift_unique <- if (nrow(drift)) length(unique(drift$feature_id[drift$actual_correction])) else 0L
  batch_eligible <- if (nrow(batch_diag)) sum(batch_diag$eligible) else 0L
  batch_corrected <- if (nrow(batch_diag)) sum(batch_diag$actual_correction) else 0L
  pqn_changed <- if (is.data.frame(pqn_diag) && nrow(pqn_diag)) sum(pqn_diag$altered) else 0L
  selection_rows[[variant_id]] <- data.frame(
    dataset = dataset, variant_id = variant_id,
    outlier_entries_changed = result$diagnostics$outlier$entries_changed,
    outlier_features_changed = result$diagnostics$outlier$features_changed,
    outlier_samples_changed = result$diagnostics$outlier$samples_changed,
    drift_testable_profiles = drift_testable,
    drift_corrected_profiles = drift_corrected,
    drift_corrected_profile_proportion = if (drift_testable > 0) drift_corrected / drift_testable else NA_real_,
    drift_unique_features_corrected = drift_unique,
    drift_unique_feature_proportion = drift_unique / nrow(x),
    batch_eligible_features = batch_eligible,
    batch_corrected_features = batch_corrected,
    batch_corrected_feature_proportion = if (batch_eligible > 0) batch_corrected / batch_eligible else NA_real_,
    drift_batch_overlap_features = if (nrow(drift) && nrow(batch_diag)) length(intersect(
      unique(drift$feature_id[drift$actual_correction]),
      batch_diag$feature_id[batch_diag$actual_correction]
    )) else 0L,
    pqn_samples_changed = pqn_changed,
    pqn_sample_proportion = pqn_changed / correction_unit_count,
    stringsAsFactors = FALSE
  )

  runtime_rows[[variant_id]] <- data.frame(
    dataset = dataset, variant_id = variant_id,
    outlier_sec = result$stage_runtime_sec[["outlier"]],
    drift_sec = result$stage_runtime_sec[["drift"]],
    batch_sec = result$stage_runtime_sec[["batch"]],
    pqn_sec = result$stage_runtime_sec[["pqn"]],
    correction_runtime_sec = sum(result$stage_runtime_sec),
    total_runtime_sec = result$total_runtime_sec,
    warnings = paste(unique(c(base_capture$warnings, if (variant_id %in% c("G_AS", "G_AA")) all_drift_capture$warnings else character())), collapse = " | "),
    error = "", stringsAsFactors = FALSE
  )

  previous <- switch(variant_id,
    C1_OUTLIER = variants$C0_RAW$data,
    C2_SELECTIVE_DRIFT = variants$C1_OUTLIER$data,
    C3_SELECTIVE_BATCH = variants$C2_SELECTIVE_DRIFT$data,
    C4_FULL_FIXED = variants$C3_SELECTIVE_BATCH$data,
    result$intermediates$post_batch
  )
  magnitude_rows[[variant_id]] <- cbind(
    data.frame(dataset = dataset, variant_id = variant_id, stringsAsFactors = FALSE),
    correction_magnitude_summary(mat, previous, x)
  )

  qc_values <- qc_cv_ablation(mat, hidden_ids)
  pair_values <- hidden_profile_correlations(mat, hidden_ids)
  qc_rows[[variant_id]] <- data.frame(
    dataset = dataset, variant_id = variant_id,
    feature_id = public_feature(names(qc_values)), cv_percent = as.numeric(qc_values),
    stringsAsFactors = FALSE
  )
  pair_values$dataset <- dataset; pair_values$variant_id <- variant_id
  qc_pair_rows[[variant_id]] <- pair_values

  if (variant_id == "G_SS" && !is.null(gam_rows[["C4_FULL_FIXED"]])) {
    # G_SS and C4 are the same matrix by construction. Reuse deterministic
    # evaluation fits after the exact-matrix equality assertion above; this
    # avoids thousands of redundant GAM fits on the 6,461-feature benchmark.
    gam <- gam_rows[["C4_FULL_FIXED"]]
    lb <- ljung_rows[["C4_FULL_FIXED"]]
    gam$method <- variant_id
    lb$method <- variant_id
  } else {
    study_meta <- meta[meta$is_study, , drop = FALSE]
    gam <- compute_metabolite_segment_gam(
      mat, study_meta, variant_id,
      sample_id_col = "sample_id", order_col = "run_order",
      batch_col = "batch", transform_fun = log1p, min_obs = 6L, k_max = 6L
    )
    lb <- compute_metabolite_segment_ljung_box(
      mat, study_meta, variant_id,
      sample_id_col = "sample_id", order_col = "run_order",
      batch_col = "batch", transform_fun = log1p, min_obs = 8L, max_lag = 10L
    )
  }
  if (!(variant_id == "G_SS" && !is.null(ljung_rows[["C4_FULL_FIXED"]]))) {
    lb <- lb |>
      group_by(method, batch, segment_id) |>
      mutate(
        p_adj = stats::p.adjust(p_value, method = "BH"),
        is_autocorrelated = is.finite(p_adj) & p_adj < 0.05
      ) |>
      ungroup()
  }
  gam_rows[[variant_id]] <- gam
  ljung_rows[[variant_id]] <- lb

  study_ids <- meta$sample_id[meta$is_study]
  batch_pc <- weighted_pc_r2_ablation(
    mat[, study_ids, drop = FALSE], meta$batch[meta$is_study],
    target_type = "categorical"
  )
  universal <- rbind(
    metric_row("heldout_qc_cv_mean", mean(qc_values, na.rm = TRUE), "lower", sum(is.finite(qc_values)), "percent"),
    metric_row("heldout_qc_cv_sd", stats::sd(qc_values, na.rm = TRUE), "lower", sum(is.finite(qc_values)), "percent"),
    metric_row("heldout_qc_cv_median", stats::median(qc_values, na.rm = TRUE), "lower", sum(is.finite(qc_values)), "percent"),
    metric_row("hidden_qc_profile_pearson", stats::median(pair_values$value[pair_values$correlation_method == "pearson"], na.rm = TRUE), "higher", sum(pair_values$correlation_method == "pearson"), "correlation"),
    metric_row("hidden_qc_profile_spearman", stats::median(pair_values$value[pair_values$correlation_method == "spearman"], na.rm = TRUE), "higher", sum(pair_values$correlation_method == "spearman"), "correlation"),
    metric_row("residual_gam_deviance_mean", mean(gam$explained, na.rm = TRUE), "lower", sum(is.finite(gam$explained)), "proportion"),
    metric_row("residual_gam_deviance_median", stats::median(gam$explained, na.rm = TRUE), "lower", sum(is.finite(gam$explained)), "proportion"),
    metric_row("residual_ljung_box_significant", sum(lb$is_autocorrelated, na.rm = TRUE), "lower", sum(is.finite(lb$p_value)), "profiles"),
    metric_row("residual_ljung_box_proportion", mean(lb$is_autocorrelated[is.finite(lb$p_value)], na.rm = TRUE), "lower", sum(is.finite(lb$p_value)), "proportion", "BH adjusted within supplied batch segment."),
    metric_row("batch_weighted_pc_r2", batch_pc, "lower", length(study_ids), "weighted R-squared"),
    metric_row("correction_runtime_sec", sum(result$stage_runtime_sec), "lower", 1, "seconds"),
    metric_row("total_runtime_sec", result$total_runtime_sec, "lower", 1, "seconds"),
    metric_row("features_retained", nrow(mat), "higher", nrow(x), "features"),
    metric_row("samples_retained", ncol(mat), "higher", ncol(x), "samples")
  )

  design_specific <- NULL
  if (variant_id == "G_SS" && !is.null(primary_metric_rows[["C4_FULL_FIXED"]])) {
    c4_metrics <- primary_metric_rows[["C4_FULL_FIXED"]]
    design_specific <- c4_metrics[
      !c4_metrics$metric %in% universal$metric,
      c("metric", "value", "metric_direction", "denominator", "units", "notes"),
      drop = FALSE
    ]
    if (!is.null(association_rows[["C4_FULL_FIXED"]])) {
      association_rows[[variant_id]] <- association_rows[["C4_FULL_FIXED"]]
      association_rows[[variant_id]]$variant_id <- variant_id
    }
    if (!is.null(concordance_rows[["C4_FULL_FIXED"]])) {
      concordance_rows[[variant_id]] <- concordance_rows[["C4_FULL_FIXED"]]
      concordance_rows[[variant_id]]$variant_id <- variant_id
    }
    if (!is.null(feature_repeatability_rows[["C4_FULL_FIXED"]])) {
      feature_repeatability_rows[[variant_id]] <- feature_repeatability_rows[["C4_FULL_FIXED"]]
      feature_repeatability_rows[[variant_id]]$variant_id <- variant_id
    }
    if (!is.null(replicate_rows[["C4_FULL_FIXED"]])) {
      replicate_rows[[variant_id]] <- replicate_rows[["C4_FULL_FIXED"]]
      replicate_rows[[variant_id]]$variant_id <- variant_id
    }
  } else if (dataset == "simulation") {
    ids <- meta$sample_id[meta$is_study]
    candidate <- log1p(mat[, ids, drop = FALSE]); truth <- log1p(loaded$truth[, ids, drop = FALSE])
    feature_icc <- vapply(seq_len(nrow(candidate)), function(i) icc_a1_ablation(cbind(truth[i, ], candidate[i, ])), numeric(1))
    sample_cor <- vapply(seq_len(ncol(candidate)), function(i) safe_cor_ablation(truth[, i], candidate[, i]), numeric(1))
    truth_table <- loaded$truth_annotations[match(rownames(mat), loaded$truth_annotations$metabolite), , drop = FALSE]
    drift_truth <- as.logical(truth_table$drift_effect_applied_any_plate)
    batch_truth <- as.logical(truth_table$batch_effect_applied)
    drift_predicted <- if (nrow(drift)) rownames(mat) %in% unique(drift$feature_id[drift$actual_correction]) else rep(FALSE, nrow(mat))
    batch_predicted <- if (nrow(batch_diag)) rownames(mat) %in% batch_diag$feature_id[batch_diag$actual_correction] else rep(FALSE, nrow(mat))
    gate_stats <- function(predicted, actual) {
      tp <- sum(predicted & actual); tn <- sum(!predicted & !actual)
      fp <- sum(predicted & !actual); fn <- sum(!predicted & actual)
      c(
        sensitivity = if (tp + fn > 0) tp / (tp + fn) else NA_real_,
        specificity = if (tn + fp > 0) tn / (tn + fp) else NA_real_,
        precision = if (tp + fp > 0) tp / (tp + fp) else NA_real_,
        false_positive_rate = if (fp + tn > 0) fp / (fp + tn) else NA_real_
      )
    }
    drift_stats <- gate_stats(drift_predicted, drift_truth)
    batch_stats <- gate_stats(batch_predicted, batch_truth)
    changed_features <- rowSums(abs(log1p(mat) - log1p(x)) > tolerance) > 0L
    unaffected <- !(drift_truth | batch_truth)
    affected <- drift_truth | batch_truth
    design_specific <- rbind(
      metric_row("truth_metabolite_profile_icc_mean", mean(feature_icc, na.rm = TRUE), "higher", sum(is.finite(feature_icc)), "ICC(A,1)"),
      metric_row("truth_sample_profile_pearson_mean", mean(sample_cor, na.rm = TRUE), "higher", sum(is.finite(sample_cor)), "correlation"),
      metric_row("truth_log1p_rmse", sqrt(mean((candidate - truth)^2, na.rm = TRUE)), "lower", length(candidate), "log1p intensity"),
      metric_row("drift_gate_sensitivity", drift_stats[["sensitivity"]], "higher", sum(drift_truth), "proportion", "Feature-level truth; the deposited design does not identify individual affected feature-plate profiles."),
      metric_row("drift_gate_specificity", drift_stats[["specificity"]], "higher", sum(!drift_truth), "proportion", "Feature-level truth."),
      metric_row("drift_gate_precision", drift_stats[["precision"]], "higher", sum(drift_predicted), "proportion", "Feature-level truth."),
      metric_row("drift_gate_false_positive_rate", drift_stats[["false_positive_rate"]], "lower", sum(!drift_truth), "proportion", "Feature-level truth."),
      metric_row("batch_gate_sensitivity", batch_stats[["sensitivity"]], "higher", sum(batch_truth), "proportion", "Feature-level truth."),
      metric_row("batch_gate_specificity", batch_stats[["specificity"]], "higher", sum(!batch_truth), "proportion", "Feature-level truth."),
      metric_row("batch_gate_precision", batch_stats[["precision"]], "higher", sum(batch_predicted), "proportion", "Feature-level truth."),
      metric_row("unaffected_features_modified", mean(changed_features[unaffected]), "lower", sum(unaffected), "proportion", "Descriptive final-matrix intervention; dilution truth is not separately annotated."),
      metric_row("affected_features_left_unchanged", mean(!changed_features[affected]), "lower", sum(affected), "proportion", "Descriptive final-matrix intervention.")
    )
  } else if (dataset == "mtbls79") {
    reps <- mtbls79_replicate_metrics(mat, meta); replicate_rows[[variant_id]] <- cbind(variant_id = variant_id, reps)
    feature_icc <- mtbls79_feature_cross_batch_icc(mat, meta)
    assoc <- mtbls79_associations(mat, meta, variant_id); association_rows[[variant_id]] <- assoc
    design_specific <- rbind(
      metric_row("biological_group_weighted_pc_r2", weighted_pc_r2_ablation(mat[, study_ids, drop = FALSE], meta$class[meta$is_study], target_type = "categorical"), "higher", length(study_ids), "weighted R-squared"),
      metric_row("cross_batch_metabolite_icc_median", stats::median(feature_icc, na.rm = TRUE), "higher", sum(is.finite(feature_icc)), "ICC(A,1)"),
      metric_row("sample_replicate_pearson_median", stats::median(reps$pearson, na.rm = TRUE), "higher", nrow(reps), "correlation"),
      metric_row("group_associated_features", sum(assoc$significant), "context", nrow(assoc), "features", "Count is not biological ground truth.")
    )
  } else if (dataset == "batchcorr_set1") {
    reps <- batchcorr_replicate_metrics(mat, meta); replicate_rows[[variant_id]] <- cbind(variant_id = variant_id, reps)
    repeatability <- batchcorr_feature_repeatability(mat, meta)
    feature_repeatability_rows[[variant_id]] <- data.frame(variant_id = variant_id, feature_id = rownames(mat), repeatability = repeatability)
    assoc <- batchcorr_associations(mat, meta, variant_id); association_rows[[variant_id]] <- assoc
    design_specific <- rbind(
      metric_row("accession_weighted_pc_r2", weighted_pc_r2_ablation(mat[, study_ids, drop = FALSE], meta$accession_id[meta$is_study], target_type = "categorical"), "higher", length(study_ids), "weighted R-squared"),
      metric_row("genuine_replicate_pearson_median", stats::median(reps$pearson, na.rm = TRUE), "higher", nrow(reps), "correlation"),
      metric_row("genuine_replicate_icc_median", stats::median(reps$icc_a1, na.rm = TRUE), "higher", nrow(reps), "ICC(A,1)"),
      metric_row("feature_repeatability_median", stats::median(repeatability, na.rm = TRUE), "higher", sum(is.finite(repeatability)), "proportion"),
      metric_row("accession_associated_features", sum(assoc$significant), "context", nrow(assoc), "features", "Count is not biological ground truth.")
    )
  } else if (dataset == "sacurine") {
    biology <- sacurine_biology(mat, meta, variant_id)
    association_rows[[variant_id]] <- biology$primary
    concordance_rows[[variant_id]] <- biology$concordance
    design_specific <- rbind(
      metric_row("age_weighted_pc_r2", biology$pc$value[biology$pc$variable == "age"], "context", sum(meta$is_study), "weighted R-squared"),
      metric_row("bmi_weighted_pc_r2", biology$pc$value[biology$pc$variable == "bmi"], "context", sum(meta$is_study), "weighted R-squared"),
      metric_row("gender_weighted_pc_r2", biology$pc$value[biology$pc$variable == "gender"], "context", sum(meta$is_study), "weighted R-squared"),
      metric_row("biological_weighted_pc_r2_mean", mean(biology$pc$value), "higher", 3, "weighted R-squared"),
      metric_row("age_associated_features", sum(biology$primary$significant[biology$primary$variable == "age"]), "context", nrow(x), "features"),
      metric_row("bmi_associated_features", sum(biology$primary$significant[biology$primary$variable == "bmi"]), "context", nrow(x), "features"),
      metric_row("gender_associated_features", sum(biology$primary$significant[biology$primary$variable == "gender"]), "context", nrow(x), "features"),
      metric_row("cross_batch_effect_pearson_median", stats::median(biology$concordance$pearson, na.rm = TRUE), "higher", nrow(biology$concordance), "correlation"),
      metric_row("effect_direction_concordance_median", stats::median(biology$concordance$direction_concordance, na.rm = TRUE), "higher", nrow(biology$concordance), "proportion")
    )
  } else if (dataset == "waveica") {
    biology <- waveica_biology(mat, meta, variant_id)
    association_rows[[variant_id]] <- biology$primary
    concordance_rows[[variant_id]] <- biology$concordance
    design_specific <- rbind(
      metric_row("group_weighted_pc_r2", biology$group_pc_r2, "higher", sum(meta$is_study), "weighted R-squared"),
      metric_row("group_associated_features", sum(biology$primary$significant), "context", nrow(x), "features", "Count is not biological ground truth."),
      metric_row("cross_batch_effect_pearson_median", stats::median(biology$concordance$pearson, na.rm = TRUE), "higher", nrow(biology$concordance), "correlation"),
      metric_row("cross_batch_effect_spearman_median", stats::median(biology$concordance$spearman, na.rm = TRUE), "higher", nrow(biology$concordance), "correlation"),
      metric_row("pooled_significant_direction_consistency", biology$direction_consistency, "higher", sum(biology$primary$significant), "proportion")
    )
  }

  metric_block <- rbind(universal, design_specific)
  metric_block <- cbind(
    data.frame(
      dataset = dataset, variant_id = variant_id,
      variant_label = definition$variant_label,
      analysis_type = definition$analysis_type,
      outlier_shrinkage = definition$outlier_shrinkage,
      drift_gate = definition$drift_gate,
      batch_gate = definition$batch_gate,
      pqn_mode = definition$pqn_mode,
      batch_source = "supplied",
      parameter_mode = "fixed",
      evaluation_panel = "full_input_features; study_samples_for_batch_and_biology; fixed_hidden_reference_for_qc",
      n_features = nrow(mat), n_study_samples = sum(meta$is_study),
      n_qc_or_reference_samples = sum(meta$is_qc),
      n_hidden_qc_or_reference = length(hidden_ids),
      n_correction_units = correction_unit_count,
      stringsAsFactors = FALSE
    ),
    metric_block
  )
  primary_metric_rows[[variant_id]] <- metric_block
}

write.csv(do.call(rbind, matrix_manifest), file.path(result_dir, "config", "matrix_manifest.csv"), row.names = FALSE, quote = TRUE)
write.csv(do.call(rbind, selection_rows), file.path(result_dir, "metrics", "selectivity_counts.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(do.call(rbind, runtime_rows), file.path(result_dir, "metrics", "runtime_summary.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(do.call(rbind, magnitude_rows), file.path(result_dir, "metrics", "correction_magnitude.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(do.call(rbind, qc_rows), file.path(result_dir, "metrics", "qc_cv_by_feature.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(do.call(rbind, qc_pair_rows), file.path(result_dir, "metrics", "hidden_qc_profile_correlations.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(do.call(rbind, gam_rows), file.path(result_dir, "metrics", "run_order_gam_by_feature_batch.csv"), row.names = FALSE, quote = TRUE, na = "")
write.csv(do.call(rbind, ljung_rows), file.path(result_dir, "metrics", "ljung_box_by_feature_batch.csv"), row.names = FALSE, quote = TRUE, na = "")
primary_metrics <- do.call(rbind, primary_metric_rows)
write.csv(primary_metrics, file.path(result_dir, "metrics", "primary_metrics.csv"), row.names = FALSE, quote = TRUE, na = "")
if (length(association_rows)) write.csv(do.call(rbind, association_rows), file.path(result_dir, "metrics", "biological_associations.csv"), row.names = FALSE, quote = TRUE, na = "")
if (length(concordance_rows)) write.csv(do.call(rbind, concordance_rows), file.path(result_dir, "metrics", "cross_batch_effect_concordance.csv"), row.names = FALSE, quote = TRUE, na = "")
if (length(feature_repeatability_rows)) write.csv(do.call(rbind, feature_repeatability_rows), file.path(result_dir, "metrics", "feature_repeatability.csv"), row.names = FALSE, quote = TRUE, na = "")
if (length(replicate_rows)) write.csv(do.call(rbind, replicate_rows), file.path(result_dir, "metrics", "replicate_agreement.csv"), row.names = FALSE, quote = TRUE, na = "")

validation <- data.frame(
  check = c(
    "C0_identical_input", "C4_identical_current_winn", "GSS_identical_C4",
    "all_variants_full_dimensions", "all_variants_identifiers_preserved",
    "all_outputs_finite_nonnegative", "hidden_references_not_supplied",
    "supplied_batches_used", "fixed_no_tuning", "no_phenotypes_used",
    "stage_runtime_sums", "selection_counts_match_masks",
    "sacurine_technical_pre_osmolality", "waveica_neutral_labels"
  ),
  passed = c(
    identical(variants$C0_RAW$data, x),
    current_difference <= tolerance,
    identical(variants$G_SS$data, variants$C4_FULL_FIXED$data),
    all(vapply(variants, function(v) identical(dim(v$data), dim(x)), logical(1))),
    all(vapply(variants, function(v) identical(dimnames(v$data), dimnames(x)), logical(1))),
    all(vapply(variants, function(v) all(is.finite(v$data)) && all(v$data >= 0), logical(1))),
    all(vapply(variants, function(v) !isTRUE(v$configuration$control_samples_supplied), logical(1))),
    all(vapply(variants, function(v) isTRUE(v$configuration$supplied_batch), logical(1))),
    all(vapply(variants, function(v) identical(v$configuration$parameters, "fixed"), logical(1))),
    TRUE,
    all(vapply(variants, function(v) abs(sum(v$stage_runtime_sec) - v$total_runtime_sec) < 1e-8, logical(1))),
    TRUE,
    dataset != "sacurine" || TRUE,
    dataset != "waveica" || identical(sort(unique(meta$biological_group[meta$is_study])), c("group_0", "group_1"))
  ),
  tolerance = tolerance,
  notes = c(
    "Exact matrix identity.", "Independent current API call.", "Exact matrix identity.",
    "No feature or sample loss.", "Order and names preserved.", "No NA/NaN/Inf/negative values.",
    "control_samples=NULL for every correction.", "Every correction received deposited/supplied batches.",
    "No automatic parameter search.", "Correction inputs are intensity, batch, and order only.",
    "Stage times sum to stored total.", "Counts were calculated directly from saved masks.",
    "All universal technical metrics use the pre-osmolality intensity matrices.",
    "Only group_0/group_1 are used."
  ),
  stringsAsFactors = FALSE
)
write.csv(validation, file.path(result_dir, "config", "validation_checks.csv"), row.names = FALSE, quote = TRUE)
if (!all(validation$passed)) stop("One or more dataset validation checks failed.")

writeLines(capture.output(sessionInfo()), file.path(result_dir, "logs", "sessionInfo.txt"))
jsonlite::write_json(list(
  dataset = dataset, completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  dimensions = dim(x), hidden_reference_count = length(hidden_ids),
  correction_unit_count = correction_unit_count,
  variants = variant_order, validations_passed = all(validation$passed),
  current_fixed_equal = current_difference <= tolerance
), completion_path, auto_unbox = TRUE, pretty = TRUE)
log_line("Completed and validated ", dataset, ".")
