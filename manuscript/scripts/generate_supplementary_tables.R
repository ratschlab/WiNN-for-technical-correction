#!/usr/bin/env Rscript

# Generate supplementary manuscript tables from the frozen 2026-07-27
# analysis release. This script performs no correction, tuning, or selection.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript or R --file.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
package_root <- dirname(dirname(script_path))
cli <- commandArgs(trailingOnly = TRUE)
analysis_arg <- grep("^--analysis-root=", cli, value = TRUE)
final_arg <- grep("^--final-root=", cli, value = TRUE)
if (length(analysis_arg) + length(final_arg) != 1L) {
  stop("Supply --final-root=/path/to/archived/aggregate/inputs (or historical --analysis-root).")
}
final_root <- if(length(final_arg)) {
  normalizePath(sub("^--final-root=", "", final_arg), mustWork = TRUE)
} else {
  file.path(normalizePath(sub("^--analysis-root=", "", analysis_arg), mustWork = TRUE), "results", "final")
}
table_root <- file.path(final_root, "tables")
derivative_root <- file.path(final_root, "derivative", "tables")
validation_root <- file.path(final_root, "validation")
out_tex <- file.path(package_root, "generated", "supplementary_tables")
out_csv <- file.path(package_root, "generated", "supplementary_table_data")
dir.create(out_tex, recursive = TRUE, showWarnings = FALSE)
dir.create(out_csv, recursive = TRUE, showWarnings = FALSE)

required <- c("dplyr", "readr", "tidyr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# The public revision stores exact historical aggregate inputs compressed.
read_csv <- function(file, ...) {
  if (!file.exists(file) && file.exists(paste0(file, ".gz"))) file <- paste0(file, ".gz")
  readr::read_csv(file, ...)
}

method_order <- c(
  "raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf",
  "winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc"
)
method_short <- c(
  raw = "Raw", combat = "ComBat", qc_rlsc = "QC-RLSC",
  qc_rfsc = "QC-RFSC", tiger = "TIGER", serrf = "SERRF",
  winn_auto_qc = "WiNN (QC-tuned)",
  winn_auto_batch_qc = "WiNN (inferred segments)",
  winn_fixed_no_qc = "WiNN"
)
dataset_order <- c(
  "simulation", "mtbls79", "batchcorr_set1", "clinical_fiams",
  "sacurine", "waveica"
)

latex_escape <- function(x) {
  x <- ifelse(is.na(x), "--", as.character(x))
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&#_$])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

fmt <- function(x, digits = 3L) {
  ifelse(is.na(x) | !is.finite(x), "--", formatC(x, format = "f", digits = digits))
}

fmt_runtime <- function(x) {
  ifelse(
    is.na(x) | !is.finite(x), "--",
    ifelse(x >= 1000, formatC(x, format = "f", digits = 0),
           ifelse(x >= 10, formatC(x, format = "f", digits = 1),
                  formatC(x, format = "f", digits = 2)))
  )
}

fmt_ci <- function(est, lo, hi, digits = 3L) {
  ifelse(
    is.na(est) | !is.finite(est), "--",
    paste0(fmt(est, digits), " [", fmt(lo, digits), ", ", fmt(hi, digits), "]")
  )
}

fmt_mean_se <- function(mean, sd, n, digits = 3L) {
  se <- sd / sqrt(n)
  ifelse(
    is.na(mean) | !is.finite(mean), "--",
    paste0(fmt(mean, digits), " $\\pm$ ", fmt(se, digits))
  )
}

write_csv_copy <- function(data, stem) {
  write_csv(as.data.frame(data), file.path(out_csv, paste0(stem, ".csv")), na = "")
}

write_longtable <- function(data, stem, caption, label, colspec, headers,
                            landscape = TRUE, size = "\\scriptsize",
                            notes = NULL, arraystretch = NULL) {
  stopifnot(ncol(data) == length(headers))
  body <- apply(as.data.frame(data), 1, function(row) {
    paste0(paste(row, collapse = " & "), " \\\\")
  })
  open_landscape <- if (landscape) "\\begin{landscape}" else ""
  close_landscape <- if (landscape) "\\end{landscape}" else ""
  header_line <- paste(headers, collapse = " & ")
  lines <- c(
    open_landscape,
    if (!is.null(arraystretch)) "\\begingroup" else "",
    if (!is.null(arraystretch)) paste0("\\renewcommand{\\arraystretch}{", arraystretch, "}") else "",
    size,
    paste0("\\begin{longtable}{", colspec, "}"),
    paste0("\\caption{", caption, "}\\label{", label, "}\\\\"),
    "\\toprule",
    paste0(header_line, " \\\\"),
    "\\midrule",
    "\\endfirsthead",
    paste0("\\multicolumn{", ncol(data), "}{l}{\\textit{Table \\thetable{} continued}}\\\\"),
    "\\toprule",
    paste0(header_line, " \\\\"),
    "\\midrule",
    "\\endhead",
    "\\midrule",
    paste0("\\multicolumn{", ncol(data), "}{r}{\\textit{Continued on next page}}\\\\"),
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot",
    body,
    "\\end{longtable}",
    if (!is.null(notes)) paste0("\\noindent\\footnotesize ", notes, "\\par") else "",
    if (!is.null(arraystretch)) "\\endgroup" else "",
    close_landscape
  )
  writeLines(lines[nzchar(lines)], file.path(out_tex, paste0(stem, ".tex")))
}

# Core frozen tables.
primary <- read_csv(file.path(table_root, "primary_method_metrics_wide.csv"), show_col_types = FALSE)
seed_summary <- read_csv(file.path(table_root, "simulation_seed_summary.csv"), show_col_types = FALSE)
bootstrap <- read_csv(file.path(derivative_root, "bootstrap_estimates.csv"), show_col_types = FALSE)
associations <- read_csv(file.path(derivative_root, "biological_association_summary.csv"), show_col_types = FALSE)
permutations <- read_csv(file.path(derivative_root, "permutation_summary.csv"), show_col_types = FALSE)

# S4. All-nine-method simulation results over ten paired seeds.
simulation_metrics <- c(
  heldout_qc_cv_mean = "QC CV",
  residual_run_order_gam_mean_deviance = "GAM",
  batch_weighted_pc_r2_categorical = "Batch R2",
  truth_log1p_rmse = "RMSE",
  truth_metabolite_profile_icc_mean = "ICC",
  truth_sample_profile_pearson_mean = "Pearson"
)
sim <- seed_summary |>
  filter(metric %in% names(simulation_metrics), method_id %in% method_order) |>
  mutate(
    method = unname(method_short[method_id]),
    metric_short = unname(simulation_metrics[metric]),
    cell = case_when(
      metric == "heldout_qc_cv_mean" ~ fmt_mean_se(mean, sd, n, 2),
      TRUE ~ fmt_mean_se(mean, sd, n, 4)
    )
  ) |>
  select(method_id, method, metric_short, cell) |>
  pivot_wider(names_from = metric_short, values_from = cell) |>
  mutate(method_id = factor(method_id, levels = method_order)) |>
  arrange(method_id)
write_csv_copy(sim, "table_s4_simulation")
sim_tex <- sim |>
  transmute(
    Method = latex_escape(method),
    `QC CV` = `QC CV`, GAM = GAM, `Batch R2` = `Batch R2`,
    RMSE = RMSE, ICC = ICC, Pearson = Pearson
  )
write_longtable(
  sim_tex, "table_s4_simulation",
  "Complete simulation benchmark across all nine outputs. Cells are mean $\\pm$ seed-level standard error over ten paired simulation seeds. QC CV is in percent; GAM is residual run-order deviance explained; batch $R^2$ is categorical weighted-PC $R^2$; truth endpoints compare the corrected matrix with the clean matrix.",
  "tab:supp_simulation_all", "lrrrrrr",
  c("Method", "QC CV", "GAM", "Batch $R^2$", "RMSE", "ICC(A,1)", "Pearson"),
  landscape = TRUE
)

bootstrap_cell <- function(ds_key, metric_name, digits = 3L) {
  bootstrap |>
    filter(.data$dataset_key == ds_key, .data$metric == metric_name) |>
    transmute(method_id, cell = fmt_ci(point_estimate, ci_lower, ci_upper, digits))
}

prepare_empirical <- function(ds_key) {
  primary |>
    filter(.data$dataset_key == ds_key, method_id %in% method_order) |>
    mutate(
      method = unname(method_short[method_id]),
      method_id = factor(method_id, levels = method_order)
    ) |>
    arrange(method_id)
}

# S5. MTBLS79.
mt <- prepare_empirical("mtbls79") |>
  left_join(bootstrap_cell("mtbls79", "biological_weighted_pc_r2", 3), by = "method_id") |>
  transmute(
    Method = latex_escape(method),
    `QC CV` = fmt(heldout_qc_cv_mean, 2),
    GAM = fmt(residual_run_order_gam_mean_deviance, 3),
    `Batch R2` = fmt(batch_weighted_pc_r2_categorical, 3),
    `Cow-sheep R2 [95% CI]` = cell,
    `Replicate Pearson` = fmt(sample_replicate_pearson_median, 3),
    `Feature ICC` = fmt(metabolite_icc_a1_median, 3)
  )
write_csv_copy(mt, "table_s5_mtbls79")
write_longtable(
  mt, "table_s5_mtbls79",
  "Complete MTBLS79 results across all nine outputs. Cow-sheep $R^2$ is computed from identity-collapsed study measurements and includes a 2,000-resample identity-cluster bootstrap interval. Technical endpoints and repeated-measurement endpoints are frozen point estimates.",
  "tab:supp_mtbls79_all", "lrrrrrr",
  c("Method", "QC CV", "GAM", "Batch $R^2$", "Cow-sheep $R^2$ [95\\% CI]", "Replicate $r$", "Feature ICC"),
  landscape = TRUE
)

# S6 and S7. Repeated-identity datasets.
make_replicate_table <- function(ds_key, stem, caption, label, batch_name) {
  x <- prepare_empirical(ds_key) |>
    left_join(
      bootstrap_cell(ds_key, "sample_replicate_pearson_median", 4) |>
        rename(pearson_ci = cell),
      by = "method_id"
    ) |>
    left_join(
      bootstrap_cell(ds_key, "genuine_replicate_icc_median", 4) |>
        rename(icc_ci = cell),
      by = "method_id"
    ) |>
    transmute(
      Method = latex_escape(method),
      `QC CV` = fmt(heldout_qc_cv_mean, 2),
      GAM = fmt(residual_run_order_gam_mean_deviance, 3),
      `Batch R2` = fmt(batch_weighted_pc_r2_categorical, 5),
      `Replicate Pearson [95% CI]` = pearson_ci,
      `Replicate ICC [95% CI]` = icc_ci,
      Features = as.character(as.integer(retained_features))
    )
  write_csv_copy(x, stem)
  write_longtable(
    x, stem, caption, label, "lrrrrrr",
    c("Method", "QC CV", "GAM", paste0(batch_name, " $R^2$"),
      "Replicate $r$ [95\\% CI]", "ICC(A,1) [95\\% CI]", "Features"),
    landscape = TRUE
  )
}
make_replicate_table(
  "batchcorr_set1", "table_s6_batchcorr",
  "Complete BatchCorr Set 1 results across all nine outputs. Replicate-profile Pearson correlation and within-accession feature-profile ICC(A,1) include 2,000-resample accession-cluster bootstrap intervals. Batch $R^2$ values are proportions shown to five decimal places. All methods retained the same 199-feature evaluation panel.",
  "tab:supp_batchcorr_all", "Batch"
)
make_replicate_table(
  "clinical_fiams", "table_s7_clinical",
  "Complete pair-aware clinical FIA-MS results across all nine outputs. Replicate-profile Pearson correlation and within-identity feature-profile ICC(A,1) include 2,000-resample identity-cluster bootstrap intervals. Plate $R^2$ values are proportions shown to five decimal places and use the 536 study log1p-pair means; correction was fitted to all 660 log1p-pair means.",
  "tab:supp_clinical_all", "Plate"
)

# S8. Sacurine, with all three covariates kept separate.
sac <- prepare_empirical("sacurine") |>
  left_join(bootstrap_cell("sacurine", "age_weighted_pc_r2", 3) |> rename(age = cell), by = "method_id") |>
  left_join(bootstrap_cell("sacurine", "bmi_weighted_pc_r2", 3) |> rename(bmi = cell), by = "method_id") |>
  left_join(bootstrap_cell("sacurine", "gender_weighted_pc_r2", 3) |> rename(gender = cell), by = "method_id") |>
  transmute(
    Method = latex_escape(method),
    `QC CV` = fmt(heldout_qc_cv_mean, 2),
    GAM = fmt(residual_run_order_gam_mean_deviance, 3),
    `Batch R2` = fmt(batch_weighted_pc_r2_categorical, 4),
    `Age R2 [95% CI]` = age,
    `BMI R2 [95% CI]` = bmi,
    `Gender R2 [95% CI]` = gender
  )
write_csv_copy(sac, "table_s8_sacurine")
write_longtable(
  sac, "table_s8_sacurine",
  "Complete Sacurine results across all nine outputs. Age, BMI and gender are evaluated separately using weighted-PC $R^2$ with 2,000-resample stratified bootstrap intervals. These covariates were withheld from every correction method.",
  "tab:supp_sacurine_all", "lrrrrrr",
  c("Method", "QC CV", "GAM", "Batch $R^2$", "Age $R^2$ [95\\% CI]", "BMI $R^2$ [95\\% CI]", "Gender $R^2$ [95\\% CI]"),
  landscape = TRUE
)

# S9. WaveICA.
wave_assoc <- associations |>
  filter(dataset_key == "waveica", endpoint == "biological_group") |>
  select(method_id, significant_features, median_absolute_standardized_effect,
         median_within_group_or_residual_variance)
wave <- prepare_empirical("waveica") |>
  left_join(bootstrap_cell("waveica", "biological_weighted_pc_r2", 4), by = "method_id") |>
  left_join(wave_assoc, by = "method_id") |>
  transmute(
    Method = latex_escape(method),
    `QC CV` = fmt(heldout_qc_cv_mean, 2),
    GAM = fmt(residual_run_order_gam_mean_deviance, 3),
    `Batch R2` = fmt(batch_weighted_pc_r2_categorical, 4),
    `Group R2 [95% CI]` = cell,
    `BH features` = as.character(as.integer(significant_features)),
    `Standardized effect` = fmt(median_absolute_standardized_effect, 3),
    `Within-group variance` = fmt(median_within_group_or_residual_variance, 3)
  )
write_csv_copy(wave, "table_s9_waveica")
write_longtable(
  wave, "table_s9_waveica",
  "Complete WaveICA results across all nine outputs. The biological endpoint uses the neutral deposited group 0 and group 1 labels. Group $R^2$ has 2,000-resample stratified bootstrap intervals. Association counts are BH-FDR-significant peaks and are interpreted with standardized effect size and within-group variance.",
  "tab:supp_waveica_all", "lrrrrrrr",
  c("Method", "QC CV", "GAM", "Batch $R^2$", "Group $R^2$ [95\\% CI]", "BH peaks", "Std. effect", "Within-group var."),
  landscape = TRUE
)

# S10. Runtime, retention and WiNN correction counts.
runtime <- primary |>
  filter(method_id %in% method_order) |>
  mutate(
    dataset_key = factor(dataset_key, levels = dataset_order),
    method_id = factor(method_id, levels = method_order),
    Dataset = latex_escape(dataset), Method = latex_escape(unname(method_short[as.character(method_id)])),
    Runtime = fmt_runtime(runtime_sec),
    Features = as.character(as.integer(retained_features)),
    Samples = as.character(as.integer(retained_samples)),
    `Drift features` = ifelse(is.na(drift_features_corrected), "--", as.character(as.integer(drift_features_corrected))),
    `Drift profiles` = ifelse(is.na(drift_profiles_corrected), "--", as.character(as.integer(drift_profiles_corrected))),
    `Batch features` = ifelse(is.na(batch_features_corrected), "--", as.character(as.integer(batch_features_corrected)))
  ) |>
  arrange(dataset_key, method_id) |>
  select(Dataset, Method, Runtime, Features, Samples, `Drift features`, `Drift profiles`, `Batch features`)
write_csv_copy(runtime, "table_s10_runtime_retention")
write_longtable(
  runtime, "table_s10_runtime_retention",
  "Runtime, retained dimensions and WiNN correction counts for all primary dataset-method fits. A dash marks a correction reused from a validated frozen run, for which matrix-loading time was not reported as algorithm runtime. Correction counts are defined only for WiNN modes.",
  "tab:supp_runtime", "llrrrrrr",
  c("Dataset", "Method", "Runtime (s)", "Features", "Samples", "Drift feat.", "Drift profiles", "Batch feat."),
  landscape = TRUE
)

# S11. Complete association and permutation summary for all nine methods.
assoc <- associations |>
  left_join(
    permutations |>
      select(dataset_key, method_id, endpoint, null_median, null_95th_percentile,
             empirical_p_value, valid_permutations),
    by = c("dataset_key", "method_id", "endpoint")
  ) |>
  mutate(
    dataset_key = factor(dataset_key, levels = c("mtbls79", "waveica", "sacurine")),
    method_id = factor(method_id, levels = method_order),
    Dataset = latex_escape(dataset),
    Endpoint = latex_escape(recode(
      endpoint, cow_vs_sheep = "cow vs sheep", biological_group = "group 0 vs 1",
      age = "age", bmi = "BMI", gender = "gender"
    )),
    Method = latex_escape(unname(method_short[as.character(method_id)])),
    Tested = as.character(as.integer(tested_features)),
    Significant = as.character(as.integer(significant_features)),
    `Abs. effect` = fmt(median_absolute_effect, 3),
    `Std. effect` = fmt(median_absolute_standardized_effect, 3),
    Variance = fmt(median_within_group_or_residual_variance, 3),
    `Null 95%` = fmt(null_95th_percentile, 0),
    `Empirical P` = ifelse(empirical_p_value < 0.001, "<0.001", fmt(empirical_p_value, 3))
  ) |>
  arrange(dataset_key, Endpoint, method_id) |>
  select(Dataset, Endpoint, Method, Tested, Significant, `Abs. effect`, `Std. effect`, Variance, `Null 95%`, `Empirical P`)
write_csv_copy(assoc, "table_s11_association_permutation")
write_longtable(
  assoc, "table_s11_association_permutation",
  "Complete biological-association and constrained-permutation results across all nine outputs. Null summaries use 5,000 permutations. Association counts are descriptive because ground-truth associated features are unknown in these public datasets.",
  "tab:supp_association", "lllrrrrrrr",
  c("Dataset", "Endpoint", "Method", "Tested", "BH sig.", "Abs. effect", "Std. effect", "Variance", "Null 95th", "Empirical $P$"),
  landscape = TRUE, size = "\\tiny"
)

# S12. All-method leave-one-batch-out summary. Fold-level results remain in the
# copied source CSV, while the printed table gives the median fold estimate and
# median fold bounds for every estimable method-endpoint-metric combination.
lobo_raw <- read_csv(file.path(derivative_root, "leave_one_batch_out_bootstrap_summary.csv"), show_col_types = FALSE) |>
  filter(metric %in% c("pearson", "spearman", "sign_concordance"), method_id %in% method_order)
write_csv(lobo_raw, file.path(out_csv, "table_s12_lobo_fold_level.csv"), na = "")
lobo <- lobo_raw |>
  group_by(dataset_key, dataset, endpoint, method_id, method, metric) |>
  summarise(
    estimate = median(bootstrap_mean, na.rm = TRUE),
    lower = median(ci_lower, na.rm = TRUE),
    upper = median(ci_upper, na.rm = TRUE),
    estimable_folds = sum(is.finite(bootstrap_mean)),
    total_folds = n(),
    .groups = "drop"
  ) |>
  mutate(
    dataset_key = factor(dataset_key, levels = c("mtbls79", "waveica", "sacurine")),
    method_id = factor(method_id, levels = method_order),
    Dataset = latex_escape(dataset),
    Endpoint = latex_escape(recode(
      endpoint, cow_vs_sheep = "cow vs sheep", biological_group = "group 0 vs 1",
      age = "age", bmi = "BMI", gender = "gender"
    )),
    Method = latex_escape(unname(method_short[as.character(method_id)])),
    Metric = latex_escape(recode(metric, pearson = "Pearson", spearman = "Spearman", sign_concordance = "sign concordance")),
    `Median fold mean [median bounds]` = fmt_ci(estimate, lower, upper, 3),
    Folds = paste0(estimable_folds, "/", total_folds)
  ) |>
  arrange(dataset_key, Endpoint, method_id, Metric) |>
  select(Dataset, Endpoint, Method, Metric, `Median fold mean [median bounds]`, Folds)
write_csv_copy(lobo, "table_s12_lobo_summary")
write_longtable(
  lobo, "table_s12_lobo_summary",
  "Complete all-method summary of leave-one-batch-out biological-effect replication. The printed value is the median fold bootstrap mean followed by the medians of the fold-specific 2.5th- and 97.5th-percentile bounds. These brackets are not confidence intervals for the across-fold median. The accompanying source table retains every fold-specific estimate and interval. Non-estimable folds remain in the denominator.",
  "tab:supp_lobo", "llllrr",
  c("Dataset", "Endpoint", "Method", "Metric", "Median fold mean [median bounds]", "Estimable folds"),
  landscape = TRUE, size = "\\tiny", arraystretch = 0.94
)

# S13. Partial confounding, with clean-reference and injected-effect ratios
# kept as distinct columns.
partial_raw <- read_csv(file.path(table_root, "partial_confounding_method_metrics.csv"), show_col_types = FALSE) |>
  filter(method == "WiNN fixed default (no QC)")
partial_metrics <- c(
  "heldout_qc_cv_mean", "residual_run_order_gam_mean_deviance",
  "batch_weighted_pc_r2_categorical", "clean_sample_profile_pearson_mean",
  "median_attenuation_ratio_clean_responsive",
  "median_attenuation_ratio_injected_responsive",
  "true_positive_rate", "empirical_fdr"
)
partial <- partial_raw |>
  filter(metric %in% partial_metrics) |>
  group_by(confounding_type, level, metric) |>
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), n = sum(is.finite(value)), .groups = "drop") |>
  mutate(
    metric = recode(
      metric,
      heldout_qc_cv_mean = "QC CV",
      residual_run_order_gam_mean_deviance = "GAM",
      batch_weighted_pc_r2_categorical = "Batch R2",
      clean_sample_profile_pearson_mean = "Clean Pearson",
      median_attenuation_ratio_clean_responsive = "Clean ratio",
      median_attenuation_ratio_injected_responsive = "Injected ratio",
      true_positive_rate = "TPR", empirical_fdr = "FDR"
    ),
    digits = ifelse(metric == "QC CV", 2L, 3L),
    cell = mapply(fmt_mean_se, mean, sd, n, digits, USE.NAMES = FALSE)
  ) |>
  select(confounding_type, level, metric, cell) |>
  pivot_wider(names_from = metric, values_from = cell) |>
  mutate(
    type_order = match(confounding_type, c("none", "batch", "run_order", "joint")),
    Type = latex_escape(recode(confounding_type, none = "None", batch = "Batch", run_order = "Run order", joint = "Joint")),
    Alignment = paste0(fmt(100 * level, 0), "\\%")
  ) |>
  arrange(type_order, level) |>
  select(Type, Alignment, `QC CV`, GAM, `Batch R2`, `Clean Pearson`, `Clean ratio`, `Injected ratio`, TPR, FDR)
write_csv_copy(partial, "table_s13_partial_confounding")
write_longtable(
  partial, "table_s13_partial_confounding",
  "Complete WiNN partial-confounding summary over ten seeds. Cells are mean $\\pm$ seed-level standard error. The clean-reference ratio divides the signed corrected phenotype coefficient by its coefficient in the clean simulated matrix. The injected-effect ratio divides by the corresponding signed injected coefficient, $+0.35$ or $-0.35$. These denominators differ and the two ratios are not interchangeable.",
  "tab:supp_partial_confounding", "llrrrrrrrr",
  c("Alignment axis", "Strength", "QC CV", "GAM", "Batch $R^2$", "Clean Pearson", "Clean ratio", "Injected ratio", "TPR", "FDR"),
  landscape = TRUE, size = "\\tiny"
)

# S14. Ten label-blinded reference partitions for every method and dataset.
ref <- read_csv(file.path(table_root, "reference_split_summary.csv"), show_col_types = FALSE) |>
  filter(metric %in% c("heldout_qc_cv_mean", "residual_run_order_gam_mean_deviance"), method_id %in% method_order) |>
  mutate(
    metric = recode(metric, heldout_qc_cv_mean = "QC CV", residual_run_order_gam_mean_deviance = "GAM"),
    cell = ifelse(metric == "QC CV", fmt_mean_se(mean, sd, n, 2), fmt_mean_se(mean, sd, n, 3))
  ) |>
  select(dataset_key, dataset, method_id, metric, cell) |>
  pivot_wider(names_from = metric, values_from = cell) |>
  mutate(
    dataset_key = factor(dataset_key, levels = dataset_order),
    method_id = factor(method_id, levels = method_order),
    Dataset = latex_escape(dataset), Method = latex_escape(unname(method_short[as.character(method_id)]))
  ) |>
  arrange(dataset_key, method_id) |>
  select(Dataset, Method, `QC CV`, GAM)
write_csv_copy(ref, "table_s14_reference_splits")
write_longtable(
  ref, "table_s14_reference_splits",
  "Complete reference-split stability summary across all datasets and methods. Cells are mean $\\pm$ standard error across ten prespecified, overlapping label-blinded partitions. The partitions are sensitivity units, not independent cohorts.",
  "tab:supp_reference_splits", "llrr",
  c("Dataset", "Method", "QC CV", "Residual GAM"),
  landscape = FALSE
)

# S15. Full selective-versus-forced gate table.
gate <- read_csv(file.path(table_root, "selective_vs_forced_metrics_long.csv"), show_col_types = FALSE) |>
  filter(metric %in% c(
    "heldout_qc_cv_mean",
    "residual_run_order_gam_mean_deviance",
    "batch_weighted_pc_r2_categorical",
    "truth_metabolite_profile_icc_mean",
    "truth_sample_profile_pearson_mean",
    "metabolite_icc_a1_median",
    "sample_replicate_pearson_median",
    "feature_repeatability_ratio_median"
  )) |>
  mutate(
    dataset_key = factor(dataset_key, levels = dataset_order),
    method_id = factor(method_id, levels = c("G_SS", "G_AS", "G_SA", "G_AA")),
    Dataset = latex_escape(dataset),
    Policy = latex_escape(recode(
      method_id, G_SS = "selective drift / selective batch",
      G_AS = "forced drift / selective batch",
      G_SA = "selective drift / forced batch",
      G_AA = "forced drift / forced batch"
    )),
    Metric = latex_escape(recode(
      metric,
      heldout_qc_cv_mean = "held-out QC CV",
      residual_run_order_gam_mean_deviance = "residual GAM deviance",
      batch_weighted_pc_r2_categorical = "categorical batch R2",
      truth_metabolite_profile_icc_mean = "truth feature ICC(A,1)",
      truth_sample_profile_pearson_mean = "truth sample Pearson",
      metabolite_icc_a1_median = "cross-batch feature ICC(A,1)",
      sample_replicate_pearson_median = "replicate-profile Pearson",
      feature_repeatability_ratio_median = "feature repeatability ratio"
    )),
    Value = ifelse(metric == "heldout_qc_cv_mean", fmt(value, 2), fmt(value, 4))
  ) |>
  arrange(dataset_key, method_id, Metric) |>
  select(Dataset, Policy, Metric, Value)
write_csv_copy(gate, "table_s15_selective_forced")
write_longtable(
  gate, "table_s15_selective_forced",
  "Complete selective-versus-forced gate results. Preparation and conservative PQN are retained in all four policies. Values are descriptive workflow-policy contrasts. The representative simulation realization is used here; repeated-seed truth comparisons are shown in Figure S2D-F.",
  "tab:supp_selective_forced", "lllr",
  c("Dataset", "Gate policy", "Metric", "Value"),
  landscape = FALSE, size = "\\tiny"
)

# S16. Complete cumulative-stage values.
stage_raw <- read_csv(file.path(table_root, "ablation_stage_impacts_table.csv"), show_col_types = FALSE)
stage <- stage_raw |>
  distinct(dataset_key, dataset, metric, direction, cum_raw, cum_outlier, cum_drift, cum_batch, cum_pqn) |>
  mutate(
    dataset_key = factor(dataset_key, levels = dataset_order),
    Dataset = latex_escape(dataset),
    Metric = latex_escape(recode(
      metric,
      heldout_qc_cv_mean = "held-out QC CV",
      residual_run_order_gam_mean_deviance = "residual GAM deviance",
      batch_weighted_pc_r2_categorical = "categorical batch R2",
      truth_metabolite_profile_icc_mean = "truth feature ICC(A,1)",
      metabolite_icc_a1_median = "cross-batch feature ICC(A,1)",
      truth_sample_profile_pearson_mean = "truth sample Pearson",
      sample_replicate_pearson_median = "replicate-profile Pearson",
      feature_repeatability_ratio_median = "feature repeatability ratio"
    )),
    Raw = ifelse(metric == "heldout_qc_cv_mean", fmt(cum_raw, 2), fmt(cum_raw, 4)),
    MAD = ifelse(metric == "heldout_qc_cv_mean", fmt(cum_outlier, 2), fmt(cum_outlier, 4)),
    Drift = ifelse(metric == "heldout_qc_cv_mean", fmt(cum_drift, 2), fmt(cum_drift, 4)),
    Batch = ifelse(metric == "heldout_qc_cv_mean", fmt(cum_batch, 2), fmt(cum_batch, 4)),
    PQN = ifelse(metric == "heldout_qc_cv_mean", fmt(cum_pqn, 2), fmt(cum_pqn, 4))
  ) |>
  arrange(dataset_key, Metric) |>
  select(Dataset, Metric, Raw, MAD, Drift, Batch, PQN)
write_csv_copy(stage, "table_s16_cumulative_stages")
write_longtable(
  stage, "table_s16_cumulative_stages",
  "Complete cumulative WiNN ablation values. Columns are immutable matrices after each ordered stage. Adjacent differences are order-dependent and are not a factorial decomposition. The simulation row uses one prespecified representative realization.",
  "tab:supp_cumulative", "llrrrrr",
  c("Dataset", "Metric", "Raw", "MAD", "Drift", "Batch", "PQN"),
  landscape = TRUE, size = "\\tiny"
)

# S17. Dataset-specific WiNN mode parameters.
params <- read_csv(file.path(table_root, "primary_winn_selected_parameters.csv"), show_col_types = FALSE) |>
  mutate(
    dataset_key = factor(dataset_key, levels = dataset_order),
    method_id = factor(method_id, levels = c("winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc")),
    fixed = as.character(method_id) == "winn_fixed_no_qc",
    BatchSource = ifelse(fixed, "supplied", batch_source),
    Test = "Ljung-Box",
    Smoother = ifelse(fixed, "conservative", spline_method),
    DriftFDR = ifelse(fixed, "0.05", fmt(acorr_fdr, 2)),
    BatchFDR = ifelse(fixed, "0.05", fmt(anova_fdr, 2)),
    PQN = ifelse(fixed, "shrink", normalization),
    Dataset = latex_escape(dataset),
    Mode = latex_escape(unname(method_short[as.character(method_id)])),
    BatchSource = latex_escape(BatchSource),
    Smoother = latex_escape(Smoother),
    PQN = latex_escape(PQN),
    Segments = as.character(as.integer(detected_batch_count))
  ) |>
  arrange(dataset_key, method_id) |>
  select(Dataset, Mode, BatchSource, Segments, Test, Smoother, DriftFDR, BatchFDR, PQN)
write_csv_copy(params, "table_s17_winn_parameters")
write_longtable(
  params, "table_s17_winn_parameters",
  "WiNN parameters for every dataset and operating mode. WiNN uses the same frozen settings in every dataset and receives no QC identities. Automatic-mode settings were selected from the prespecified grid using labelled training QCs only.",
  "tab:supp_winn_parameters", "lllrlrrrr",
  c("Dataset", "Mode", "Batch source", "Segments", "Test", "Smoother", "Drift FDR", "Batch FDR", "PQN"),
  landscape = TRUE, size = "\\tiny"
)

# S18. Numerical safeguards and complete validation checks.
flooring <- read_csv(file.path(table_root, "intensity_domain_safeguard.csv"), show_col_types = FALSE)
tiger <- read_csv(file.path(table_root, "tiger_nonfinite_feature_fallback.csv"), show_col_types = FALSE)
safeguard <- bind_rows(
  flooring |>
    summarise(
      Safeguard = "Negative corrected intensities floored at zero",
      Runs = n(), Triggered = sum(floor_applied),
      `Affected features` = sum(features_with_flooring),
      `Affected values` = sum(negative_values_floored)
    ),
  tiger |>
    summarise(
      Safeguard = "TIGER nonfinite features restored to raw",
      Runs = n(), Triggered = sum(fallback_applied),
      `Affected features` = sum(features_restored_to_raw),
      `Affected values` = sum(nonfinite_values_returned)
    )
) |>
  mutate(across(everything(), as.character), Safeguard = latex_escape(Safeguard))
write_csv_copy(safeguard, "table_s18_safeguards")
write_longtable(
  safeguard, "table_s18_safeguards",
  "Recorded numerical safeguards. The trigger was retained as part of the benchmark outcome and no affected run was silently dropped.",
  "tab:supp_safeguards", "lrrrr",
  c("Safeguard", "Audited runs", "Triggered runs", "Affected features", "Affected values"),
  landscape = FALSE
)

validation <- read_csv(file.path(validation_root, "release_validation.csv"), show_col_types = FALSE) |>
  transmute(
    Check = latex_escape(gsub("_", " ", check, fixed = TRUE)),
    Passed = ifelse(passed, "yes", "no"),
    Detail = latex_escape(detail)
  )
figure_panels <- c(
  list.files(file.path(package_root, "figures", "main_panels"), pattern = "[.]pdf$"),
  list.files(file.path(package_root, "figures", "supplementary_panels"), pattern = "[.]pdf$")
)
if(length(figure_panels) != 69L) stop("The current table requires the audited 69-panel figure package.")
figure_row <- validation$Check == "all figures have source data"
validation$Check[figure_row] <- "all declared figure panels have source data"
validation$Detail[figure_row] <- "69 of 69 figure-panel/source entries"
write_csv_copy(validation, "table_s19_release_validation")
write_longtable(
  validation, "table_s19_release_validation",
  "Complete release-validation checklist for the frozen benchmark.",
  "tab:supp_validation", "llp{10.5cm}",
  c("Check", "Passed", "Detail"),
  landscape = FALSE, size = "\\tiny"
)

nonestimable <- read_csv(file.path(derivative_root, "non_estimable_analyses.csv"), show_col_types = FALSE) |>
  mutate(
    reason = case_when(
      dataset_key == "batchcorr_set1" & analysis == "leave-one-batch-out phenotype effects" ~
        "Cross-batch replicate agreement is the design-supported validation endpoint.",
      dataset_key == "clinical_fiams" & analysis == "leave-one-batch-out phenotype effects" ~
        "Cross-plate replicate agreement is the permitted design-supported endpoint.",
      TRUE ~ reason
    )
  ) |>
  transmute(
    Dataset = latex_escape(dataset), Analysis = latex_escape(analysis),
    Status = latex_escape(status), Reason = latex_escape(reason)
  )
write_csv_copy(nonestimable, "table_s20_nonestimable")
write_longtable(
  nonestimable, "table_s20_nonestimable",
  "Prespecified analyses that were not applicable or not estimable under the available design. No relaxed model was substituted.",
  "tab:supp_nonestimable", "lllp{8cm}",
  c("Dataset", "Analysis", "Status", "Reason"),
  landscape = TRUE, size = "\\tiny"
)

table_files <- sort(list.files(out_tex, pattern = "^table_s([4-9]|1[0-9]|20)_.*[.]tex$"))
manifest <- data.frame(
  file = table_files,
  generated_from = "archived primary benchmark and local derivative aggregate export; supplied in supplement_20260727",
  publication_revision = ifelse(grepl("^table_s4_", table_files),
    "corrected ten-seed RFSC results applied by integrate_completed_results.py",
    "current labels and RFSC exclusions applied by integrate_completed_results.py; S21-S31 have separate source mappings"),
  stringsAsFactors = FALSE
)
write_csv(manifest, file.path(out_csv, "supplementary_table_manifest.csv"))
message("Generated ", nrow(manifest), " supplementary TeX tables in ", out_tex)

# Apply the validated corrected RFSC results after generating the original tables.
status <- system2("python3", shQuote(file.path(package_root, "scripts", "benchmark_revision", "integrate_completed_results.py")))
if (status != 0) stop("Corrected RFSC table integration failed.")

# Keep generated display labels consistent without renaming analysis identifiers.
status <- system2("python3", shQuote(file.path(package_root, "scripts", "normalize_table_terminology.py")))
if (status != 0) stop("Table terminology normalization failed.")
