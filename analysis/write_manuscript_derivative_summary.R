#!/usr/bin/env Rscript

# Generate manuscript-facing tables from frozen primary and derivative results.
# This script reads saved metrics only. It never sources correction, tuning, or
# selection code and cannot create a corrected matrix.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
final_root <- file.path(repo_root, "results", "final")
derivative_root <- file.path(final_root, "derivative")
output_dir <- file.path(derivative_root, "manuscript")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("Missing package: digest", call. = FALSE)

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required saved result: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_lines <- function(lines, filename) {
  writeLines(lines, file.path(output_dir, filename), useBytes = TRUE)
}

fmt <- function(x, digits = 3L) {
  ifelse(is.na(x), "--", formatC(x, digits = digits, format = "f"))
}

method_ids <- c("raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf", "winn_fixed_no_qc")
method_labels <- c(
  raw = "Raw", combat = "ComBat", qc_rlsc = "QC-RLSC", qc_rfsc = "QC-RFSC",
  tiger = "TIGER", serrf = "SERRF", winn_fixed_no_qc = "WiNN fixed"
)

primary <- read_required(file.path(final_root, "tables", "primary_method_metrics_wide.csv"))
primary_long <- read_required(file.path(final_root, "tables", "primary_method_metrics_long.csv"))
simulation <- read_required(file.path(final_root, "tables", "simulation_seed_summary.csv"))
bootstrap <- read_required(file.path(derivative_root, "tables", "bootstrap_estimates.csv"))
permutation <- read_required(file.path(derivative_root, "tables", "permutation_summary.csv"))
associations <- read_required(file.path(derivative_root, "tables", "biological_association_summary.csv"))
lobo <- read_required(file.path(derivative_root, "tables", "leave_one_batch_out_bootstrap_summary.csv"))
ablation_impacts <- read_required(file.path(final_root, "tables", "ablation_stage_impacts.csv"))
reference_split_metrics <- read_required(file.path(final_root, "tables", "reference_split_method_metrics.csv"))
partial_confounding <- read_required(file.path(final_root, "tables", "partial_confounding_summary.csv"))

sim_metrics <- c(
  "heldout_qc_cv_mean", "residual_run_order_gam_mean_deviance",
  "truth_metabolite_profile_icc_mean", "truth_sample_profile_pearson_mean"
)
sim <- simulation[simulation$method_id %in% method_ids & simulation$metric %in% sim_metrics, ]
if (any(sim$n != 10L) || nrow(sim) != length(method_ids) * length(sim_metrics)) {
  stop("The frozen simulation summary is not a complete 10-seed panel.", call. = FALSE)
}
sim_wide <- reshape(
  sim[c("method_id", "metric", "mean", "sd")], idvar = "method_id",
  timevar = "metric", direction = "wide"
)
sim_wide <- sim_wide[match(method_ids, sim_wide$method_id), ]
write.csv(
  sim[order(match(sim$method_id, method_ids), match(sim$metric, sim_metrics)), ],
  file.path(output_dir, "main_simulation_rows.csv"), row.names = FALSE, na = ""
)
simulation_rows <- vapply(seq_len(nrow(sim_wide)), function(index) {
  row <- sim_wide[index, ]
  label <- method_labels[[row$method_id]]
  if (row$method_id == "winn_fixed_no_qc") label <- paste0("\\textbf{", label, "}")
  values <- vapply(sim_metrics, function(metric) {
    paste0(
      fmt(row[[paste0("mean.", metric)]], if (metric == "heldout_qc_cv_mean") 2L else 3L),
      "$\\pm$",
      fmt(row[[paste0("sd.", metric)]], if (metric == "heldout_qc_cv_mean") 2L else 3L)
    )
  }, character(1))
  paste(c(label, values), collapse = " & ")
}, character(1))
write_lines(c(paste0(simulation_rows, " \\\\"), "\\bottomrule"), "main_simulation_rows.tex")

make_primary_rows <- function(dataset_key, biological_source, biological_metrics, filename) {
  technical <- primary[primary$dataset_key == dataset_key & primary$method_id %in% method_ids, ]
  technical <- technical[match(method_ids, technical$method_id), ]
  if (nrow(technical) != length(method_ids)) stop("Incomplete primary panel for ", dataset_key, call. = FALSE)
  values <- biological_source[
    biological_source$dataset_key == dataset_key & biological_source$method_id %in% method_ids &
      biological_source$metric %in% biological_metrics, , drop = FALSE
  ]
  rows <- vapply(seq_len(nrow(technical)), function(index) {
    row <- technical[index, ]
    label <- method_labels[[row$method_id]]
    if (row$method_id == "winn_fixed_no_qc") label <- paste0("\\textbf{", label, "}")
    biological <- vapply(biological_metrics, function(metric) {
      found <- values[values$method_id == row$method_id & values$metric == metric, ]
      if (nrow(found) != 1L) stop("Missing biological metric ", dataset_key, "/", row$method_id, "/", metric, call. = FALSE)
      source_column <- if ("point_estimate" %in% names(found)) "point_estimate" else "value"
      fmt(found[[source_column]], 3L)
    }, character(1))
    paste(c(
      label, fmt(row$heldout_qc_cv_mean, 2L), fmt(row$residual_run_order_gam_mean_deviance, 3L), biological
    ), collapse = " & ")
  }, character(1))
  write_lines(c(paste0(rows, " \\\\"), "\\bottomrule"), filename)
}

make_primary_rows(
  "mtbls79", primary_long, c("metabolite_icc_a1_median", "sample_replicate_pearson_median"),
  "main_mtbls79_rows.tex"
)
make_primary_rows(
  "clinical_fiams", bootstrap,
  c("genuine_replicate_icc_median", "sample_replicate_pearson_median"),
  "main_clinical_rows.tex"
)
write.csv(
  primary_long[
    primary_long$dataset_key == "mtbls79" & primary_long$method_id %in% method_ids &
      primary_long$metric %in% c(
        "heldout_qc_cv_mean", "residual_run_order_gam_mean_deviance",
        "metabolite_icc_a1_median", "sample_replicate_pearson_median"
      ),
  ],
  file.path(output_dir, "main_mtbls79_rows.csv"), row.names = FALSE, na = ""
)
write.csv(
  rbind(
    primary_long[
      primary_long$dataset_key == "clinical_fiams" & primary_long$method_id %in% method_ids &
        primary_long$metric %in% c("heldout_qc_cv_mean", "residual_run_order_gam_mean_deviance"),
    ],
    transform(
      bootstrap[
        bootstrap$dataset_key == "clinical_fiams" & bootstrap$method_id %in% method_ids &
          bootstrap$metric %in% c("genuine_replicate_icc_median", "sample_replicate_pearson_median"),
        c("metric", "point_estimate", "direction", "method_id", "method", "dataset_key", "dataset")
      ],
      value = point_estimate, denominator = NA_real_, units = NA_character_, notes = NA_character_
    )[names(primary_long)]
  ),
  file.path(output_dir, "main_clinical_rows.csv"), row.names = FALSE, na = ""
)

added_specs <- list(
  batchcorr_set1 = c("genuine_replicate_icc_median"),
  sacurine = c("age_weighted_pc_r2", "bmi_weighted_pc_r2", "gender_weighted_pc_r2"),
  waveica = c("biological_weighted_pc_r2")
)
dataset_labels <- c(
  simulation = "Simulation", mtbls79 = "MTBLS79", clinical_fiams = "Clinical FIA-MS",
  batchcorr_set1 = "BatchCorr Set 1", sacurine = "Sacurine", waveica = "WaveICA"
)
added_rows <- vapply(names(added_specs), function(dataset_key) {
  tech <- primary[
    primary$dataset_key == dataset_key & primary$method_id %in% c("raw", "winn_fixed_no_qc"),
  ]
  tech <- tech[match(c("raw", "winn_fixed_no_qc"), tech$method_id), ]
  bio <- bootstrap[
    bootstrap$dataset_key == dataset_key & bootstrap$method_id %in% c("raw", "winn_fixed_no_qc") &
      bootstrap$metric %in% added_specs[[dataset_key]],
  ]
  bio_values <- vapply(added_specs[[dataset_key]], function(metric) {
    value <- bio[bio$metric == metric, ]
    value <- value[match(c("raw", "winn_fixed_no_qc"), value$method_id), ]
    paste0(fmt(value$point_estimate[1], 3L), "$\\rightarrow$", fmt(value$point_estimate[2], 3L))
  }, character(1))
  preservation <- paste(bio_values, collapse = "; ")
  paste(c(
    dataset_labels[[dataset_key]],
    paste0(fmt(tech$heldout_qc_cv_mean[1], 2L), "$\\rightarrow$", fmt(tech$heldout_qc_cv_mean[2], 2L)),
    paste0(fmt(tech$residual_run_order_gam_mean_deviance[1], 3L), "$\\rightarrow$", fmt(tech$residual_run_order_gam_mean_deviance[2], 3L)),
    paste0(fmt(tech$batch_weighted_pc_r2_categorical[1], 4L), "$\\rightarrow$", fmt(tech$batch_weighted_pc_r2_categorical[2], 4L)),
    preservation,
    paste0(tech$retained_features[2], "/", tech$retained_features[1])
  ), collapse = " & ")
}, character(1))
write_lines(c(paste0(added_rows, " \\\\"), "\\bottomrule"), "main_added_public_rows.tex")

long_primary_metrics <- function(dataset_keys, selected_methods, metrics) {
  rows <- primary[primary$dataset_key %in% dataset_keys & primary$method_id %in% selected_methods, ]
  do.call(rbind, lapply(seq_len(nrow(rows)), function(index) {
    do.call(rbind, lapply(metrics, function(metric) data.frame(
      dataset_key = rows$dataset_key[index], dataset = rows$dataset[index],
      method_id = rows$method_id[index], method = rows$method[index],
      domain = "technical", metric = metric, value = rows[[metric]][index],
      stringsAsFactors = FALSE
    )))
  }))
}
added_source <- rbind(
  long_primary_metrics(
    names(added_specs), c("raw", "winn_fixed_no_qc"),
    c("heldout_qc_cv_mean", "residual_run_order_gam_mean_deviance", "batch_weighted_pc_r2_categorical", "retained_features")
  ),
  transform(
    bootstrap[
      bootstrap$dataset_key %in% names(added_specs) & bootstrap$method_id %in% c("raw", "winn_fixed_no_qc") &
        bootstrap$metric %in% unique(unlist(added_specs)),
      c("dataset_key", "dataset", "method_id", "method", "metric", "point_estimate")
    ],
    domain = "biological", value = point_estimate
  )[c("dataset_key", "dataset", "method_id", "method", "domain", "metric", "value")]
)
write.csv(added_source, file.path(output_dir, "main_added_public_rows.csv"), row.names = FALSE, na = "")

key_specs <- data.frame(
  dataset_key = c(
    "simulation", "simulation", "simulation", "mtbls79", "batchcorr_set1", "batchcorr_set1",
    "sacurine", "sacurine", "sacurine", "waveica", "clinical_fiams", "clinical_fiams"
  ),
  metric = c(
    "truth_log1p_rmse", "truth_metabolite_profile_icc_mean", "truth_sample_profile_pearson_mean",
    "biological_weighted_pc_r2", "genuine_replicate_icc_median", "cross_batch_top1_accuracy",
    "age_weighted_pc_r2", "bmi_weighted_pc_r2", "gender_weighted_pc_r2",
    "biological_weighted_pc_r2", "genuine_replicate_icc_median", "cross_batch_top1_accuracy"
  ),
  endpoint = c(
    "Truth log1p RMSE", "Truth metabolite ICC(A,1)", "Truth sample Pearson",
    "Cow--sheep weighted-PC $R^2$", "Replicate ICC(A,1)", "Cross-batch top-1 matching",
    "Age weighted-PC $R^2$", "BMI weighted-PC $R^2$", "Gender weighted-PC $R^2$",
    "Deposited-group weighted-PC $R^2$", "Replicate ICC(A,1)", "Cross-plate top-1 matching"
  ), stringsAsFactors = FALSE
)
key_results <- do.call(rbind, lapply(seq_len(nrow(key_specs)), function(index) {
  specification <- key_specs[index, , drop = FALSE]
  rownames(specification) <- NULL
  selected <- bootstrap[
    bootstrap$dataset_key == specification$dataset_key & bootstrap$metric == specification$metric &
      bootstrap$method_id %in% c("raw", "winn_fixed_no_qc"),
  ]
  selected <- selected[match(c("raw", "winn_fixed_no_qc"), selected$method_id), ]
  if (nrow(selected) != 2L || any(is.na(selected$method_id))) {
    stop("Incomplete key result for ", specification$dataset_key, "/", specification$metric, call. = FALSE)
  }
  data.frame(
    specification, method_id = selected$method_id, method = selected$method,
    estimate = selected$point_estimate, ci_lower = selected$ci_lower, ci_upper = selected$ci_upper,
    stringsAsFactors = FALSE
  )
}))
write.csv(key_results, file.path(output_dir, "manuscript_key_biological_results.csv"), row.names = FALSE, na = "")

key_rows <- vapply(seq_len(nrow(key_specs)), function(index) {
  specification <- key_specs[index, ]
  selected <- key_results[
    key_results$dataset_key == specification$dataset_key & key_results$metric == specification$metric,
  ]
  raw <- selected[selected$method_id == "raw", ]
  fixed <- selected[selected$method_id == "winn_fixed_no_qc", ]
  display_ci <- function(row) paste0(fmt(row$estimate, 3L), " (", fmt(row$ci_lower, 3L), "--", fmt(row$ci_upper, 3L), ")")
  paste(c(
    dataset_labels[[specification$dataset_key]],
    specification$endpoint, display_ci(raw), display_ci(fixed)
  ), collapse = " & ")
}, character(1))
write_lines(c(paste0(key_rows, " \\\\"), "\\bottomrule"), "supp_derivative_biology_rows.tex")

lookup_boot <- function(dataset_key, method_id, metric) {
  found <- bootstrap[
    bootstrap$dataset_key == dataset_key & bootstrap$method_id == method_id & bootstrap$metric == metric,
  ]
  if (nrow(found) != 1L) stop("Missing bootstrap value: ", dataset_key, "/", method_id, "/", metric, call. = FALSE)
  found$point_estimate
}

lookup_primary <- function(dataset_key, method_id, column) {
  found <- primary[primary$dataset_key == dataset_key & primary$method_id == method_id, ]
  if (nrow(found) != 1L) stop("Missing primary value: ", dataset_key, "/", method_id, call. = FALSE)
  found[[column]]
}

noqc_specs <- list(
  mtbls79 = list(metric = "metabolite_icc_a1_median", label = "metabolite ICC(A,1)", source = "primary"),
  clinical_fiams = list(metric = "genuine_replicate_icc_median", label = "replicate ICC(A,1)", source = "bootstrap"),
  batchcorr_set1 = list(metric = "genuine_replicate_icc_median", label = "replicate ICC(A,1)", source = "bootstrap"),
  sacurine = list(metric = "age_weighted_pc_r2", label = "age weighted-PC $R^2$", source = "bootstrap"),
  waveica = list(metric = "biological_weighted_pc_r2", label = "group weighted-PC $R^2$", source = "bootstrap")
)
noqc_rows <- vapply(names(noqc_specs), function(dataset_key) {
  spec <- noqc_specs[[dataset_key]]
  technical <- primary[
    primary$dataset_key == dataset_key & primary$method_id %in% c("raw", "winn_fixed_no_qc"),
  ]
  technical <- technical[match(c("raw", "winn_fixed_no_qc"), technical$method_id), ]
  preservation <- if (spec$source == "primary") {
    c(
      lookup_primary(dataset_key, "raw", spec$metric),
      lookup_primary(dataset_key, "winn_fixed_no_qc", spec$metric)
    )
  } else {
    c(
      lookup_boot(dataset_key, "raw", spec$metric),
      lookup_boot(dataset_key, "winn_fixed_no_qc", spec$metric)
    )
  }
  reduction <- 100 * (technical$heldout_qc_cv_mean[1] - technical$heldout_qc_cv_mean[2]) /
    technical$heldout_qc_cv_mean[1]
  paste(c(
    dataset_labels[[dataset_key]],
    paste0(fmt(technical$heldout_qc_cv_mean[1], 1L), "$\\rightarrow$", fmt(technical$heldout_qc_cv_mean[2], 1L)),
    fmt(reduction, 1L),
    paste0(fmt(technical$batch_weighted_pc_r2_categorical[1], 3L), "$\\rightarrow$", fmt(technical$batch_weighted_pc_r2_categorical[2], 3L)),
    paste0(fmt(preservation[1], 3L), "$\\rightarrow$", fmt(preservation[2], 3L)),
    spec$label,
    paste0(technical$retained_features[2], "/", technical$retained_features[1])
  ), collapse = " & ")
}, character(1))

sim_raw_cv <- simulation[simulation$method_id == "raw" & simulation$metric == "heldout_qc_cv_mean", "mean"]
sim_fixed_cv <- simulation[simulation$method_id == "winn_fixed_no_qc" & simulation$metric == "heldout_qc_cv_mean", "mean"]
sim_raw_batch <- simulation[simulation$method_id == "raw" & simulation$metric == "batch_weighted_pc_r2_categorical", "mean"]
sim_fixed_batch <- simulation[simulation$method_id == "winn_fixed_no_qc" & simulation$metric == "batch_weighted_pc_r2_categorical", "mean"]
sim_raw_icc <- simulation[simulation$method_id == "raw" & simulation$metric == "truth_metabolite_profile_icc_mean", "mean"]
sim_fixed_icc <- simulation[simulation$method_id == "winn_fixed_no_qc" & simulation$metric == "truth_metabolite_profile_icc_mean", "mean"]
simulation_noqc <- paste(c(
  "Simulation", paste0(fmt(sim_raw_cv, 1L), "$\\rightarrow$", fmt(sim_fixed_cv, 1L)),
  fmt(100 * (sim_raw_cv - sim_fixed_cv) / sim_raw_cv, 1L),
  paste0(fmt(sim_raw_batch, 3L), "$\\rightarrow$", fmt(sim_fixed_batch, 3L)),
  paste0(fmt(sim_raw_icc, 3L), "$\\rightarrow$", fmt(sim_fixed_icc, 3L)),
  "truth ICC(A,1)", "1000/1000"
), collapse = " & ")
write_lines(c(paste0(c(noqc_rows, simulation_noqc), " \\\\"), "\\bottomrule"), "main_noqc_rows.tex")
noqc_source <- rbind(
  long_primary_metrics(
    names(noqc_specs), c("raw", "winn_fixed_no_qc"),
    c("heldout_qc_cv_mean", "batch_weighted_pc_r2_categorical", "retained_features")
  ),
  added_source[
    added_source$domain == "biological" & added_source$dataset_key %in% names(noqc_specs),
  ]
)
write.csv(noqc_source, file.path(output_dir, "main_noqc_rows.csv"), row.names = FALSE, na = "")

full_method_ids <- c(
  "raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf",
  "winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc"
)
full_labels <- c(
  raw = "Raw", combat = "ComBat", qc_rlsc = "QC-RLSC", qc_rfsc = "QC-RFSC",
  tiger = "TIGER", serrf = "SERRF", winn_auto_qc = "WiNN auto (QC)",
  winn_auto_batch_qc = "WiNN auto-batch (QC)", winn_fixed_no_qc = "WiNN fixed"
)
full_dataset_order <- c("simulation", "mtbls79", "batchcorr_set1", "sacurine", "waveica", "clinical_fiams")
preservation_display <- function(dataset_key, method_id) {
  if (dataset_key == "simulation") {
    return(simulation[
      simulation$method_id == method_id & simulation$metric == "truth_metabolite_profile_icc_mean", "mean"
    ])
  }
  if (dataset_key == "mtbls79") return(lookup_primary(dataset_key, method_id, "metabolite_icc_a1_median"))
  if (dataset_key %in% c("batchcorr_set1", "clinical_fiams")) {
    return(lookup_boot(dataset_key, method_id, "genuine_replicate_icc_median"))
  }
  if (dataset_key == "waveica") return(lookup_boot(dataset_key, method_id, "biological_weighted_pc_r2"))
  if (dataset_key == "sacurine") {
    return(paste(vapply(c("age_weighted_pc_r2", "bmi_weighted_pc_r2", "gender_weighted_pc_r2"), function(metric) {
      fmt(lookup_boot(dataset_key, method_id, metric), 3L)
    }, character(1)), collapse = "/"))
  }
  stop("Unknown dataset: ", dataset_key, call. = FALSE)
}

full_records <- do.call(rbind, lapply(full_dataset_order, function(dataset_key) {
  do.call(rbind, lapply(full_method_ids, function(method_id) {
    if (dataset_key == "simulation") {
      metric_value <- function(metric) simulation[
        simulation$method_id == method_id & simulation$metric == metric, "mean"
      ]
      cv <- metric_value("heldout_qc_cv_mean")
      gam <- metric_value("residual_run_order_gam_mean_deviance")
      batch <- metric_value("batch_weighted_pc_r2_categorical")
      runtime <- NA_real_
      coverage <- 100
    } else {
      row <- primary[primary$dataset_key == dataset_key & primary$method_id == method_id, ]
      if (nrow(row) != 1L) stop("Incomplete full panel: ", dataset_key, "/", method_id, call. = FALSE)
      cv <- row$heldout_qc_cv_mean
      gam <- row$residual_run_order_gam_mean_deviance
      batch <- row$batch_weighted_pc_r2_categorical
      runtime <- row$runtime_sec
      input_features <- primary$retained_features[primary$dataset_key == dataset_key & primary$method_id == "raw"]
      coverage <- 100 * row$retained_features / input_features
    }
    preservation <- preservation_display(dataset_key, method_id)
    if (is.numeric(preservation)) preservation <- fmt(preservation, 3L)
    data.frame(
      dataset_key = dataset_key, dataset = dataset_labels[[dataset_key]], method_id = method_id,
      method = full_labels[[method_id]], heldout_qc_cv_mean = cv,
      residual_run_order_gam_mean_deviance = gam,
      batch_weighted_pc_r2_categorical = batch, preservation = preservation,
      runtime_sec = runtime, feature_coverage_percent = coverage, stringsAsFactors = FALSE
    )
  }))
}))
write.csv(full_records, file.path(output_dir, "full_benchmark_tables.csv"), row.names = FALSE, na = "")
full_rows <- vapply(seq_len(nrow(full_records)), function(index) {
  row <- full_records[index, ]
  paste0(
    row$dataset, " & ", row$method, " & ", fmt(row$heldout_qc_cv_mean, 2L), " & ",
    fmt(row$residual_run_order_gam_mean_deviance, 3L), " & ",
    fmt(row$batch_weighted_pc_r2_categorical, 3L), " & ", row$preservation, " & ",
    fmt(row$runtime_sec, 1L), " & ", fmt(row$feature_coverage_percent, 1L), "\\% \\\\"
  )
}, character(1))
full_tex <- c(
  "\\begin{landscape}",
  "\\begin{longtable}{llrrrrrr}",
  "\\caption{Complete nine-method benchmark from the frozen WiNN 0.1.6 release. Simulation values are means across the ten completed frozen seeds; empirical values use the common frozen panel.}\\label{tab:supp_full_benchmark}\\\\",
  "\\toprule Dataset & Method & CV (\\%) & GAM dev. & Batch $R^2$ & Preservation & Runtime (s) & Coverage \\\\ \\midrule",
  "\\endfirsthead",
  "\\toprule Dataset & Method & CV (\\%) & GAM dev. & Batch $R^2$ & Preservation & Runtime (s) & Coverage \\\\ \\midrule",
  "\\endhead",
  full_rows,
  "\\bottomrule",
  "\\end{longtable}",
  "\\noindent Preservation is truth metabolite ICC(A,1) for simulation, cross-batch metabolite ICC(A,1) for MTBLS79, genuine replicate ICC(A,1) for BatchCorr and the clinical cohort, age/BMI/gender weighted-PC $R^2$ for Sacurine (shown in that order), and deposited-group weighted-PC $R^2$ for WaveICA. A dash marks an unrepeated fit time, not a zero-second runtime.",
  "\\end{landscape}"
)
write_lines(full_tex, "full_benchmark_tables.tex")

ablation_specs <- data.frame(
  table_id = c(
    "heldout_qc_cv", "gam_deviance", "batch_r2", "metabolite_icc_a1",
    "sample_pearson", "feature_repeatability_ratio"
  ),
  caption = c(
    "Held-out reference mean CV (percentage points; lower is better)",
    "Mean residual GAM deviance explained (lower is better)",
    "Categorical batch weighted-PC $R^2$ (lower is better)",
    "Metabolite ICC(A,1) (higher is better)",
    "Sample Pearson agreement (higher is better)",
    "Feature repeatability ratio (higher is better)"
  ), stringsAsFactors = FALSE
)
ablation_dataset_order <- c(
  "simulation", "mtbls79", "clinical_fiams", "batchcorr_set1", "sacurine", "waveica"
)
ablation_tex <- unlist(lapply(seq_len(nrow(ablation_specs)), function(spec_index) {
  specification <- ablation_specs[spec_index, ]
  selected <- ablation_impacts[ablation_impacts$table_id == specification$table_id, ]
  rows <- vapply(ablation_dataset_order, function(dataset_key) {
    values <- selected[selected$dataset_key == dataset_key, ]
    values <- values[order(values$stage_index), ]
    if (nrow(values) != 4L) stop("Incomplete ablation table: ", specification$table_id, "/", dataset_key, call. = FALSE)
    raw <- unique(values$cum_raw)
    if (length(raw) != 1L) stop("Ambiguous cumulative raw value in ablation table.", call. = FALSE)
    display <- c(raw, values$raw_change)
    paste0(values$dataset[1], " & ", paste(fmt(display, 3L), collapse = " & "), " \\\\")
  }, character(1))
  c(
    paste0(
      "\\begin{table*}[!htbp]\\centering\\scriptsize\\caption{", specification$caption,
      ". Values after Raw are consecutive changes in the fixed order; dashes denote an endpoint unsupported by design. The simulation row is one prespecified representative realization.}"
    ),
    "\\begin{tabular}{lrrrrr}\\toprule Dataset & Raw & Outlier $\\Delta$ & Drift $\\Delta$ & Batch $\\Delta$ & PQN $\\Delta$ \\\\ \\midrule",
    rows,
    "\\bottomrule\\end{tabular}\\end{table*}"
  )
}), use.names = FALSE)
write_lines(ablation_tex, "ablation_tables.tex")

batch_trajectory <- ablation_impacts[
  ablation_impacts$table_id == "batch_r2" & ablation_impacts$stage_index == 1L,
  c("dataset_key", "dataset", "cum_raw", "cum_outlier", "cum_drift", "cum_batch", "cum_pqn")
]
batch_trajectory <- batch_trajectory[match(ablation_dataset_order, batch_trajectory$dataset_key), ]
if (nrow(batch_trajectory) != length(ablation_dataset_order) || any(is.na(batch_trajectory$dataset_key))) {
  stop("Incomplete categorical batch-R2 cumulative trajectory.", call. = FALSE)
}
trajectory_rows <- vapply(seq_len(nrow(batch_trajectory)), function(index) {
  paste0(
    batch_trajectory$dataset[index], " & ",
    paste(fmt(unlist(batch_trajectory[index, c("cum_raw", "cum_outlier", "cum_drift", "cum_batch", "cum_pqn")]), 4L), collapse = " & "),
    " \\\\"
  )
}, character(1))
trajectory_tex <- c(
  "\\begin{table*}[!htbp]\\centering\\scriptsize\\caption{Factor-coded categorical batch weighted-PC $R^2$ through the frozen fixed cumulative workflow.}\\label{tab:supp_batch_trajectory}",
  "\\begin{tabular}{lrrrrr}\\toprule Dataset & Raw & Outlier & Selective drift & Residual batch & PQN \\\\ \\midrule",
  trajectory_rows,
  "\\bottomrule\\end{tabular}\\end{table*}"
)
write_lines(trajectory_tex, "supp_categorical_r2_trajectory.tex")
trajectory_source <- reshape(
  batch_trajectory,
  varying = c("cum_raw", "cum_outlier", "cum_drift", "cum_batch", "cum_pqn"),
  v.names = "value", timevar = "stage", times = c("Raw", "Outlier", "Selective drift", "Residual batch", "PQN"),
  direction = "long"
)
rownames(trajectory_source) <- NULL
write.csv(
  trajectory_source[c("dataset_key", "dataset", "stage", "value")],
  file.path(output_dir, "supp_categorical_r2_trajectory.csv"), row.names = FALSE, na = ""
)

reference_cv <- reference_split_metrics[
  reference_split_metrics$metric == "heldout_qc_cv_mean" &
    reference_split_metrics$method_id %in% c("raw", "winn_fixed_no_qc"),
  c("dataset_key", "dataset", "split_id", "method_id", "value")
]
reference_wide <- reshape(
  reference_cv, idvar = c("dataset_key", "dataset", "split_id"),
  timevar = "method_id", direction = "wide"
)
reference_wide$fixed_minus_raw <- reference_wide$value.winn_fixed_no_qc - reference_wide$value.raw
reference_summary <- do.call(rbind, lapply(split(reference_wide, reference_wide$dataset_key), function(values) {
  data.frame(
    dataset_key = values$dataset_key[1], dataset = values$dataset[1], splits = nrow(values),
    mean_cv_change = mean(values$fixed_minus_raw), sd_cv_change = stats::sd(values$fixed_minus_raw),
    stringsAsFactors = FALSE
  )
}))
reference_summary <- reference_summary[match(
  c("simulation", "mtbls79", "clinical_fiams", "batchcorr_set1", "sacurine", "waveica"),
  reference_summary$dataset_key
), ]
reference_rows <- vapply(seq_len(nrow(reference_summary)), function(index) paste0(
  reference_summary$dataset[index], " & ", reference_summary$splits[index], " & ",
  fmt(reference_summary$mean_cv_change[index], 2L), " & ", fmt(reference_summary$sd_cv_change[index], 2L), " \\\\"
), character(1))
write_lines(c(
  "\\begin{table*}[!htbp]\\centering\\scriptsize\\caption{Fixed WiNN sensitivity across 10 overlapping reference splits. Negative CV change favors WiNN over raw.}\\label{tab:supp_reference_splits}",
  "\\begin{tabular}{lrrr}\\toprule Dataset & Splits & Mean CV change & SD \\\\ \\midrule",
  reference_rows, "\\bottomrule\\end{tabular}\\end{table*}"
), "supp_reference_split_stability.tex")
write.csv(reference_summary, file.path(output_dir, "supp_reference_split_stability.csv"), row.names = FALSE, na = "")

confounding_metrics <- c(
  "heldout_qc_cv_mean", "batch_weighted_pc_r2_categorical",
  "median_attenuation_ratio_injected_responsive", "true_positive_rate", "empirical_fdr"
)
confounding_selected <- partial_confounding[
  partial_confounding$method == "WiNN fixed default (no QC)" &
    partial_confounding$metric %in% confounding_metrics &
    ((partial_confounding$confounding_type == "none" & partial_confounding$level == 0) |
      (partial_confounding$confounding_type != "none" & partial_confounding$level %in% c(0.5, 0.75, 0.9, 1))),
  c("confounding_type", "level", "metric", "mean")
]
confounding_wide <- reshape(
  confounding_selected, idvar = c("confounding_type", "level"),
  timevar = "metric", direction = "wide"
)
confounding_order <- order(
  match(confounding_wide$confounding_type, c("batch", "joint", "none", "run_order")),
  confounding_wide$level
)
confounding_wide <- confounding_wide[confounding_order, ]
confounding_rows <- vapply(seq_len(nrow(confounding_wide)), function(index) paste0(
  gsub("_", "\\\\_", confounding_wide$confounding_type[index]), " & ", fmt(confounding_wide$level[index], 2L), " & ",
  fmt(confounding_wide$mean.heldout_qc_cv_mean[index], 2L), " & ",
  fmt(confounding_wide$mean.batch_weighted_pc_r2_categorical[index], 3L), " & ",
  fmt(confounding_wide$mean.median_attenuation_ratio_injected_responsive[index], 3L), " & ",
  fmt(confounding_wide$mean.true_positive_rate[index], 3L), " & ",
  fmt(confounding_wide$mean.empirical_fdr[index], 3L), " \\\\"
), character(1))
write_lines(c(
  "\\begin{table*}[!htbp]\\centering\\scriptsize\\caption{WiNN fixed partial-confounding boundary (10-seed means at selected scenarios).}\\label{tab:supp_confounding}",
  "\\begin{tabular}{llrrrrr}\\toprule Alignment & Strength & CV (\\%) & Batch $R^2$ & Attenuation & TPR & FDR \\\\ \\midrule",
  confounding_rows, "\\bottomrule\\end{tabular}\\end{table*}"
), "supp_partial_confounding.tex")
write.csv(confounding_wide, file.path(output_dir, "supp_partial_confounding.csv"), row.names = FALSE, na = "")

association_keys <- merge(
  associations[associations$method_id %in% c("raw", "winn_fixed_no_qc"), ],
  permutation[permutation$method_id %in% c("raw", "winn_fixed_no_qc"),
    c("dataset_key", "method_id", "endpoint", "null_median", "null_95th_percentile", "empirical_p_value")],
  by = c("dataset_key", "method_id", "endpoint"), all.x = TRUE, sort = FALSE
)
write.csv(association_keys, file.path(output_dir, "manuscript_association_calibration.csv"), row.names = FALSE, na = "")

lobo_summary <- aggregate(
  cbind(bootstrap_mean, ci_lower, ci_upper) ~ dataset_key + dataset + endpoint + method_id + method + metric,
  data = lobo[lobo$method_id %in% c("raw", "winn_fixed_no_qc"), ],
  FUN = stats::median, na.rm = TRUE
)
names(lobo_summary)[names(lobo_summary) == "bootstrap_mean"] <- "median_fold_bootstrap_mean"
names(lobo_summary)[names(lobo_summary) == "ci_lower"] <- "median_fold_ci_lower"
names(lobo_summary)[names(lobo_summary) == "ci_upper"] <- "median_fold_ci_upper"
write.csv(lobo_summary, file.path(output_dir, "manuscript_lobo_summary.csv"), row.names = FALSE, na = "")

manifest <- data.frame(
  artifact = list.files(output_dir, full.names = FALSE),
  generated_from_frozen_results = TRUE,
  correction_methods_executed = FALSE,
  tuning_executed = FALSE,
  selection_executed = FALSE,
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(output_dir, "manuscript_summary_manifest.csv"), row.names = FALSE)

# Refresh both inventories after writing manuscript-facing derivatives. The
# manifest files exclude themselves, so their contents are deterministic for a
# fixed set of artifacts.
derivative_files <- list.files(derivative_root, recursive = TRUE, full.names = TRUE)
derivative_files <- derivative_files[
  file.info(derivative_files)$isdir %in% FALSE & basename(derivative_files) != "derivative_file_manifest.csv"
]
derivative_manifest <- data.frame(
  path = substring(derivative_files, nchar(derivative_root) + 2L),
  bytes = as.numeric(file.info(derivative_files)$size),
  sha256 = vapply(derivative_files, digest::digest, character(1), file = TRUE, algo = "sha256"),
  stringsAsFactors = FALSE
)
write.csv(derivative_manifest, file.path(derivative_root, "derivative_file_manifest.csv"), row.names = FALSE)

final_files <- list.files(final_root, recursive = TRUE, full.names = TRUE)
final_files <- final_files[
  file.info(final_files)$isdir %in% FALSE & basename(final_files) != "file_manifest.csv"
]
final_manifest <- data.frame(
  path = substring(final_files, nchar(final_root) + 2L),
  bytes = as.numeric(file.info(final_files)$size),
  sha256 = vapply(final_files, digest::digest, character(1), file = TRUE, algo = "sha256"),
  stringsAsFactors = FALSE
)
write.csv(final_manifest, file.path(final_root, "file_manifest.csv"), row.names = FALSE)
cat("Wrote manuscript-facing derivative summaries without executing correction code.\n")
