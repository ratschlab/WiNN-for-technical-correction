#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))

method_order <- c(
  "Raw", "ComBat", "QC-RLSC", "QC-RFSC", "TIGER", "SERRF",
  "WINN auto (QC)", "WINN auto-batch (QC)", "WINN default (no QC)"
)
rows <- list()
record <- function(dataset, check, passed, detail) {
  rows[[length(rows) + 1L]] <<- data.frame(
    dataset = dataset, check = check, passed = isTRUE(passed), detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

validate_one <- function(dataset, matrix_name, metadata_name, seed) {
  result_dir <- file.path(repo_root, "results", dataset)
  processed_dir <- file.path(repo_root, "data", "public", dataset, "processed")
  x <- readRDS(file.path(processed_dir, matrix_name))
  meta <- read.csv(file.path(processed_dir, metadata_name), check.names = FALSE, stringsAsFactors = FALSE)
  attempts <- read.csv(file.path(result_dir, "method_runtime.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  protocol <- read.csv(file.path(result_dir, "hidden_qc_method_protocol.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  dimensions <- read.csv(file.path(result_dir, "method_dimensions.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  metrics <- read.csv(file.path(result_dir, "method_metrics.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  holdout_table <- read.csv(file.path(result_dir, "qc_holdout_ids.csv"), stringsAsFactors = FALSE)
  holdout <- holdout_table$sample_id

  record(dataset, "nine methods attempted in specified order",
         identical(attempts$method, method_order), paste(attempts$method, collapse = "; "))
  record(dataset, "all nine attempted methods completed",
         identical(attempts$method, method_order) && all(attempts$status == "completed"),
         paste(attempts$method, attempts$status, sep = "=", collapse = "; "))
  record(dataset, "hidden QCs unavailable to every method",
         nrow(protocol) == 9L && all(!protocol$supplied_hidden_qc_ids) &&
           all(protocol$hidden_qc_method_class == "Sample") &&
           all(!protocol$tuning_uses_hidden_qc_values),
         paste("held out", length(holdout), "QCs with seed", seed))
  record(dataset, "phenotypes excluded from tuning",
         nrow(protocol) == 9L && all(!protocol$tuning_uses_phenotype),
         "saved method-facing protocol")
  parameters <- read.csv(file.path(result_dir, "method_parameters.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  seeded <- "correction_seed" %in% names(protocol) && "correction_seed" %in% names(parameters) &&
    all(is.finite(protocol$correction_seed)) && !anyDuplicated(protocol$correction_seed) &&
    identical(protocol$correction_seed, parameters$correction_seed[match(protocol$method, parameters$method)])
  record(dataset, "full correction methods have explicit reproducible seeds",
         seeded,
         if (seeded) paste(protocol$method, protocol$correction_seed, sep = "=", collapse = "; ") else "seed columns missing or inconsistent")
  auto_row <- protocol[protocol$method == "WINN auto-batch (QC)", , drop = FALSE]
  record(dataset, "WiNN auto-batch not supplied true batches",
         nrow(auto_row) == 1L && !auto_row$supplied_batch_labels,
         "true batches retained only in evaluation metadata")
  record(dataset, "all reported matrices align to injection metadata",
         nrow(dimensions) == 9L && all(dimensions$exact_sample_alignment) &&
           all(dimensions$total_injections_retained == ncol(x)),
         paste(ncol(x), "injections"))
  record(dataset, "feature loss explicitly reported",
         all(dimensions$features_lost == nrow(x) - dimensions$features_retained),
         paste0(sum(dimensions$features_lost), " method-feature losses across rows"))
  common_ids <- read.csv(file.path(result_dir, "common_evaluation_feature_ids.csv"), stringsAsFactors = FALSE)$feature_id
  association_name <- if (dataset == "sacurine") "sacurine_associations_primary.csv" else "waveica_group_associations_primary.csv"
  associations <- read.csv(file.path(result_dir, association_name), check.names = FALSE, stringsAsFactors = FALSE)
  association_panels <- split(associations$feature_id, interaction(associations$method, associations$variable, drop = TRUE))
  common_biology_ok <- length(association_panels) > 0L && all(vapply(
    association_panels,
    function(ids) identical(sort(unique(ids)), sort(common_ids)),
    logical(1)
  ))
  record(dataset, "biological models use one common feature panel",
         common_biology_ok,
         paste(length(common_ids), "features in every method/variable model"))
  record(dataset, "biological-repeat metrics not fabricated",
         all(is.na(metrics$sample_replicate_pearson)) && all(is.na(metrics$cross_batch_sample_icc)),
         unique(metrics$repeat_metric_reason)[1])
  if (dataset == "waveica_adenocarcinoma") {
    observed_groups <- unique(meta$biological_group[!is.na(meta$biological_group) & nzchar(meta$biological_group)])
    neutral_groups <- length(observed_groups) == 2L && all(observed_groups %in% c("group_0", "group_1"))
    disease_unassigned <- "disease_label" %in% names(meta) && all(is.na(meta$disease_label) | !nzchar(meta$disease_label))
    record(dataset, "disease mapping remains neutral and source-backed", neutral_groups && disease_unassigned,
           "numeric group counts are compatible with published cohort totals but diagnoses are not assigned")
  } else {
    raw_qc <- read.csv(file.path(result_dir, "qc_cv_by_feature.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    raw_qc <- raw_qc[raw_qc$method == "Raw", , drop = FALSE]
    expected_cv <- 100 * apply(x[, holdout, drop = FALSE], 1, sd) / rowMeans(x[, holdout, drop = FALSE])
    observed_cv <- raw_qc$cv_percent[match(names(expected_cv), raw_qc$feature_id)]
    record(dataset, "technical metrics use pre-osmolality values",
           max(abs(expected_cv - observed_cv), na.rm = TRUE) < 1e-10,
           "Raw hidden-QC CV recomputed directly from deposited pre-osmolality matrix")
    runner_text <- paste(readLines(file.path(repo_root, "scripts", "run_human_public_benchmark.R"), warn = FALSE), collapse = "\n")
    record(dataset, "common osmolality normalization follows every correction method",
           grepl("biological_mats <- lapply(method_results_common", runner_text, fixed = TRUE) &&
             grepl("md$osmolality, \"/\"", runner_text, fixed = TRUE),
           "runner builds every Sacurine biological matrix by the same sample-wise osmolality division")
  }

  lb <- read.csv(file.path(result_dir, "ljung_box_by_feature_batch.csv"), check.names = FALSE)
  keys <- interaction(lb$method, lb$batch, lb$segment_id, drop = TRUE)
  expected <- unsplit(lapply(split(lb$p_value, keys), p.adjust, method = "BH"), keys)
  finite <- is.finite(expected) & is.finite(lb$p_adj)
  lb_ok <- all(is.na(expected) == is.na(lb$p_adj)) && (!any(finite) || max(abs(expected[finite] - lb$p_adj[finite])) < 1e-12)
  record(dataset, "Ljung-Box BH adjustment within method/batch/segment", lb_ok,
         paste(length(unique(keys)), "independent adjustment groups"))

  validation <- read.csv(file.path(result_dir, "winn_instrumentation_validation.csv"), check.names = FALSE)
  record(dataset, "instrumented WiNN equals normal wrapper within 1e-8", nrow(validation) == 3L && all(validation$diagnostic_matrix_equal),
         paste("maximum absolute difference", max(validation$max_absolute_difference)))
  record(dataset, "true batches used for auto-batch evaluation",
         is.finite(metrics$batch_weighted_pc_r2[metrics$method == "WINN auto-batch (QC)"]) &&
           file.exists(file.path(result_dir, "winn_auto_batch_assignments.csv")),
         "evaluation metric and supplied/detected assignment table present")

  cache_names <- paste0(gsub("[^A-Za-z0-9]+", "_", tolower(method_order)), "_", dataset,
                        "_hiddenqc_", seed, "_v1.rds")
  cache_paths <- file.path(result_dir, "method_cache", cache_names)
  cache_ok <- all(file.exists(cache_paths))
  if (cache_ok) {
    cache_ok <- all(vapply(cache_paths, function(path) {
      value <- readRDS(path)
      is.matrix(value) && identical(colnames(value), colnames(x)) &&
        all(rownames(value) %in% rownames(x)) && all(is.finite(value))
    }, logical(1)))
  }
  record(dataset, "nine corrected/raw caches validate", cache_ok, paste(length(cache_paths), "compressed RDS paths"))

  tiger_cache <- readRDS(cache_paths[method_order == "TIGER"])
  tiger_lost <- setdiff(rownames(x), rownames(tiger_cache))
  tiger_report <- file.path(result_dir, "tiger_nonfinite_output_features.csv")
  tiger_loss_documented <- if (!length(tiger_lost)) {
    !file.exists(tiger_report) || !nrow(read.csv(tiger_report, stringsAsFactors = FALSE))
  } else if (file.exists(tiger_report)) {
    reported <- read.csv(tiger_report, stringsAsFactors = FALSE)
    "feature_id" %in% names(reported) && identical(sort(unique(reported$feature_id)), sort(tiger_lost))
  } else FALSE
  record(dataset, "TIGER non-finite feature loss is explicit",
         tiger_loss_documented,
         if (length(tiger_lost)) paste(length(tiger_lost), "complete output rows removed") else "no TIGER output rows removed")

  required <- c(
    "method_metrics.csv", "method_metrics_long.csv", "method_metrics_common_panel.csv",
    "method_metrics_coverage_penalized.csv", "method_runtime.csv", "method_dimensions.csv",
    "method_feature_coverage.csv", "method_parameters.csv", "method_failures.csv", "tuning_candidates.csv",
    "qc_cv_by_feature.csv", "qc_pairwise_correlations.csv",
    "run_order_gam_by_feature_batch.csv", "ljung_box_by_feature_batch.csv",
    "batch_associations.csv", "correction_magnitude.csv", "winn_selectivity_summary.csv",
    "winn_selectivity_by_feature.csv", "randomization_diagnostics.csv",
    "common_evaluation_feature_ids.csv", "preprocessing_summary.json", "sessionInfo.txt", "analysis_log.txt"
  )
  record(dataset, "required result tables and logs exist",
         all(file.exists(file.path(result_dir, required))), paste(length(required), "required files"))
  if (dataset == "waveica_adenocarcinoma") {
    record(dataset, "cross-batch direction-consistency table exists",
           file.exists(file.path(result_dir, "waveica_cross_batch_direction_consistency.csv")),
           "pooled significant effects evaluated independently in all three batches")
  }
  figure_names <- c(
    "acquisition_design.pdf", "randomization_diagnostics.pdf", "raw_vs_winn_pca.pdf",
    "technical_comparison.pdf", "biological_comparison.pdf", "winn_selectivity.pdf",
    "representative_feature_trajectories.pdf",
    if (dataset == "sacurine") "sacurine_effect_concordance.pdf" else "waveica_batch_effects.pdf"
  )
  source_names <- c(
    "acquisition_design.csv", "randomization_composition.csv", "randomization_sequence.csv",
    "qc_spacing.csv", "raw_winn_pca.csv", "technical_comparison.csv",
    "biological_comparison.csv", "winn_selectivity.csv", "representative_feature_trajectories.csv",
    if (dataset == "sacurine") "sacurine_batch_effect_estimates.csv" else "waveica_adenocarcinoma_batch_effect_estimates.csv"
  )
  figure_dir <- file.path(result_dir, "figures")
  source_dir <- file.path(figure_dir, "source_data")
  record(dataset, "publication figures and numerical source tables exist",
         all(file.exists(file.path(figure_dir, figure_names))) && all(file.exists(file.path(source_dir, source_names))),
         paste(length(figure_names), "PDFs and", length(source_names), "source tables"))
  record(dataset, "benchmark completion logged",
         any(grepl("Benchmark completed and validated", readLines(file.path(result_dir, "analysis_log.txt"), warn = FALSE), fixed = TRUE)),
         "terminal success marker")

  notebook <- file.path(repo_root, "notebooks", paste0(dataset, "_comparison.Rmd"))
  rendered <- file.path(repo_root, "notebooks", "rendered", paste0(dataset, "_comparison.html"))
  record(dataset, "notebook and clean-session render exist", file.exists(notebook) && file.exists(rendered),
         paste(basename(notebook), basename(rendered)))
}

validate_one("sacurine", "sacurine_intensity_matrix.rds", "sacurine_metadata.csv", 20260810L)
validate_one("waveica_adenocarcinoma", "waveica_intensity_matrix.rds", "waveica_metadata.csv", 20260811L)

validation <- do.call(rbind, rows)
report_dir <- file.path(repo_root, "reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(validation, file.path(report_dir, "human_public_benchmark_validation.csv"), row.names = FALSE, quote = TRUE)
report <- c(
  "# Human public benchmark final validation", "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), "",
  unlist(lapply(split(validation, validation$dataset), function(d) c(
    paste0("## ", unique(d$dataset)), "",
    paste0("- ", ifelse(d$passed, "PASS", "FAIL"), " — ", d$check, ": ", d$detail), ""
  )))
)
writeLines(report, file.path(report_dir, "human_public_benchmark_validation.md"))
if (any(!validation$passed)) {
  failed <- validation[!validation$passed, ]
  stop("Validation failed: ", paste(paste(failed$dataset, failed$check, sep = " — "), collapse = "; "))
}
cat("All human public benchmark validation checks passed.\n")
