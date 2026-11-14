#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
source(file.path(repo_root, "analysis", "ablation_stage_helpers.R"), local = TRUE)

source_paths <- list.files(
  file.path(repo_root, "results", "final", "tables"),
  pattern = "^ablation_cumulative_values_.*[.]csv$", full.names = TRUE
)
if (!length(source_paths)) stop("No cumulative ablation tables were found.", call. = FALSE)
cumulative <- do.call(rbind, lapply(source_paths, function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}))
rownames(cumulative) <- NULL
cumulative <- normalize_ablation_cumulative_metric_names(cumulative)

impacts <- derive_ablation_stage_impacts(cumulative)
validate_ablation_stage_impacts(cumulative, impacts)

expected_rows <- nrow(cumulative) * 4L
stopifnot(nrow(impacts) == expected_rows)
stopifnot(all(!is.na(impacts$metric) & nzchar(impacts$metric)))
stopifnot(!anyDuplicated(impacts[c("dataset_key", "table_id", "metric", "stage")]))
stopifnot(identical(
  sort(unique(impacts$stage)),
  sort(ablation_stage_spec()$stage)
))

first <- impacts[impacts$stage_index == 1L, , drop = FALSE]
stopifnot(all.equal(
  first$raw_change,
  first$cum_outlier - first$cum_raw,
  tolerance = 1e-14, check.attributes = FALSE
))

identical_pair <- data.frame(
  dataset_key = "test", dataset = "Test", table_id = "test", metric = "test_metric",
  direction = "lower", Raw = 10, Outlier_shrinkage = 8, Selective_drift = 8,
  Selective_batch = 8, PQN = 8, stringsAsFactors = FALSE
)
identical_impacts <- derive_ablation_stage_impacts(identical_pair)
stopifnot(identical_impacts$raw_change[identical_impacts$stage != "Outlier_shrinkage"] == 0)
stopifnot(identical_impacts$direction_normalized_benefit[1] == 2)

descriptive <- identical_pair
descriptive$direction <- "descriptive"
stopifnot(all(is.na(derive_ablation_stage_impacts(descriptive)$direction_normalized_benefit)))

simulation_qc <- impacts[
  impacts$dataset_key == "simulation" & impacts$table_id == "heldout_qc_cv",
  , drop = FALSE
]
# Check the current SIM01 trajectory against the separately saved evaluator
# output, rather than requiring a new gate to reproduce old scientific values.
evaluated <- read.csv(file.path(repo_root, "results", "evaluation", "ablations",
                                "simulation", "method_metrics.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
variant_ids <- c("C0_RAW", "C1_OUTLIER", "C2_SELECTIVE_DRIFT", "C3_SELECTIVE_BATCH", "C4_FULL_FIXED")
evaluated <- evaluated[evaluated$metric == "heldout_qc_cv_mean" &
                       evaluated$method_id %in% variant_ids, , drop = FALSE]
stopifnot(nrow(evaluated) == 5L, !anyDuplicated(evaluated$method_id))
current_values <- evaluated$value[match(variant_ids, evaluated$method_id)]
stopifnot(all(is.finite(current_values)))
expected_benefit <- -diff(current_values)
stopifnot(all.equal(
  simulation_qc$direction_normalized_benefit,
  expected_benefit, tolerance = 1e-12, check.attributes = FALSE
))
# The input and MAD preparation are unchanged by drift recalibration.
stopifnot(abs(simulation_qc$direction_normalized_benefit[1] - 2.9709103854638) < 1e-12)

source(file.path(repo_root, "scripts", "weighted_pc_r2.R"), local = TRUE)
set.seed(42)
batch_matrix <- matrix(rexp(120), nrow = 10)
numeric_batch_error <- try(
  weighted_pc_r2_explicit(batch_matrix, rep(1:3, each = 4)),
  silent = TRUE
)
categorical_batch <- weighted_pc_r2_explicit(
  batch_matrix, factor(rep(1:3, each = 4)), target_type = "categorical"
)
stopifnot(inherits(numeric_batch_error, "try-error"), is.finite(categorical_batch))

saved_path <- file.path(repo_root, "results", "final", "tables", "ablation_stage_impacts.csv")
if (file.exists(saved_path)) {
  saved <- read.csv(saved_path, stringsAsFactors = FALSE, check.names = FALSE)
  validate_ablation_stage_impacts(cumulative, saved)
}

figure_path <- file.path(repo_root, "results", "final", "figures", "source_data", "ablation_stage_impacts.csv")
if (file.exists(saved_path) && file.exists(figure_path)) {
  saved <- read.csv(saved_path, stringsAsFactors = FALSE, check.names = FALSE)
  figure <- read.csv(figure_path, stringsAsFactors = FALSE, check.names = FALSE)
  plotted <- saved[saved$table_id %in% c("heldout_qc_cv", "gam_deviance", "batch_r2"), , drop = FALSE]
  key <- c("dataset_key", "table_id", "metric", "stage")
  plotted <- plotted[do.call(order, plotted[key]), , drop = FALSE]
  figure <- figure[do.call(order, figure[key]), , drop = FALSE]
  core <- c(key, "cumulative_previous", "cumulative_current", "raw_change", "direction_normalized_benefit")
  stopifnot(isTRUE(all.equal(
    plotted[core], figure[core],
    tolerance = 0, check.attributes = FALSE
  )))
}

cat("Ablation stage-impact and categorical batch-R2 tests passed.\n")
