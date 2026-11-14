#!/usr/bin/env Rscript

# Paired leave-one-batch-out uncertainty from frozen matrices only.

args <- commandArgs(trailingOnly = TRUE)
dataset_key <- sub("^--dataset=", "", args[grep("^--dataset=", args)])
if (length(dataset_key) != 1L || !dataset_key %in% c("waveica", "mtbls79", "sacurine")) {
  stop("Supply --dataset=waveica, --dataset=mtbls79, or --dataset=sacurine.", call. = FALSE)
}
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
for (package in c("digest", "jsonlite", "matrixStats")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package, call. = FALSE)
}
source(file.path(repo_root, "analysis", "derivative_biology_helpers.R"), local = TRUE)
config <- jsonlite::read_json(
  file.path(repo_root, "analysis", "config", "derivative_biology_config.json"), simplifyVector = TRUE
)
n_bootstrap <- as.integer(config$common_processing$bootstrap_replicates)
test_mode <- identical(Sys.getenv("DERIVATIVE_TEST_MODE", "false"), "true")
if (test_mode) n_bootstrap <- as.integer(Sys.getenv("DERIVATIVE_TEST_BOOTSTRAP", "10"))
cores <- max(1L, as.integer(Sys.getenv("DERIVATIVE_CORES", "1")))
if (.Platform$OS.type == "windows") cores <- 1L
method_registry <- derivative_method_registry()
# Partial re-run: see run_derivative_biology.R. Each LOBO draw is seeded per task
# before the method loop, so spliced rows share every draw with recomputed rows.
full_method_registry <- method_registry
compute_methods <- trimws(strsplit(Sys.getenv("DERIVATIVE_COMPUTE_METHODS", ""), ",", fixed = TRUE)[[1]])
compute_methods <- compute_methods[nzchar(compute_methods)]
partial_mode <- length(compute_methods) > 0L
previous_root <- Sys.getenv("DERIVATIVE_PREVIOUS_ROOT", "")
if (partial_mode) {
  if (!all(compute_methods %in% method_registry$method_id)) stop("Unknown DERIVATIVE_COMPUTE_METHODS.", call. = FALSE)
  method_registry <- method_registry[method_registry$method_id %in% compute_methods, , drop = FALSE]
}
pre_manifest <- read.csv(
  file.path(repo_root, "results", "frozen_endpoint_free_inputs", "provenance", "matrix_checksums_pre.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
frozen <- derivative_load_matrices(repo_root, dataset_key, pre_manifest)
meta <- derivative_load_metadata(repo_root, dataset_key, frozen$sample_ids)
study <- which(meta$is_study)
md <- meta[study, , drop = FALSE]
log_matrices <- lapply(frozen$matrices, function(value) log1p(value[, study, drop = FALSE]))
output_dir <- file.path(repo_root, "results", "final", "derivative", dataset_key)

bind_rows_base <- function(values) {
  values <- Filter(function(value) !is.null(value) && nrow(value), values)
  if (!length(values)) return(data.frame())
  do.call(rbind, values)
}

sample_strata_indices <- function(strata) {
  unlist(lapply(unique(strata), function(stratum) {
    index <- which(strata == stratum)
    sample(index, length(index), replace = TRUE)
  }), use.names = FALSE)
}

summarise_lobo <- function(replicates) {
  keys <- unique(replicates[c("heldout_batch", "endpoint", "method_id", "method", "metric")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(index) {
    key <- keys[index, , drop = FALSE]
    selected <- replicates[
      replicates$heldout_batch == key$heldout_batch &
        replicates$endpoint == key$endpoint &
        replicates$method_id == key$method_id &
        replicates$metric == key$metric,
      "value"
    ]
    data.frame(
      key, bootstrap_mean = mean(selected, na.rm = TRUE),
      ci_lower = stats::quantile(selected, 0.025, na.rm = TRUE, names = FALSE),
      ci_upper = stats::quantile(selected, 0.975, na.rm = TRUE, names = FALSE),
      valid_bootstrap = sum(is.finite(selected)), stringsAsFactors = FALSE
    )
  }))
}

lobo_seed <- as.integer(config$common_processing$bootstrap_seed) + 9100L +
  match(dataset_key, c("waveica", "mtbls79", "sacurine"))

if (dataset_key == "waveica") {
  group <- factor(md$biological_group, levels = c("group_0", "group_1"))
  folds <- levels(factor(md$batch))
  worker <- function(task_index) {
    set.seed(lobo_seed + task_index)
    fold_index <- ((task_index - 1L) %% length(folds)) + 1L
    bootstrap_index <- ((task_index - 1L) %/% length(folds)) + 1L
    heldout <- folds[fold_index]
    training_original <- which(md$batch != heldout)
    testing_original <- which(md$batch == heldout)
    training_draw <- training_original[sample_strata_indices(interaction(md$batch[training_original], group[training_original], drop = TRUE))]
    testing_draw <- testing_original[sample_strata_indices(group[testing_original])]
    bind_rows_base(lapply(seq_len(nrow(method_registry)), function(method_index) {
      method_id <- method_registry$method_id[method_index]
      training_fit <- two_group_model(log_matrices[[method_id]][, training_draw, drop = FALSE], group[training_draw], "group_0", "group_1")
      testing_fit <- two_group_model(log_matrices[[method_id]][, testing_draw, drop = FALSE], group[testing_draw], "group_0", "group_1")
      if (is.null(training_fit) || is.null(testing_fit)) return(NULL)
      values <- effect_concordance(training_fit$effect, testing_fit$effect)
      data.frame(
        bootstrap = bootstrap_index, heldout_batch = heldout, endpoint = "biological_group",
        method_id = method_id, method = method_registry$method[method_index],
        metric = names(values)[1:4], value = as.numeric(values[1:4]), stringsAsFactors = FALSE
      )
    }))
  }
  tasks <- seq_len(n_bootstrap * length(folds))
  replicates <- bind_rows_base(parallel::mclapply(tasks, worker, mc.cores = cores, mc.set.seed = TRUE, mc.preschedule = TRUE))
} else if (dataset_key == "mtbls79") {
  group <- as.character(md$class)
  unit <- as.character(md$sample_rep)
  folds <- unique(as.character(md$batch))
  fold_data <- lapply(folds, function(heldout) {
    training <- md$batch != heldout; testing <- md$batch == heldout
    common_units <- intersect(unique(unit[training]), unique(unit[testing]))
    classes <- vapply(common_units, function(value) unique(group[unit == value])[1], character(1))
    list(
      heldout = heldout, units = common_units, classes = classes,
      training_columns = lapply(common_units, function(value) which(training & unit == value)),
      testing_columns = lapply(common_units, function(value) which(testing & unit == value))
    )
  })
  worker <- function(task_index) {
    set.seed(lobo_seed + task_index)
    fold_index <- ((task_index - 1L) %% length(folds)) + 1L
    bootstrap_index <- ((task_index - 1L) %/% length(folds)) + 1L
    fold <- fold_data[[fold_index]]
    draw_index <- sample_strata_indices(fold$classes)
    bind_rows_base(lapply(seq_len(nrow(method_registry)), function(method_index) {
      method_id <- method_registry$method_id[method_index]
      training_matrix <- vapply(draw_index, function(index) {
        rowMeans(log_matrices[[method_id]][, fold$training_columns[[index]], drop = FALSE])
      }, numeric(nrow(log_matrices[[method_id]])))
      testing_matrix <- vapply(draw_index, function(index) {
        rowMeans(log_matrices[[method_id]][, fold$testing_columns[[index]], drop = FALSE])
      }, numeric(nrow(log_matrices[[method_id]])))
      classes <- fold$classes[draw_index]
      training_fit <- two_group_model(training_matrix, classes, "S", "C")
      testing_fit <- two_group_model(testing_matrix, classes, "S", "C")
      if (is.null(training_fit) || is.null(testing_fit)) return(NULL)
      values <- effect_concordance(training_fit$effect, testing_fit$effect)
      data.frame(
        bootstrap = bootstrap_index, heldout_batch = fold$heldout, endpoint = "cow_vs_sheep",
        method_id = method_id, method = method_registry$method[method_index],
        metric = names(values)[1:4], value = as.numeric(values[1:4]), stringsAsFactors = FALSE
      )
    }))
  }
  tasks <- seq_len(n_bootstrap * length(folds))
  replicates <- bind_rows_base(parallel::mclapply(tasks, worker, mc.cores = cores, mc.set.seed = TRUE, mc.preschedule = TRUE))
} else {
  md$age_z <- as.numeric(scale(md$age)); md$bmi_z <- as.numeric(scale(md$bmi))
  md$gender <- factor(md$gender); md$sampling_factor <- factor(md$sampling)
  biological_log <- lapply(frozen$matrices, function(value) {
    log1p(sweep(value[, study, drop = FALSE], 2L, md$osmolality, "/"))
  })
  folds <- unique(as.character(md$batch))
  endpoints <- c("age", "bmi", "gender")
  worker <- function(task_index) {
    set.seed(lobo_seed + task_index)
    fold_index <- ((task_index - 1L) %% length(folds)) + 1L
    bootstrap_index <- ((task_index - 1L) %/% length(folds)) + 1L
    heldout <- folds[fold_index]
    training_original <- which(md$batch != heldout)
    testing_original <- which(md$batch == heldout)
    training_draw <- training_original[sample_strata_indices(md$gender[training_original])]
    testing_draw <- testing_original[sample_strata_indices(md$gender[testing_original])]
    training_md <- droplevels(md[training_draw, , drop = FALSE])
    testing_md <- droplevels(md[testing_draw, , drop = FALSE])
    training_design <- tryCatch(stats::model.matrix(~ age_z + bmi_z + gender + sampling_factor, data = training_md), error = function(e) NULL)
    testing_design <- tryCatch(stats::model.matrix(~ age_z + bmi_z + gender + sampling_factor, data = testing_md), error = function(e) NULL)
    if (is.null(training_design) || is.null(testing_design)) return(NULL)
    bind_rows_base(lapply(seq_len(nrow(method_registry)), function(method_index) {
      method_id <- method_registry$method_id[method_index]
      bind_rows_base(lapply(endpoints, function(endpoint) {
        coefficient <- switch(endpoint,
          age = "age_z", bmi = "bmi_z", gender = grep("^gender", colnames(training_design), value = TRUE)[1]
        )
        testing_coefficient <- if (endpoint == "gender") grep("^gender", colnames(testing_design), value = TRUE)[1] else coefficient
        if (!length(coefficient) || !length(testing_coefficient) ||
            !coefficient %in% colnames(training_design) || !testing_coefficient %in% colnames(testing_design)) return(NULL)
        training_fit <- ols_feature_model(
          biological_log[[method_id]][, training_draw, drop = FALSE], training_design,
          match(coefficient, colnames(training_design))
        )
        testing_fit <- ols_feature_model(
          biological_log[[method_id]][, testing_draw, drop = FALSE], testing_design,
          match(testing_coefficient, colnames(testing_design))
        )
        if (is.null(training_fit) || is.null(testing_fit)) return(NULL)
        values <- effect_concordance(training_fit$table$coefficient, testing_fit$table$coefficient)
        data.frame(
          bootstrap = bootstrap_index, heldout_batch = heldout, endpoint = endpoint,
          method_id = method_id, method = method_registry$method[method_index],
          metric = names(values)[1:4], value = as.numeric(values[1:4]), stringsAsFactors = FALSE
        )
      }))
    }))
  }
  tasks <- seq_len(n_bootstrap * length(folds))
  replicates <- bind_rows_base(parallel::mclapply(tasks, worker, mc.cores = cores, mc.set.seed = TRUE, mc.preschedule = TRUE))
}

if (partial_mode) {
  previous_path <- file.path(previous_root, dataset_key, "leave_one_batch_out_bootstrap_replicates.csv")
  previous <- read.csv(previous_path, stringsAsFactors = FALSE, check.names = FALSE)
  previous <- previous[!previous$method_id %in% compute_methods, , drop = FALSE]
  if (!nrow(previous) || !setequal(names(previous), names(replicates))) stop("Previous LOBO replicates are not reusable.", call. = FALSE)
  replicates <- rbind(previous, replicates[names(previous)])
  replicates <- replicates[order(match(replicates$method_id, full_method_registry$method_id), seq_len(nrow(replicates))), , drop = FALSE]
  write.csv(data.frame(
    dataset_key = dataset_key, computed_methods = paste(compute_methods, collapse = ";"),
    reused_methods = paste(setdiff(full_method_registry$method_id, compute_methods), collapse = ";"),
    reused_file = previous_path, reused_file_sha256 = digest::digest(previous_path, file = TRUE, algo = "sha256"),
    stringsAsFactors = FALSE
  ), file.path(output_dir, "lobo_partial_rerun_manifest.csv"), row.names = FALSE)
}
write.csv(replicates, file.path(output_dir, "leave_one_batch_out_bootstrap_replicates.csv"), row.names = FALSE, na = "")
write.csv(summarise_lobo(replicates), file.path(output_dir, "leave_one_batch_out_bootstrap_summary.csv"), row.names = FALSE, na = "")
cat("Completed derivative-only leave-one-batch-out bootstrap for ", dataset_key, ".\n", sep = "")
