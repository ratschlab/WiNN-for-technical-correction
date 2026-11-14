#!/usr/bin/env Rscript

# Downstream evaluation only. This script cannot execute correction, tuning, or
# parameter-selection code; it reads frozen corrected matrices and saved metrics.

args <- commandArgs(trailingOnly = TRUE)
dataset_arg <- sub("^--dataset=", "", args[grep("^--dataset=", args)])
if (length(dataset_arg) != 1L) stop("Supply exactly one --dataset=<key> argument.", call. = FALSE)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
required_packages <- c("digest", "jsonlite", "matrixStats", "irlba")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
source(file.path(repo_root, "analysis", "derivative_biology_helpers.R"), local = TRUE)

config_path <- file.path(repo_root, "analysis", "config", "derivative_biology_config.json")
config <- jsonlite::read_json(config_path, simplifyVector = TRUE)
config_sha <- digest::digest(config_path, file = TRUE, algo = "sha256")
n_bootstrap <- as.integer(config$common_processing$bootstrap_replicates)
n_permutations <- as.integer(config$common_processing$permutation_replicates)
test_mode <- identical(Sys.getenv("DERIVATIVE_TEST_MODE", "false"), "true")
if (test_mode) {
  n_bootstrap <- as.integer(Sys.getenv("DERIVATIVE_TEST_BOOTSTRAP", "10"))
  n_permutations <- as.integer(Sys.getenv("DERIVATIVE_TEST_PERMUTATIONS", "20"))
}
cores <- max(1L, as.integer(Sys.getenv("DERIVATIVE_CORES", "1")))
if (.Platform$OS.type == "windows") cores <- 1L

dataset_registry <- derivative_dataset_registry()
if (!dataset_arg %in% dataset_registry$dataset_key) stop("Unsupported dataset: ", dataset_arg, call. = FALSE)
dataset_label <- dataset_registry$dataset[match(dataset_arg, dataset_registry$dataset_key)]
method_registry <- derivative_method_registry()
output_dir <- file.path(repo_root, "results", "final", "derivative", dataset_arg)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Partial re-run (e.g. after a WiNN-only package update). DERIVATIVE_COMPUTE_METHODS
# lists the methods whose statistics are recomputed; every other method's rows are
# copied from DERIVATIVE_PREVIOUS_ROOT/<dataset>/ (a previous frozen derivative
# release evaluated on the same frozen matrices, common panel and seeds). Bootstrap
# and permutation draws are generated from fixed seeds before any per-method
# computation, so the spliced rows share every draw with the recomputed ones and
# paired differences remain exactly paired. DERIVATIVE_RESUME_FROM=bootstrap reuses
# the point, association and permutation outputs already written to output_dir by
# an interrupted run of this release and continues with the bootstrap.
full_method_registry <- method_registry
compute_methods <- trimws(strsplit(Sys.getenv("DERIVATIVE_COMPUTE_METHODS", ""), ",", fixed = TRUE)[[1]])
compute_methods <- compute_methods[nzchar(compute_methods)]
partial_mode <- length(compute_methods) > 0L
previous_root <- Sys.getenv("DERIVATIVE_PREVIOUS_ROOT", "")
resume_from_bootstrap <- identical(Sys.getenv("DERIVATIVE_RESUME_FROM", ""), "bootstrap")
if (partial_mode) {
  if (!all(compute_methods %in% method_registry$method_id)) stop("Unknown DERIVATIVE_COMPUTE_METHODS.", call. = FALSE)
  if (!dir.exists(file.path(previous_root, dataset_arg))) stop("DERIVATIVE_PREVIOUS_ROOT lacks ", dataset_arg, call. = FALSE)
  if (!dataset_arg %in% c("waveica", "batchcorr_set1", "clinical_fiams")) {
    stop("Partial re-run is implemented for waveica, batchcorr_set1 and clinical_fiams only.", call. = FALSE)
  }
  method_registry <- method_registry[method_registry$method_id %in% compute_methods, , drop = FALSE]
}
if (resume_from_bootstrap && dataset_arg != "waveica") stop("DERIVATIVE_RESUME_FROM is implemented for waveica only.", call. = FALSE)
reused_files <- character()
splice_previous <- function(value, filename) {
  if (!partial_mode) return(value)
  path <- file.path(previous_root, dataset_arg, filename)
  previous <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  previous <- previous[!previous$method_id %in% compute_methods, , drop = FALSE]
  if (!nrow(previous)) stop(filename, ": previous release has no rows to reuse.", call. = FALSE)
  if (nrow(value)) {
    if (!setequal(names(previous), names(value))) stop(filename, ": columns differ from the previous release.", call. = FALSE)
    value <- value[names(previous)]
  }
  reused_files <<- unique(c(reused_files, path))
  out <- rbind(previous, value)
  out[order(match(out$method_id, full_method_registry$method_id), seq_len(nrow(out))), , drop = FALSE]
}

write_csv <- function(value, filename) {
  write.csv(value, file.path(output_dir, filename), row.names = FALSE, na = "")
}

bind_rows_base <- function(values) {
  values <- Filter(function(value) !is.null(value) && nrow(value), values)
  if (!length(values)) return(data.frame())
  columns <- unique(unlist(lapply(values, names), use.names = FALSE))
  values <- lapply(values, function(value) {
    missing <- setdiff(columns, names(value))
    for (column in missing) value[[column]] <- NA
    value[columns]
  })
  do.call(rbind, values)
}

summarise_pca_bootstrap <- function(replicates, point, directions) {
  summaries <- list(); differences <- list()
  for (metric in names(directions)) {
    values <- replicates[replicates$metric == metric, , drop = FALSE]
    registry <- full_method_registry
    wide <- matrix(NA_real_, nrow = n_bootstrap, ncol = nrow(registry))
    for (method_index in seq_len(nrow(registry))) {
      selected <- values[values$method_id == registry$method_id[method_index], ]
      wide[selected$bootstrap, method_index] <- selected$value
    }
    if (anyNA(wide)) stop(metric, ": bootstrap replicates are incomplete across methods.", call. = FALSE)
    point_values <- point$value[match(
      paste(registry$method_id, metric), paste(point$method_id, point$metric)
    )]
    summaries[[metric]] <- bootstrap_summary(
      wide, point_values, metric, directions[[metric]], registry
    )
    differences[[metric]] <- paired_winn_differences(
      wide, point_values, metric, directions[[metric]], registry
    )
  }
  list(summary = bind_rows_base(summaries), differences = bind_rows_base(differences))
}

run_pca_draws <- function(log_matrices, draws, metric_function) {
  worker <- function(bootstrap_index) {
    draw <- draws[[bootstrap_index]]
    do.call(rbind, lapply(seq_len(nrow(method_registry)), function(method_index) {
      method_id <- method_registry$method_id[method_index]
      values <- metric_function(log_matrices[[method_id]], draw)
      data.frame(
        bootstrap = bootstrap_index, method_id = method_id,
        method = method_registry$method[method_index], metric = names(values),
        value = as.numeric(values), stringsAsFactors = FALSE
      )
    }))
  }
  bind_rows_base(parallel::mclapply(seq_len(n_bootstrap), worker, mc.cores = cores, mc.preschedule = TRUE))
}

permutation_outputs <- function(null_by_method, observed_by_method, endpoint) {
  null <- bind_rows_base(lapply(seq_len(nrow(method_registry)), function(index) {
    method_id <- method_registry$method_id[index]
    data.frame(
      permutation = seq_along(null_by_method[[method_id]]),
      method_id = method_id, method = method_registry$method[index], endpoint = endpoint,
      significant_features = null_by_method[[method_id]], stringsAsFactors = FALSE
    )
  }))
  summary <- bind_rows_base(lapply(seq_len(nrow(method_registry)), function(index) {
    method_id <- method_registry$method_id[index]
    permutation_summary(
      null_by_method[[method_id]], observed_by_method[[method_id]], method_id,
      method_registry$method[index], endpoint
    )
  }))
  list(null = null, summary = summary)
}

pre_manifest <- read.csv(
  file.path(repo_root, "results", "frozen_endpoint_free_inputs", "provenance", "matrix_checksums_pre.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)

if (dataset_arg == "simulation") {
  seed_metrics <- read.csv(
    file.path(repo_root, "results", "final", "tables", "simulation_seed_method_metrics.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  metrics <- c(
    "truth_log1p_rmse", "truth_metabolite_profile_icc_mean", "truth_sample_profile_pearson_mean",
    "residual_run_order_gam_mean_deviance", "residual_ljung_box_proportion_significant",
    "batch_weighted_pc_r2_categorical", "heldout_qc_cv_mean"
  )
  directions <- c(
    truth_log1p_rmse = "lower", truth_metabolite_profile_icc_mean = "higher",
    truth_sample_profile_pearson_mean = "higher", residual_run_order_gam_mean_deviance = "lower",
    residual_ljung_box_proportion_significant = "lower", batch_weighted_pc_r2_categorical = "lower",
    heldout_qc_cv_mean = "lower"
  )
  seed_metrics <- seed_metrics[seed_metrics$metric %in% metrics, , drop = FALSE]
  seed_ids <- sort(unique(seed_metrics$seed_id))
  set.seed(as.integer(config$common_processing$simulation_seed_bootstrap_seed))
  draws <- replicate(n_bootstrap, sample(seed_ids, length(seed_ids), replace = TRUE), simplify = FALSE)
  point <- stats::aggregate(value ~ method_id + method + metric, seed_metrics, mean)
  replicates <- bind_rows_base(lapply(seq_along(draws), function(index) {
    draw <- draws[[index]]
    selected <- do.call(rbind, lapply(seq_along(draw), function(copy_index) {
      value <- seed_metrics[seed_metrics$seed_id == draw[copy_index], , drop = FALSE]
      value$copy_index <- copy_index
      value
    }))
    summary <- stats::aggregate(value ~ method_id + method + metric, selected, mean)
    summary$bootstrap <- index
    summary
  }))
  result <- summarise_pca_bootstrap(replicates, point, directions)
  write_csv(seed_metrics, "simulation_seed_metrics.csv")
  write_csv(replicates, "bootstrap_replicates.csv")
  write_csv(result$summary, "bootstrap_estimates.csv")
  write_csv(result$differences, "bootstrap_paired_differences.csv")
  write_csv(data.frame(
    dataset_key = dataset_arg, completed_seed_count = length(seed_ids),
    completed_seed_ids = paste(seed_ids, collapse = ";"), stringsAsFactors = FALSE
  ), "simulation_seed_inventory.csv")
} else {
  frozen <- derivative_load_matrices(repo_root, dataset_arg, pre_manifest)
  matrices <- frozen$matrices
  meta <- derivative_load_metadata(repo_root, dataset_arg, frozen$sample_ids)
  common_panel <- data.frame(
    dataset_key = dataset_arg, n_features = length(frozen$feature_ids),
    n_samples = length(frozen$sample_ids), feature_panel_sha256 = digest::digest(frozen$feature_ids, algo = "sha256"),
    sample_panel_sha256 = digest::digest(frozen$sample_ids, algo = "sha256"),
    all_method_dimensions_identical = TRUE, stringsAsFactors = FALSE
  )
  write_csv(common_panel, "common_panel.csv")

  if (dataset_arg == "waveica") {
    study <- which(meta$is_study)
    md <- meta[study, , drop = FALSE]
    log_matrices <- lapply(matrices, function(value) log1p(value[, study, drop = FALSE]))
    group <- factor(md$biological_group, levels = c("group_0", "group_1"))
    batch <- factor(md$batch)

    if (resume_from_bootstrap) {
      point <- read.csv(file.path(output_dir, "biological_point_metrics.csv"), stringsAsFactors = FALSE, check.names = FALSE)
      for (filename in c("association_features.csv", "permutation_null.csv", "permutation_summary.csv")) {
        if (!file.exists(file.path(output_dir, filename))) stop("Cannot resume: missing ", filename, call. = FALSE)
      }
      if (!setequal(unique(point$method_id), full_method_registry$method_id)) stop("Cannot resume: incomplete point metrics.", call. = FALSE)
    } else {
    associations <- list(); point_rows <- list()
    for (index in seq_len(nrow(method_registry))) {
      method_id <- method_registry$method_id[index]
      table <- two_group_model(log_matrices[[method_id]], group, "group_0", "group_1")
      table$method_id <- method_id; table$method <- method_registry$method[index]
      associations[[method_id]] <- table
      point_rows[[method_id]] <- data.frame(
        method_id = method_id, method = method_registry$method[index],
        metric = c("batch_weighted_pc_r2_categorical", "biological_weighted_pc_r2", "biological_associated_features", "within_group_variance", "standardized_effect_magnitude"),
        value = c(
          pca_r2_from_log(log_matrices[[method_id]], batch, "categorical"),
          pca_r2_from_log(log_matrices[[method_id]], group, "categorical"),
          sum(table$significant), stats::median(table$within_group_variance, na.rm = TRUE),
          stats::median(abs(table$standardized_effect), na.rm = TRUE)
        ), stringsAsFactors = FALSE
      )
    }
    association_table <- splice_previous(bind_rows_base(associations), "association_features.csv")
    point <- splice_previous(bind_rows_base(point_rows), "biological_point_metrics.csv")
    write_csv(association_table, "association_features.csv")
    write_csv(point, "biological_point_metrics.csv")

    group_permutations <- constrained_group_permutations(
      as.character(group), as.character(batch), n_permutations,
      as.integer(config$common_processing$permutation_seed)
    )
    null <- lapply(method_registry$method_id, function(method_id) {
      two_group_permutation_counts(log_matrices[[method_id]], group_permutations, "group_0", "group_1")
    })
    names(null) <- method_registry$method_id
    observed <- setNames(vapply(associations, function(value) sum(value$significant), numeric(1)), names(associations))
    permutation <- permutation_outputs(null, observed, "biological_group")
    write_csv(splice_previous(permutation$null, "permutation_null.csv"), "permutation_null.csv")
    write_csv(splice_previous(permutation$summary, "permutation_summary.csv"), "permutation_summary.csv")
    }

    units <- md$sample_id
    strata <- interaction(md$batch, group, drop = TRUE)
    draws <- cluster_bootstrap_draws(
      units, strata, n_bootstrap,
      as.integer(config$common_processing$bootstrap_seed)
    )
    index_lookup <- setNames(seq_along(units), units)
    metric_function <- function(log_matrix, draw) {
      index <- unname(index_lookup[draw])
      basis <- pca_basis_from_log(log_matrix[, index, drop = FALSE])
      c(
        batch_weighted_pc_r2_categorical = pca_r2_from_basis(basis, factor(md$batch[index]), "categorical"),
        biological_weighted_pc_r2 = pca_r2_from_basis(basis, factor(group[index], levels = levels(group)), "categorical")
      )
    }
    bootstrap_replicates <- splice_previous(run_pca_draws(log_matrices, draws, metric_function), "bootstrap_replicates.csv")
    pca_point <- point[point$metric %in% c("batch_weighted_pc_r2_categorical", "biological_weighted_pc_r2"), ]
    pca_result <- summarise_pca_bootstrap(
      bootstrap_replicates, pca_point,
      c(batch_weighted_pc_r2_categorical = "lower", biological_weighted_pc_r2 = "higher")
    )
    write_csv(bootstrap_replicates, "bootstrap_replicates.csv")
    write_csv(pca_result$summary, "bootstrap_estimates.csv")
    write_csv(pca_result$differences, "bootstrap_paired_differences.csv")

    lobo <- list()
    for (heldout in levels(batch)) {
      training <- batch != heldout; testing <- batch == heldout
      for (index in seq_len(nrow(method_registry))) {
        method_id <- method_registry$method_id[index]
        train_fit <- two_group_model(log_matrices[[method_id]][, training, drop = FALSE], group[training], "group_0", "group_1")
        test_fit <- two_group_model(log_matrices[[method_id]][, testing, drop = FALSE], group[testing], "group_0", "group_1")
        if (is.null(train_fit) || is.null(test_fit)) next
        primary <- effect_concordance(train_fit$effect, test_fit$effect)
        discovered <- effect_concordance(train_fit$effect, test_fit$effect, train_fit$significant)
        lobo[[length(lobo) + 1L]] <- data.frame(
          heldout_batch = heldout, method_id = method_id, method = method_registry$method[index],
          feature_panel = c("all_common_features", "training_fdr_features"),
          pearson = c(primary["pearson"], discovered["pearson"]),
          spearman = c(primary["spearman"], discovered["spearman"]),
          sign_concordance = c(primary["sign_concordance"], discovered["sign_concordance"]),
          attenuation = c(primary["attenuation"], discovered["attenuation"]),
          n_features = c(primary["n_features"], discovered["n_features"]),
          estimable = TRUE, stringsAsFactors = FALSE
        )
      }
    }
    write_csv(splice_previous(bind_rows_base(lobo), "leave_one_batch_out_effects.csv"), "leave_one_batch_out_effects.csv")
  } else if (dataset_arg == "mtbls79") {
    study <- which(meta$is_study)
    md <- meta[study, , drop = FALSE]
    log_injection <- lapply(matrices, function(value) log1p(value[, study, drop = FALSE]))
    units <- unique(md$sample_rep)
    unit_class <- vapply(units, function(unit) unique(md$class[md$sample_rep == unit])[1], character(1))
    collapsed <- lapply(log_injection, function(value) {
      result <- vapply(units, function(unit) rowMeans(value[, md$sample_rep == unit, drop = FALSE]), numeric(nrow(value)))
      rownames(result) <- rownames(value); colnames(result) <- units; result
    })
    associations <- list(); point_rows <- list()
    for (index in seq_len(nrow(method_registry))) {
      method_id <- method_registry$method_id[index]
      table <- two_group_model(collapsed[[method_id]], unit_class, "S", "C")
      table$method_id <- method_id; table$method <- method_registry$method[index]
      associations[[method_id]] <- table
      point_rows[[method_id]] <- data.frame(
        method_id = method_id, method = method_registry$method[index],
        metric = c("batch_weighted_pc_r2_categorical", "biological_weighted_pc_r2", "biological_associated_features", "within_group_variance", "standardized_effect_magnitude"),
        value = c(
          pca_r2_from_log(log_injection[[method_id]], factor(md$batch), "categorical"),
          pca_r2_from_log(collapsed[[method_id]], factor(unit_class), "categorical"),
          sum(table$significant), stats::median(table$within_group_variance, na.rm = TRUE),
          stats::median(abs(table$standardized_effect), na.rm = TRUE)
        ), stringsAsFactors = FALSE
      )
    }
    point <- bind_rows_base(point_rows)
    write_csv(bind_rows_base(associations), "association_features.csv")
    write_csv(point, "biological_point_metrics.csv")

    group_permutations <- constrained_group_permutations(unit_class, rep("all", length(units)), n_permutations, as.integer(config$common_processing$permutation_seed) + 79L)
    null <- lapply(method_registry$method_id, function(method_id) {
      two_group_permutation_counts(collapsed[[method_id]], group_permutations, "S", "C")
    })
    names(null) <- method_registry$method_id
    observed <- setNames(vapply(associations, function(value) sum(value$significant), numeric(1)), names(associations))
    permutation <- permutation_outputs(null, observed, "cow_vs_sheep")
    write_csv(permutation$null, "permutation_null.csv")
    write_csv(permutation$summary, "permutation_summary.csv")

    draws <- cluster_bootstrap_draws(units, unit_class, n_bootstrap, as.integer(config$common_processing$bootstrap_seed) + 79L)
    unit_lookup <- setNames(seq_along(units), units)
    injection_lookup <- split(seq_len(nrow(md)), md$sample_rep)
    metric_function <- function(log_matrix, draw) {
      collapsed_index <- unname(unit_lookup[draw])
      injection_index <- unlist(injection_lookup[draw], use.names = FALSE)
      batch_basis <- pca_basis_from_log(log_matrix[, injection_index, drop = FALSE])
      biological_matrix <- vapply(draw, function(unit) {
        rowMeans(log_matrix[, injection_lookup[[unit]], drop = FALSE])
      }, numeric(nrow(log_matrix)))
      rownames(biological_matrix) <- rownames(log_matrix)
      biological_basis <- pca_basis_from_log(biological_matrix)
      c(
        batch_weighted_pc_r2_categorical = pca_r2_from_basis(batch_basis, factor(md$batch[injection_index]), "categorical"),
        biological_weighted_pc_r2 = pca_r2_from_basis(biological_basis, factor(unit_class[collapsed_index]), "categorical")
      )
    }
    bootstrap_replicates <- run_pca_draws(log_injection, draws, metric_function)
    pca_point <- point[point$metric %in% c("batch_weighted_pc_r2_categorical", "biological_weighted_pc_r2"), ]
    pca_result <- summarise_pca_bootstrap(
      bootstrap_replicates, pca_point,
      c(batch_weighted_pc_r2_categorical = "lower", biological_weighted_pc_r2 = "higher")
    )
    write_csv(bootstrap_replicates, "bootstrap_replicates.csv")
    write_csv(pca_result$summary, "bootstrap_estimates.csv")
    write_csv(pca_result$differences, "bootstrap_paired_differences.csv")

    lobo <- list()
    for (heldout in unique(md$batch)) {
      for (index in seq_len(nrow(method_registry))) {
        method_id <- method_registry$method_id[index]
        collapse_fold <- function(selected) {
          selected_units <- unique(md$sample_rep[selected])
          matrix <- vapply(selected_units, function(unit) {
            columns <- selected & md$sample_rep == unit
            rowMeans(log_injection[[method_id]][, columns, drop = FALSE])
          }, numeric(nrow(log_injection[[method_id]])))
          list(matrix = matrix, class = vapply(selected_units, function(unit) unique(md$class[md$sample_rep == unit])[1], character(1)))
        }
        train <- collapse_fold(md$batch != heldout); test <- collapse_fold(md$batch == heldout)
        train_fit <- two_group_model(train$matrix, train$class, "S", "C")
        test_fit <- two_group_model(test$matrix, test$class, "S", "C")
        if (is.null(train_fit) || is.null(test_fit)) next
        primary <- effect_concordance(train_fit$effect, test_fit$effect)
        discovered <- effect_concordance(train_fit$effect, test_fit$effect, train_fit$significant)
        lobo[[length(lobo) + 1L]] <- data.frame(
          heldout_batch = heldout, method_id = method_id, method = method_registry$method[index],
          feature_panel = c("all_common_features", "training_fdr_features"),
          pearson = c(primary["pearson"], discovered["pearson"]),
          spearman = c(primary["spearman"], discovered["spearman"]),
          sign_concordance = c(primary["sign_concordance"], discovered["sign_concordance"]),
          attenuation = c(primary["attenuation"], discovered["attenuation"]),
          n_features = c(primary["n_features"], discovered["n_features"]), estimable = TRUE,
          stringsAsFactors = FALSE
        )
      }
    }
    write_csv(bind_rows_base(lobo), "leave_one_batch_out_effects.csv")
  } else if (dataset_arg == "sacurine") {
    study <- which(meta$is_study)
    md <- meta[study, , drop = FALSE]
    md$age_z <- as.numeric(scale(md$age)); md$bmi_z <- as.numeric(scale(md$bmi))
    md$gender <- factor(md$gender); md$sampling_factor <- factor(md$sampling)
    raw_log <- lapply(matrices, function(value) log1p(value[, study, drop = FALSE]))
    biological_log <- lapply(matrices, function(value) log1p(sweep(value[, study, drop = FALSE], 2L, md$osmolality, "/")))
    full_design <- stats::model.matrix(~ age_z + bmi_z + gender + sampling_factor, data = md)
    coefficient_names <- c(age = "age_z", bmi = "bmi_z", gender = grep("^gender", colnames(full_design), value = TRUE)[1])
    point_rows <- list(); associations <- list()
    for (method_index in seq_len(nrow(method_registry))) {
      method_id <- method_registry$method_id[method_index]
      raw_basis <- pca_basis_from_log(raw_log[[method_id]])
      biological_basis <- pca_basis_from_log(biological_log[[method_id]])
      point_rows[[method_id]] <- data.frame(
        method_id = method_id, method = method_registry$method[method_index],
        metric = c("batch_weighted_pc_r2_categorical", "age_weighted_pc_r2", "bmi_weighted_pc_r2", "gender_weighted_pc_r2"),
        value = c(
          pca_r2_from_basis(raw_basis, factor(md$batch), "categorical"),
          pca_r2_from_basis(biological_basis, md$age, "continuous"),
          pca_r2_from_basis(biological_basis, md$bmi, "continuous"),
          pca_r2_from_basis(biological_basis, md$gender, "categorical")
        ), stringsAsFactors = FALSE
      )
      associations[[method_id]] <- bind_rows_base(lapply(names(coefficient_names), function(endpoint) {
        fit <- ols_feature_model(biological_log[[method_id]], full_design, match(coefficient_names[[endpoint]], colnames(full_design)))
        table <- fit$table; table$endpoint <- endpoint; table$method_id <- method_id; table$method <- method_registry$method[method_index]; table
      }))
    }
    point <- bind_rows_base(point_rows)
    association_table <- bind_rows_base(associations)
    write_csv(point, "biological_point_metrics.csv")
    write_csv(association_table, "sacurine_covariate_associations.csv")

    permutations <- restricted_index_permutations(md$batch, n_permutations, as.integer(config$common_processing$permutation_seed) + 113L)
    all_null <- list(); all_summary <- list()
    for (endpoint in names(coefficient_names)) {
      reduced_formula <- switch(endpoint,
        age = ~ bmi_z + gender + sampling_factor,
        bmi = ~ age_z + gender + sampling_factor,
        gender = ~ age_z + bmi_z + sampling_factor
      )
      reduced_design <- stats::model.matrix(reduced_formula, data = md)
      coefficient_index <- match(coefficient_names[[endpoint]], colnames(full_design))
      null <- lapply(method_registry$method_id, function(method_id) {
        freedman_lane_counts(biological_log[[method_id]], full_design, reduced_design, coefficient_index, permutations)
      })
      names(null) <- method_registry$method_id
      observed <- setNames(vapply(method_registry$method_id, function(method_id) {
        sum(associations[[method_id]]$significant[associations[[method_id]]$endpoint == endpoint])
      }, numeric(1)), method_registry$method_id)
      output <- permutation_outputs(null, observed, endpoint)
      all_null[[endpoint]] <- output$null; all_summary[[endpoint]] <- output$summary
    }
    write_csv(bind_rows_base(all_null), "permutation_null.csv")
    write_csv(bind_rows_base(all_summary), "permutation_summary.csv")

    draws <- cluster_bootstrap_draws(md$sample_id, md$batch, n_bootstrap, as.integer(config$common_processing$bootstrap_seed) + 113L)
    lookup <- setNames(seq_len(nrow(md)), md$sample_id)
    combined_log <- lapply(method_registry$method_id, function(method_id) list(raw = raw_log[[method_id]], biological = biological_log[[method_id]]))
    names(combined_log) <- method_registry$method_id
    metric_function <- function(value, draw) {
      index <- unname(lookup[draw])
      raw_basis <- pca_basis_from_log(value$raw[, index, drop = FALSE])
      bio_basis <- pca_basis_from_log(value$biological[, index, drop = FALSE])
      c(
        batch_weighted_pc_r2_categorical = pca_r2_from_basis(raw_basis, factor(md$batch[index]), "categorical"),
        age_weighted_pc_r2 = pca_r2_from_basis(bio_basis, md$age[index], "continuous"),
        bmi_weighted_pc_r2 = pca_r2_from_basis(bio_basis, md$bmi[index], "continuous"),
        gender_weighted_pc_r2 = pca_r2_from_basis(bio_basis, factor(md$gender[index]), "categorical")
      )
    }
    bootstrap_replicates <- run_pca_draws(combined_log, draws, metric_function)
    directions <- c(
      batch_weighted_pc_r2_categorical = "lower", age_weighted_pc_r2 = "higher",
      bmi_weighted_pc_r2 = "higher", gender_weighted_pc_r2 = "higher"
    )
    pca_result <- summarise_pca_bootstrap(bootstrap_replicates, point, directions)
    write_csv(bootstrap_replicates, "bootstrap_replicates.csv")
    write_csv(pca_result$summary, "bootstrap_estimates.csv")
    write_csv(pca_result$differences, "bootstrap_paired_differences.csv")

    lobo <- list()
    for (heldout in unique(md$batch)) {
      for (method_index in seq_len(nrow(method_registry))) {
        method_id <- method_registry$method_id[method_index]
        for (endpoint in names(coefficient_names)) {
          train_md <- droplevels(md[md$batch != heldout, , drop = FALSE])
          test_md <- droplevels(md[md$batch == heldout, , drop = FALSE])
          train_design <- stats::model.matrix(~ age_z + bmi_z + gender + sampling_factor, data = train_md)
          test_design <- stats::model.matrix(~ age_z + bmi_z + gender + sampling_factor, data = test_md)
          coefficient <- coefficient_names[[endpoint]]
          if (!coefficient %in% colnames(train_design) || !coefficient %in% colnames(test_design)) next
          train_fit <- ols_feature_model(biological_log[[method_id]][, md$batch != heldout, drop = FALSE], train_design, match(coefficient, colnames(train_design)))
          test_fit <- ols_feature_model(biological_log[[method_id]][, md$batch == heldout, drop = FALSE], test_design, match(coefficient, colnames(test_design)))
          if (is.null(train_fit) || is.null(test_fit)) next
          primary <- effect_concordance(train_fit$table$coefficient, test_fit$table$coefficient)
          discovered <- effect_concordance(train_fit$table$coefficient, test_fit$table$coefficient, train_fit$table$significant)
          lobo[[length(lobo) + 1L]] <- data.frame(
            heldout_batch = heldout, endpoint = endpoint, method_id = method_id, method = method_registry$method[method_index],
            feature_panel = c("all_common_features", "training_fdr_features"),
            pearson = c(primary["pearson"], discovered["pearson"]),
            spearman = c(primary["spearman"], discovered["spearman"]),
            sign_concordance = c(primary["sign_concordance"], discovered["sign_concordance"]),
            attenuation = c(primary["attenuation"], discovered["attenuation"]),
            n_features = c(primary["n_features"], discovered["n_features"]), estimable = TRUE,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    write_csv(bind_rows_base(lobo), "leave_one_batch_out_effects.csv")
  } else if (dataset_arg %in% c("batchcorr_set1", "clinical_fiams")) {
    study <- which(meta$is_study)
    md <- meta[study, , drop = FALSE]
    identity_field <- if (dataset_arg == "batchcorr_set1") "accession_id" else "sample"
    identities <- as.character(md[[identity_field]])
    log_matrices <- lapply(matrices, function(value) log1p(value[, study, drop = FALSE]))
    identity_metrics <- list(); matching_metrics <- list(); point_rows <- list()
    for (method_index in seq_len(nrow(method_registry))) {
      method_id <- method_registry$method_id[method_index]
      replicate_table <- replicate_identity_metrics(log_matrices[[method_id]], md, identity_field)
      match_table <- cross_batch_matching(log_matrices[[method_id]], md, identity_field)
      repeatability <- feature_repeatability_ratio(log_matrices[[method_id]], identities)
      replicate_table$method_id <- method_id; replicate_table$method <- method_registry$method[method_index]
      match_table$identity$method_id <- method_id; match_table$identity$method <- method_registry$method[method_index]
      if (dataset_arg == "clinical_fiams") {
        replicate_table$identity <- hash_private_identity(replicate_table$identity)
        match_table$identity$identity <- hash_private_identity(match_table$identity$identity)
      }
      identity_metrics[[method_id]] <- replicate_table
      matching_metrics[[method_id]] <- match_table$identity
      basis <- pca_basis_from_log(log_matrices[[method_id]])
      point_rows[[method_id]] <- data.frame(
        method_id = method_id, method = method_registry$method[method_index],
        metric = c(
          "batch_weighted_pc_r2_categorical", "sample_replicate_pearson_median",
          "genuine_replicate_icc_median", "feature_repeatability_ratio_median",
          "cross_batch_top1_accuracy", "cross_batch_mean_reciprocal_rank"
        ),
        value = c(
          pca_r2_from_basis(basis, factor(md$batch), "categorical"),
          stats::median(replicate_table$pearson, na.rm = TRUE),
          stats::median(replicate_table$icc_a1, na.rm = TRUE),
          stats::median(repeatability, na.rm = TRUE),
          mean(match_table$query$top1, na.rm = TRUE),
          mean(match_table$query$reciprocal_rank, na.rm = TRUE)
        ), stringsAsFactors = FALSE
      )
    }
    point <- splice_previous(bind_rows_base(point_rows), "biological_point_metrics.csv")
    write_csv(point, "biological_point_metrics.csv")
    write_csv(splice_previous(bind_rows_base(identity_metrics), "replicate_metrics_by_identity.csv"), "replicate_metrics_by_identity.csv")
    write_csv(splice_previous(bind_rows_base(matching_metrics), "matching_metrics_by_identity.csv"), "matching_metrics_by_identity.csv")

    units <- unique(identities)
    draws <- cluster_bootstrap_draws(units, NULL, n_bootstrap, as.integer(config$common_processing$bootstrap_seed) + ifelse(dataset_arg == "batchcorr_set1", 199L, 1594L))
    unit_index <- split(seq_along(identities), identities)
    private_id_map <- if (dataset_arg == "clinical_fiams") {
      stats::setNames(hash_private_identity(units), units)
    } else NULL
    bootstrap_worker <- function(bootstrap_index) {
      draw <- draws[[bootstrap_index]]
      bind_rows_base(lapply(seq_len(nrow(method_registry)), function(method_index) {
        method_id <- method_registry$method_id[method_index]
        rep_table <- identity_metrics[[method_id]]
        match_table <- matching_metrics[[method_id]]
        rep_identity <- if (dataset_arg == "clinical_fiams") unname(private_id_map[draw]) else draw
        rep_selected <- rep_table[match(rep_identity, rep_table$identity), , drop = FALSE]
        match_selected <- match_table[match(rep_identity, match_table$identity), , drop = FALSE]
        column_index <- unlist(unit_index[draw], use.names = FALSE)
        batch_basis <- pca_basis_from_log(log_matrices[[method_id]][, column_index, drop = FALSE])
        values <- c(
          batch_weighted_pc_r2_categorical = pca_r2_from_basis(batch_basis, factor(md$batch[column_index]), "categorical"),
          sample_replicate_pearson_median = stats::median(rep_selected$pearson, na.rm = TRUE),
          genuine_replicate_icc_median = stats::median(rep_selected$icc_a1, na.rm = TRUE),
          cross_batch_top1_accuracy = mean(match_selected$top1, na.rm = TRUE),
          cross_batch_mean_reciprocal_rank = mean(match_selected$reciprocal_rank, na.rm = TRUE)
        )
        data.frame(
          bootstrap = bootstrap_index, method_id = method_id, method = method_registry$method[method_index],
          metric = names(values), value = as.numeric(values), stringsAsFactors = FALSE
        )
      }))
    }
    bootstrap_replicates <- splice_previous(bind_rows_base(parallel::mclapply(
      seq_along(draws), bootstrap_worker, mc.cores = cores, mc.preschedule = TRUE
    )), "bootstrap_replicates.csv")
    directions <- c(
      batch_weighted_pc_r2_categorical = "lower", sample_replicate_pearson_median = "higher",
      genuine_replicate_icc_median = "higher", cross_batch_top1_accuracy = "higher",
      cross_batch_mean_reciprocal_rank = "higher"
    )
    point_boot <- point[point$metric %in% names(directions), ]
    bootstrap_result <- summarise_pca_bootstrap(bootstrap_replicates, point_boot, directions)
    write_csv(bootstrap_replicates, "bootstrap_replicates.csv")
    write_csv(bootstrap_result$summary, "bootstrap_estimates.csv")
    write_csv(bootstrap_result$differences, "bootstrap_paired_differences.csv")
  }
}

manifest <- data.frame(
  dataset_key = dataset_arg, dataset = dataset_label,
  config_sha256 = config_sha, test_mode = test_mode,
  bootstrap_replicates = n_bootstrap,
  permutation_replicates = if (dataset_arg %in% c("waveica", "mtbls79", "sacurine")) n_permutations else 0L,
  correction_methods_executed = FALSE, tuning_executed = FALSE, selection_executed = FALSE,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)
write_csv(manifest, "derivative_run_manifest.csv")
if (partial_mode || resume_from_bootstrap) {
  write_csv(data.frame(
    dataset_key = dataset_arg,
    computed_methods = paste(if (partial_mode) compute_methods else full_method_registry$method_id, collapse = ";"),
    reused_methods = paste(setdiff(full_method_registry$method_id, if (partial_mode) compute_methods else full_method_registry$method_id), collapse = ";"),
    resumed_from_bootstrap = resume_from_bootstrap,
    reused_file = if (length(reused_files)) reused_files else NA_character_,
    reused_file_sha256 = if (length(reused_files)) vapply(reused_files, digest::digest, character(1), file = TRUE, algo = "sha256") else NA_character_,
    stringsAsFactors = FALSE
  ), "derivative_partial_rerun_manifest.csv")
}
utils::capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
cat("Completed derivative-only evaluation for ", dataset_arg, ".\n", sep = "")
