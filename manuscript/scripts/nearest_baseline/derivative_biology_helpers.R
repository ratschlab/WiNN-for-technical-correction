# Downstream evaluation helpers for the frozen endpoint-free WiNN release.
# This file deliberately contains no correction, tuning, or selection calls.

derivative_method_registry <- function() {
  data.frame(
    method_id = c(
      "raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf",
      "winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc"
    ),
    method = c(
      "Raw", "ComBat", "QC-RLSC", "QC-RFSC", "TIGER", "SERRF",
      "WiNN auto (QC)", "WiNN auto-batch (QC)", "WiNN fixed default (no QC)"
    ),
    method_order = 1:9,
    stringsAsFactors = FALSE
  )
}

derivative_dataset_registry <- function() {
  data.frame(
    dataset_key = c("simulation", "mtbls79", "clinical_fiams", "batchcorr_set1", "sacurine", "waveica"),
    dataset = c("Simulation", "MTBLS79", "Clinical FIA-MS", "BatchCorr Set 1", "Sacurine", "WaveICA"),
    stringsAsFactors = FALSE
  )
}

as_logical_strict <- function(value, field) {
  if (is.logical(value)) return(value)
  result <- ifelse(tolower(as.character(value)) == "true", TRUE,
    ifelse(tolower(as.character(value)) == "false", FALSE, NA))
  if (anyNA(result)) stop("Could not parse logical field: ", field, call. = FALSE)
  result
}

derivative_load_metadata <- function(repo_root, dataset_key, sample_ids) {
  if (dataset_key == "simulation") {
    meta <- read.csv(
      file.path(repo_root, "data", "simulated", "canonical", "SIM01", "sample_metadata.csv"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    meta$is_qc <- as_logical_strict(meta$is_qc, "is_qc")
    meta$is_study <- !meta$is_qc
  } else if (dataset_key == "mtbls79") {
    meta <- read.csv(
      file.path(repo_root, "data", "public", "processed", "MTBLS79_metadata.csv"),
      row.names = 1, stringsAsFactors = FALSE, check.names = FALSE
    )
    meta$sample_id <- rownames(meta)
    meta$batch <- as.character(meta$Batch)
    meta$run_order <- seq_len(nrow(meta))
    meta$class <- as.character(meta$Class)
    meta$sample_rep <- as.character(meta$Sample_Rep)
    meta$is_qc <- meta$class == "QC"
    meta$is_study <- !meta$is_qc
  } else if (dataset_key == "clinical_fiams") {
    meta <- read.csv(
      file.path(repo_root, "data", "private", "clinical_fiams", "pair_metadata.csv"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    meta$batch <- as.character(meta$batch)
    meta$is_qc <- meta$role != "clinical"
    meta$is_study <- meta$role == "clinical"
  } else if (dataset_key == "batchcorr_set1") {
    meta <- read.csv(
      file.path(repo_root, "data", "public", "batchcorr_set1", "processed", "BatchCorr_set1_metadata.csv"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    meta$batch <- as.character(meta$batch)
    meta$is_qc <- as_logical_strict(meta$is_qc, "is_qc")
    meta$is_study <- as_logical_strict(meta$is_study_sample, "is_study_sample")
    meta$accession_id <- as.character(meta$accession_id)
  } else if (dataset_key == "sacurine") {
    meta <- read.csv(
      file.path(repo_root, "data", "public", "sacurine", "processed", "sacurine_metadata.csv"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    meta$batch <- as.character(meta$batch)
    meta$is_qc <- as_logical_strict(meta$is_qc, "is_qc")
    meta$is_study <- !meta$is_qc
  } else if (dataset_key == "waveica") {
    meta <- read.csv(
      file.path(repo_root, "data", "public", "waveica_adenocarcinoma", "processed", "waveica_metadata.csv"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    meta$batch <- as.character(meta$batch)
    meta$is_qc <- as_logical_strict(meta$is_qc, "is_qc")
    meta$is_study <- !meta$is_qc
    meta$biological_group <- as.character(meta$biological_group)
  } else {
    stop("Unsupported derivative dataset: ", dataset_key, call. = FALSE)
  }
  meta$sample_id <- as.character(meta$sample_id)
  meta <- meta[match(sample_ids, meta$sample_id), , drop = FALSE]
  if (anyNA(meta$sample_id) || !identical(as.character(meta$sample_id), sample_ids)) {
    stop(dataset_key, ": frozen metadata do not align with frozen matrices.", call. = FALSE)
  }
  if (anyNA(meta$batch) || anyNA(meta$is_study)) {
    stop(dataset_key, ": missing batch or study-sample status.", call. = FALSE)
  }
  meta
}

derivative_load_matrices <- function(repo_root, dataset_key, pre_manifest = NULL) {
  registry <- derivative_method_registry()
  matrix_root <- file.path(repo_root, "results", "frozen_endpoint_free_inputs", "primary", dataset_key)
  paths <- file.path(matrix_root, registry$method_id, "corrected_matrix.rds")
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing frozen matrix: ", paste(missing, collapse = ", "), call. = FALSE)

  if (!is.null(pre_manifest)) {
    expected_path <- file.path("results", "primary", dataset_key, registry$method_id, "corrected_matrix.rds")
    expected_sha <- pre_manifest$sha256[match(expected_path, pre_manifest$path)]
    if (anyNA(expected_sha)) stop(dataset_key, ": primary matrices are absent from the frozen pre-analysis manifest.", call. = FALSE)
    observed_sha <- vapply(paths, digest::digest, character(1), file = TRUE, algo = "sha256")
    if (!identical(unname(observed_sha), unname(expected_sha))) {
      stop(dataset_key, ": a local frozen matrix does not match the endpoint-free release checksum.", call. = FALSE)
    }
  }

  matrices <- lapply(paths, readRDS)
  names(matrices) <- registry$method_id
  matrices <- lapply(matrices, function(value) {
    value <- as.matrix(value)
    storage.mode(value) <- "double"
    value
  })
  raw <- matrices[["raw"]]
  common_features <- rownames(raw)
  common_samples <- colnames(raw)
  for (method_id in names(matrices)) {
    value <- matrices[[method_id]]
    if (any(!is.finite(value)) || any(value < 0)) stop(dataset_key, "/", method_id, ": non-finite or negative matrix value.", call. = FALSE)
    common_features <- common_features[common_features %in% rownames(value)]
    common_samples <- common_samples[common_samples %in% colnames(value)]
  }
  if (!length(common_features) || !length(common_samples)) stop(dataset_key, ": empty common matrix panel.", call. = FALSE)
  matrices <- lapply(matrices, function(value) value[common_features, common_samples, drop = FALSE])
  for (method_id in names(matrices)) {
    if (!identical(dimnames(matrices[[method_id]]), dimnames(matrices[["raw"]]))) {
      stop(dataset_key, ": method matrices do not share an identical common panel.", call. = FALSE)
    }
  }
  list(
    matrices = matrices,
    methods = registry,
    feature_ids = common_features,
    sample_ids = common_samples,
    matrix_paths = paths
  )
}

pca_basis_from_log_exact <- function(log_matrix) {
  log_matrix <- as.matrix(log_matrix)
  if (ncol(log_matrix) < 3L) return(NULL)
  z <- t(log_matrix)
  means <- colMeans(z)
  standard_deviation <- apply(z, 2L, stats::sd)
  keep <- is.finite(standard_deviation) & standard_deviation > 0 & is.finite(means)
  z <- sweep(z[, keep, drop = FALSE], 2L, means[keep], "-")
  z <- sweep(z, 2L, standard_deviation[keep], "/")
  if (nrow(z) < 3L || ncol(z) < 2L) return(NULL)

  if (nrow(z) <= ncol(z)) {
    decomposition <- eigen(tcrossprod(z), symmetric = TRUE)
    values <- pmax(decomposition$values, 0)
    keep_pc <- seq_len(min(5L, sum(values > .Machine$double.eps), ncol(decomposition$vectors)))
    scores <- decomposition$vectors[, keep_pc, drop = FALSE]
  } else {
    decomposition <- eigen(crossprod(z), symmetric = TRUE)
    values <- pmax(decomposition$values, 0)
    keep_pc <- seq_len(min(5L, sum(values > .Machine$double.eps), ncol(decomposition$vectors)))
    scores <- z %*% decomposition$vectors[, keep_pc, drop = FALSE]
  }
  if (!length(keep_pc)) return(NULL)
  weights <- values[keep_pc]
  list(scores = scores, weights = weights)
}

pca_basis_from_log <- function(log_matrix) {
  log_matrix <- as.matrix(log_matrix)
  if (ncol(log_matrix) < 3L) return(NULL)
  z <- t(log_matrix)
  means <- colMeans(z)
  standard_deviation <- apply(z, 2L, stats::sd)
  keep <- is.finite(standard_deviation) & standard_deviation > 0 & is.finite(means)
  z <- sweep(z[, keep, drop = FALSE], 2L, means[keep], "-")
  z <- sweep(z, 2L, standard_deviation[keep], "/")
  if (nrow(z) < 3L || ncol(z) < 2L) return(NULL)

  n_components <- min(5L, min(dim(z)) - 1L)
  if (n_components < 1L || min(dim(z)) <= 300L || !requireNamespace("irlba", quietly = TRUE)) {
    return(pca_basis_from_log_exact(log_matrix))
  }

  # A fixed starting vector makes the truncated decomposition deterministic.
  # The exact eigendecomposition remains the fallback for small or exceptional
  # matrices. The retained weighted-PC R-squared is numerically equivalent to
  # the exact top-five decomposition but avoids repeatedly computing unused PCs
  # in the paired bootstrap.
  decomposition <- tryCatch(
    irlba::irlba(
      z, nv = n_components, nu = n_components, tol = 1e-10, maxit = 2000L,
      v = rep(1 / sqrt(ncol(z)), ncol(z))
    ),
    error = function(condition) NULL
  )
  if (is.null(decomposition) || any(!is.finite(decomposition$d))) {
    return(pca_basis_from_log_exact(log_matrix))
  }
  list(scores = decomposition$u, weights = decomposition$d^2)
}

pca_r2_from_basis <- function(basis, target, target_type) {
  target_type <- match.arg(target_type, c("categorical", "continuous"))
  if (is.null(basis) || !is.list(basis)) return(NA_real_)
  if (length(target) != nrow(basis$scores)) stop("PCA target length mismatch.", call. = FALSE)
  observed <- if (target_type == "continuous") is.finite(target) else !is.na(target)
  if (sum(observed) < 3L || !all(observed)) {
    stop("PCA targets must be complete; subset the matrix and target together before decomposition.", call. = FALSE)
  }
  scores <- basis$scores
  weights <- basis$weights
  r_squared <- vapply(seq_along(weights), function(index) {
    score <- scores[, index]
    if (target_type == "continuous") {
      if (!is.numeric(target) || stats::sd(target) <= 0) return(NA_real_)
      stats::cor(score, target)^2
    } else {
      group <- factor(target)
      if (nlevels(group) < 2L) return(NA_real_)
      grand <- mean(score)
      total <- sum((score - grand)^2)
      if (total <= 0) return(NA_real_)
      indices <- split(seq_along(group), group)
      between <- sum(vapply(indices, function(i) length(i) * (mean(score[i]) - grand)^2, numeric(1)))
      between / total
    }
  }, numeric(1))
  valid <- is.finite(r_squared) & is.finite(weights) & weights > 0
  if (!any(valid)) return(NA_real_)
  sum(weights[valid] * r_squared[valid]) / sum(weights[valid])
}

pca_r2_from_log <- function(log_matrix, target, target_type) {
  target_type <- match.arg(target_type, c("categorical", "continuous"))
  observed <- if (target_type == "continuous") is.finite(target) else !is.na(target)
  if (length(target) != ncol(log_matrix)) stop("PCA target length mismatch.", call. = FALSE)
  pca_r2_from_basis(
    pca_basis_from_log(log_matrix[, observed, drop = FALSE]),
    target[observed], target_type
  )
}

weighted_pc_r2_frozen <- function(matrix, target, target_type) {
  if (target_type == "categorical" && is.numeric(target)) {
    stop("Categorical PCA R-squared target must be supplied as a factor or character vector.", call. = FALSE)
  }
  pca_r2_from_log(log1p(matrix), target, target_type)
}

two_group_model <- function(log_matrix, group, level0, level1, fdr = 0.05) {
  group <- factor(as.character(group), levels = c(level0, level1))
  keep <- !is.na(group)
  log_matrix <- log_matrix[, keep, drop = FALSE]
  group <- droplevels(group[keep])
  i0 <- which(group == level0); i1 <- which(group == level1)
  if (length(i0) < 2L || length(i1) < 2L) return(NULL)
  mean0 <- rowMeans(log_matrix[, i0, drop = FALSE])
  mean1 <- rowMeans(log_matrix[, i1, drop = FALSE])
  variance0 <- matrixStats::rowVars(log_matrix[, i0, drop = FALSE])
  variance1 <- matrixStats::rowVars(log_matrix[, i1, drop = FALSE])
  pooled_variance <- ((length(i0) - 1L) * variance0 + (length(i1) - 1L) * variance1) /
    (length(i0) + length(i1) - 2L)
  standard_error <- sqrt(pooled_variance * (1 / length(i0) + 1 / length(i1)))
  effect <- mean1 - mean0
  t_statistic <- effect / standard_error
  p_value <- 2 * stats::pt(-abs(t_statistic), df = length(i0) + length(i1) - 2L)
  p_value[!is.finite(p_value)] <- NA_real_
  p_adjusted <- stats::p.adjust(p_value, method = "BH")
  pooled_sd <- sqrt(pooled_variance)
  data.frame(
    feature_id = rownames(log_matrix), effect = effect,
    standardized_effect = effect / pooled_sd,
    within_group_variance = pooled_variance,
    t_statistic = t_statistic, p_value = p_value, p_adj = p_adjusted,
    significant = is.finite(p_adjusted) & p_adjusted < fdr,
    stringsAsFactors = FALSE
  )
}

constrained_group_permutations <- function(group, blocks, n_permutations, seed) {
  group <- as.character(group)
  blocks <- as.character(blocks)
  set.seed(seed)
  permutations <- replicate(n_permutations, {
    result <- group
    for (block in unique(blocks)) {
      index <- which(blocks == block)
      result[index] <- sample(group[index], length(index), replace = FALSE)
    }
    result
  }, simplify = "matrix")
  permutations
}

two_group_permutation_counts <- function(log_matrix, group_permutations, level0, level1, fdr = 0.05, chunk_size = 250L) {
  n <- ncol(log_matrix)
  if (nrow(group_permutations) != n) stop("Permutation dimension mismatch.", call. = FALSE)
  sums_squared <- rowSums(log_matrix^2)
  result <- integer(ncol(group_permutations))
  chunks <- split(seq_len(ncol(group_permutations)), ceiling(seq_len(ncol(group_permutations)) / chunk_size))
  for (indices in chunks) {
    labels <- group_permutations[, indices, drop = FALSE] == level1
    n1 <- colSums(labels); n0 <- n - n1
    weights1 <- sweep(labels, 2L, n1, "/")
    weights0 <- sweep(!labels, 2L, n0, "/")
    mean1 <- log_matrix %*% weights1
    mean0 <- log_matrix %*% weights0
    ss_within <- sums_squared - sweep(mean1^2, 2L, n1, "*") - sweep(mean0^2, 2L, n0, "*")
    variance <- sweep(ss_within, 2L, n - 2L, "/")
    standard_error <- sqrt(sweep(variance, 2L, 1 / n1 + 1 / n0, "*"))
    statistic <- (mean1 - mean0) / standard_error
    p_values <- 2 * stats::pt(-abs(statistic), df = n - 2L)
    for (column in seq_along(indices)) {
      adjusted <- stats::p.adjust(p_values[, column], method = "BH")
      result[indices[column]] <- sum(is.finite(adjusted) & adjusted < fdr)
    }
  }
  result
}

bootstrap_summary <- function(values, point_values, metric, direction, method_registry) {
  rows <- do.call(rbind, lapply(seq_len(ncol(values)), function(index) {
    x <- values[, index]
    data.frame(
      method_id = method_registry$method_id[index], method = method_registry$method[index],
      metric = metric, direction = direction, point_estimate = point_values[index],
      bootstrap_mean = mean(x, na.rm = TRUE),
      ci_lower = stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE),
      ci_upper = stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE),
      valid_bootstrap = sum(is.finite(x)), stringsAsFactors = FALSE
    )
  }))
  rows
}

paired_winn_differences <- function(values, point_values, metric, direction, method_registry, winn_id = "winn_fixed_no_qc") {
  winn_index <- match(winn_id, method_registry$method_id)
  sign <- if (direction == "higher") 1 else if (direction == "lower") -1 else NA_real_
  do.call(rbind, lapply(seq_len(ncol(values)), function(index) {
    if (index == winn_index) return(NULL)
    difference <- sign * (values[, winn_index] - values[, index])
    data.frame(
      metric = metric, direction = direction,
      winn_method_id = winn_id, winn_method = method_registry$method[winn_index],
      comparator_method_id = method_registry$method_id[index], comparator_method = method_registry$method[index],
      point_direction_normalized_difference = sign * (point_values[winn_index] - point_values[index]),
      bootstrap_mean_difference = mean(difference, na.rm = TRUE),
      ci_lower = stats::quantile(difference, 0.025, na.rm = TRUE, names = FALSE),
      ci_upper = stats::quantile(difference, 0.975, na.rm = TRUE, names = FALSE),
      probability_winn_advantage = mean(difference > 0, na.rm = TRUE),
      valid_bootstrap = sum(is.finite(difference)),
      difference_definition = "Positive values favor fixed/no-QC WiNN after metric-direction normalization.",
      stringsAsFactors = FALSE
    )
  }))
}

cluster_bootstrap_draws <- function(units, strata = NULL, n_bootstrap, seed) {
  units <- as.character(units)
  if (is.null(strata)) strata <- rep("all", length(units))
  strata <- as.character(strata)
  if (length(strata) != length(units) || anyDuplicated(units)) stop("Invalid independent-unit bootstrap specification.", call. = FALSE)
  set.seed(seed)
  replicate(n_bootstrap, {
    unlist(lapply(unique(strata), function(stratum) {
      available <- units[strata == stratum]
      sample(available, length(available), replace = TRUE)
    }), use.names = FALSE)
  }, simplify = FALSE)
}

icc_a1_profiles <- function(matrix) {
  matrix <- as.matrix(matrix)
  keep <- apply(matrix, 1L, function(value) all(is.finite(value)))
  matrix <- matrix[keep, , drop = FALSE]
  n <- nrow(matrix); k <- ncol(matrix)
  if (n < 3L || k < 2L) return(NA_real_)
  row_means <- rowMeans(matrix); column_means <- colMeans(matrix); grand <- mean(matrix)
  ss_rows <- k * sum((row_means - grand)^2)
  ss_columns <- n * sum((column_means - grand)^2)
  ss_error <- sum((matrix - row_means - rep(column_means, each = n) + grand)^2)
  ms_rows <- ss_rows / (n - 1L); ms_columns <- ss_columns / (k - 1L)
  ms_error <- ss_error / ((n - 1L) * (k - 1L))
  denominator <- ms_rows + (k - 1L) * ms_error + k * (ms_columns - ms_error) / n
  if (!is.finite(denominator) || denominator <= 0) NA_real_ else (ms_rows - ms_error) / denominator
}

replicate_identity_metrics <- function(log_matrix, metadata, identity_field, batch_field = "batch") {
  identities <- as.character(metadata[[identity_field]])
  batches <- as.character(metadata[[batch_field]])
  groups <- split(seq_along(identities), identities)
  groups <- groups[lengths(groups) >= 2L]
  do.call(rbind, lapply(names(groups), function(identity) {
    index <- groups[[identity]]
    pairs <- utils::combn(index, 2L, simplify = FALSE)
    cross_pairs <- Filter(function(pair) batches[pair[1]] != batches[pair[2]], pairs)
    if (!length(cross_pairs)) cross_pairs <- pairs
    correlations <- vapply(cross_pairs, function(pair) {
      suppressWarnings(stats::cor(log_matrix[, pair[1]], log_matrix[, pair[2]]))
    }, numeric(1))
    data.frame(
      identity = identity, n_replicates = length(index), n_cross_batch_pairs = length(cross_pairs),
      pearson = stats::median(correlations, na.rm = TRUE),
      icc_a1 = icc_a1_profiles(log_matrix[, index, drop = FALSE]),
      stringsAsFactors = FALSE
    )
  }))
}

cross_batch_matching <- function(log_matrix, metadata, identity_field, batch_field = "batch") {
  identity <- as.character(metadata[[identity_field]])
  batch <- as.character(metadata[[batch_field]])
  similarity <- suppressWarnings(stats::cor(log_matrix, use = "pairwise.complete.obs"))
  diag(similarity) <- NA_real_
  queries <- lapply(seq_along(identity), function(index) {
    eligible <- which(batch != batch[index] & is.finite(similarity[index, ]))
    matching <- eligible[identity[eligible] == identity[index]]
    if (!length(eligible) || !length(matching)) return(NULL)
    ordered <- eligible[order(similarity[index, eligible], decreasing = TRUE)]
    first_match <- match(TRUE, identity[ordered] == identity[index])
    data.frame(
      query_index = index, identity = identity[index], batch = batch[index],
      top1 = as.numeric(identity[ordered[1]] == identity[index]),
      reciprocal_rank = 1 / first_match,
      best_match_similarity = max(similarity[index, matching], na.rm = TRUE),
      best_unmatched_similarity = max(similarity[index, eligible[identity[eligible] != identity[index]]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  queries <- do.call(rbind, Filter(Negate(is.null), queries))
  per_identity <- stats::aggregate(
    queries[c("top1", "reciprocal_rank", "best_match_similarity", "best_unmatched_similarity")],
    list(identity = queries$identity), mean, na.rm = TRUE
  )
  list(query = queries, identity = per_identity)
}

feature_repeatability_ratio <- function(log_matrix, identities) {
  groups <- split(seq_along(identities), as.character(identities))
  groups <- groups[lengths(groups) >= 2L]
  vapply(seq_len(nrow(log_matrix)), function(feature) {
    values <- log_matrix[feature, ]
    means <- vapply(groups, function(index) mean(values[index]), numeric(1))
    between <- stats::var(means)
    within_df <- sum(lengths(groups) - 1L)
    within <- sum(vapply(groups, function(index) stats::var(values[index]) * (length(index) - 1L), numeric(1))) / within_df
    if (!is.finite(between) || !is.finite(within) || between + within <= 0) NA_real_ else between / (between + within)
  }, numeric(1))
}

effect_concordance <- function(training_effect, heldout_effect, discovered = rep(TRUE, length(training_effect))) {
  keep <- is.finite(training_effect) & is.finite(heldout_effect) & discovered
  if (sum(keep) < 3L) {
    return(c(pearson = NA, spearman = NA, sign_concordance = NA, attenuation = NA, n_features = sum(keep)))
  }
  c(
    pearson = stats::cor(training_effect[keep], heldout_effect[keep], method = "pearson"),
    spearman = stats::cor(training_effect[keep], heldout_effect[keep], method = "spearman"),
    sign_concordance = mean(sign(training_effect[keep]) == sign(heldout_effect[keep])),
    attenuation = stats::median(abs(heldout_effect[keep]) / pmax(abs(training_effect[keep]), 1e-8)),
    n_features = sum(keep)
  )
}

ols_feature_model <- function(log_matrix, design, coefficient_index, fdr = 0.05) {
  design <- as.matrix(design)
  if (qr(design)$rank < ncol(design)) return(NULL)
  inverse <- solve(crossprod(design))
  cross <- crossprod(design, t(log_matrix))
  coefficients <- inverse %*% cross
  fitted <- design %*% coefficients
  residuals <- t(log_matrix) - fitted
  df <- nrow(design) - ncol(design)
  residual_variance <- colSums(residuals^2) / df
  standard_error <- sqrt(residual_variance * inverse[coefficient_index, coefficient_index])
  statistic <- coefficients[coefficient_index, ] / standard_error
  p_value <- 2 * stats::pt(-abs(statistic), df = df)
  p_adjusted <- stats::p.adjust(p_value, method = "BH")
  outcome_sd <- matrixStats::rowSds(log_matrix)
  list(
    table = data.frame(
      feature_id = rownames(log_matrix), coefficient = coefficients[coefficient_index, ],
      standardized_effect = coefficients[coefficient_index, ] / outcome_sd,
      residual_variance = residual_variance, t_statistic = statistic,
      p_value = p_value, p_adj = p_adjusted,
      significant = is.finite(p_adjusted) & p_adjusted < fdr,
      stringsAsFactors = FALSE
    ),
    coefficients = coefficients, fitted = fitted, residuals = residuals,
    inverse = inverse, df = df
  )
}

freedman_lane_counts <- function(log_matrix, full_design, reduced_design, coefficient_index, permutations, fdr = 0.05) {
  reduced <- ols_feature_model(log_matrix, reduced_design, 1L, fdr)
  if (is.null(reduced) || qr(full_design)$rank < ncol(full_design)) return(rep(NA_integer_, length(permutations)))
  full_inverse <- solve(crossprod(full_design))
  result <- integer(length(permutations))
  for (index in seq_along(permutations)) {
    permuted_response <- reduced$fitted + reduced$residuals[permutations[[index]], , drop = FALSE]
    cross <- crossprod(full_design, permuted_response)
    coefficients <- full_inverse %*% cross
    residual_ss <- colSums(permuted_response^2) - colSums(coefficients * cross)
    residual_variance <- pmax(residual_ss / (nrow(full_design) - ncol(full_design)), 0)
    standard_error <- sqrt(residual_variance * full_inverse[coefficient_index, coefficient_index])
    statistic <- coefficients[coefficient_index, ] / standard_error
    p_value <- 2 * stats::pt(-abs(statistic), df = nrow(full_design) - ncol(full_design))
    adjusted <- stats::p.adjust(p_value, method = "BH")
    result[index] <- sum(is.finite(adjusted) & adjusted < fdr)
  }
  result
}

restricted_index_permutations <- function(blocks, n_permutations, seed) {
  blocks <- as.character(blocks)
  set.seed(seed)
  replicate(n_permutations, {
    result <- seq_along(blocks)
    for (block in unique(blocks)) {
      index <- which(blocks == block)
      result[index] <- sample(index, length(index), replace = FALSE)
    }
    result
  }, simplify = FALSE)
}

permutation_summary <- function(null_counts, observed_count, method_id, method, endpoint) {
  data.frame(
    method_id = method_id, method = method, endpoint = endpoint,
    observed_count = observed_count,
    null_median = stats::median(null_counts, na.rm = TRUE),
    null_95th_percentile = stats::quantile(null_counts, 0.95, na.rm = TRUE, names = FALSE),
    observed_minus_null_median = observed_count - stats::median(null_counts, na.rm = TRUE),
    empirical_p_value = (1 + sum(null_counts >= observed_count, na.rm = TRUE)) /
      (1 + sum(is.finite(null_counts))),
    valid_permutations = sum(is.finite(null_counts)),
    stringsAsFactors = FALSE
  )
}

hash_private_identity <- function(value, salt = "winn-derivative-clinical-v1") {
  vapply(as.character(value), function(item) {
    substr(digest::digest(paste(salt, item, sep = "::"), algo = "sha256", serialize = FALSE), 1, 16)
  }, character(1))
}
