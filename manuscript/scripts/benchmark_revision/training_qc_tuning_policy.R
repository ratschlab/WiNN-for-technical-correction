# Frozen policy for the training-QC-only tuning sensitivity analysis.
#
# This file is deliberately independent of the benchmark notebooks.  It may be
# sourced by the eventual rerun driver, but it must not read held-out-reference
# values, phenotype/covariate labels, biological identities, or final
# preservation endpoints during candidate selection.

training_qc_profile_score <- function(mat, validation_qc_ids) {
  mat <- as.matrix(mat)
  ids <- as.character(validation_qc_ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  ids <- ids[ids %in% colnames(mat)]

  if (length(ids) < 2L) {
    return(data.frame(
      n_validation_qc = length(ids),
      n_features_input = nrow(mat),
      n_features_scored = 0L,
      mean_training_qc_profile_correlation = NA_real_,
      mean_training_qc_cv_fraction = NA_real_,
      qc_only_score = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  q <- mat[, ids, drop = FALSE]
  finite_rows <- rowSums(is.finite(q)) >= 2L
  q_score <- q[finite_rows, , drop = FALSE]

  if (!nrow(q_score)) {
    return(data.frame(
      n_validation_qc = length(ids),
      n_features_input = nrow(mat),
      n_features_scored = 0L,
      mean_training_qc_profile_correlation = NA_real_,
      mean_training_qc_cv_fraction = NA_real_,
      qc_only_score = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  cvs <- apply(q_score, 1L, function(v) {
    mu <- mean(v, na.rm = TRUE)
    if (!is.finite(mu) || abs(mu) <= .Machine$double.eps) {
      return(NA_real_)
    }
    stats::sd(v, na.rm = TRUE) / abs(mu)
  })
  correlations <- suppressWarnings(stats::cor(q_score, use = "pairwise.complete.obs"))
  correlation_values <- correlations[lower.tri(correlations)]
  mean_correlation <- mean(correlation_values, na.rm = TRUE)
  mean_cv <- mean(cvs, na.rm = TRUE)
  if (!is.finite(mean_correlation)) mean_correlation <- NA_real_
  if (!is.finite(mean_cv)) mean_cv <- NA_real_

  data.frame(
    n_validation_qc = length(ids),
    n_features_input = nrow(mat),
    n_features_scored = sum(is.finite(cvs)),
    mean_training_qc_profile_correlation = mean_correlation,
    mean_training_qc_cv_fraction = mean_cv,
    qc_only_score = mean_correlation - 0.5 * mean_cv,
    stringsAsFactors = FALSE
  )
}

aggregate_training_qc_fold_scores <- function(fold_scores) {
  required <- c(
    "candidate_id", "qc_only_score",
    "mean_training_qc_profile_correlation",
    "mean_training_qc_cv_fraction", "feature_coverage"
  )
  missing_columns <- setdiff(required, names(fold_scores))
  if (length(missing_columns)) {
    stop("Missing fold-score columns: ", paste(missing_columns, collapse = ", "))
  }

  candidate_ids <- unique(fold_scores$candidate_id)
  rows <- lapply(candidate_ids, function(candidate_id) {
    d <- fold_scores[fold_scores$candidate_id == candidate_id, , drop = FALSE]
    data.frame(
      candidate_id = candidate_id,
      n_folds_attempted = nrow(d),
      n_folds_scored = sum(is.finite(d$qc_only_score)),
      qc_only_score = mean(d$qc_only_score, na.rm = TRUE),
      mean_training_qc_profile_correlation = mean(
        d$mean_training_qc_profile_correlation,
        na.rm = TRUE
      ),
      mean_training_qc_cv_fraction = mean(
        d$mean_training_qc_cv_fraction,
        na.rm = TRUE
      ),
      feature_coverage = mean(d$feature_coverage, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  for (column in c(
    "qc_only_score", "mean_training_qc_profile_correlation",
    "mean_training_qc_cv_fraction", "feature_coverage"
  )) {
    out[[column]][!is.finite(out[[column]])] <- NA_real_
  }
  out
}

select_training_qc_candidate <- function(candidate_scores) {
  required <- c(
    "candidate_id", "qc_only_score", "feature_coverage", "default_distance"
  )
  missing_columns <- setdiff(required, names(candidate_scores))
  if (length(missing_columns)) {
    stop("Missing candidate-selection columns: ", paste(missing_columns, collapse = ", "))
  }

  d <- candidate_scores[
    is.finite(candidate_scores$qc_only_score),
    ,
    drop = FALSE
  ]
  if (!nrow(d)) return(candidate_scores[0L, , drop = FALSE])

  coverage_key <- as.numeric(d$feature_coverage)
  coverage_key[!is.finite(coverage_key)] <- -Inf
  default_key <- as.numeric(d$default_distance)
  default_key[!is.finite(default_key)] <- Inf
  candidate_key <- as.integer(d$candidate_id)
  candidate_key[is.na(candidate_key)] <- .Machine$integer.max

  # Frozen deterministic tie-breaking, in the mandated order:
  # 1. higher QC-only score;
  # 2. greater feature coverage;
  # 3. configuration closest to package/default recommendation;
  # 4. stable candidate ID.
  ordering <- order(
    -as.numeric(d$qc_only_score),
    -coverage_key,
    default_key,
    candidate_key,
    na.last = TRUE
  )
  d[ordering[1L], , drop = FALSE]
}

make_repeated_stratified_qc_folds <- function(
    qc_metadata,
    id_column = "sample_id",
    batch_column = "batch",
    order_column = "run_order",
    repeat_seeds,
    n_folds = 3L,
    protect_endpoints = TRUE) {
  required <- c(id_column, batch_column, order_column)
  missing_columns <- setdiff(required, names(qc_metadata))
  if (length(missing_columns)) {
    stop("Missing QC metadata columns: ", paste(missing_columns, collapse = ", "))
  }
  if (n_folds < 2L) stop("n_folds must be at least two.")

  qc <- qc_metadata[, required, drop = FALSE]
  names(qc) <- c("qc_id", "batch", "run_order")
  qc$qc_id <- as.character(qc$qc_id)
  qc$batch <- as.character(qc$batch)
  qc$run_order <- as.numeric(qc$run_order)
  qc <- qc[order(qc$batch, qc$run_order, qc$qc_id), , drop = FALSE]

  rows <- list()
  for (repeat_index in seq_along(repeat_seeds)) {
    repeat_seed <- as.integer(repeat_seeds[[repeat_index]])
    set.seed(repeat_seed)
    for (batch_id in unique(qc$batch)) {
      d <- qc[qc$batch == batch_id, , drop = FALSE]
      protected <- rep(FALSE, nrow(d))
      if (isTRUE(protect_endpoints) && nrow(d) >= 2L) {
        protected[c(1L, nrow(d))] <- TRUE
      }
      eligible <- which(!protected)
      shuffled <- if (length(eligible)) sample(eligible, length(eligible)) else integer()
      fold_assignment <- rep(NA_integer_, nrow(d))
      if (length(shuffled)) {
        fold_assignment[shuffled] <- rep(seq_len(n_folds), length.out = length(shuffled))
      }
      rows[[length(rows) + 1L]] <- data.frame(
        repeat_id = repeat_index,
        repeat_seed = repeat_seed,
        batch = d$batch,
        qc_id = d$qc_id,
        run_order = d$run_order,
        protected_training_endpoint = protected,
        validation_fold = fold_assignment,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
