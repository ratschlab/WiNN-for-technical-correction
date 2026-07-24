# Pure construction helpers for the paired partial-confounding simulation.
#
# This file performs no file I/O and runs no correction method.  It extends the
# canonical technical-artifact simulator by adding a phenotype to study
# injections before the unchanged batch and drift artifacts are applied.

.partial_confounding_axes <- c("batch", "run_order", "joint")
.partial_confounding_positive_strengths <- c(0.25, 0.50, 0.75, 0.90, 1.00)
.partial_confounding_required_metadata <- c(
  "sample_id", "plate", "order_in_plate", "run_order", "sample_type"
)
.partial_confounding_correction_fields <- c(
  "sample_id", "plate", "order_in_plate", "run_order", "sample_type"
)
.partial_confounding_forbidden_correction_fields <- c(
  "phenotype", "phenotype_numeric", "is_case", "responsive",
  "effect_direction", "phenotype_effect_log"
)

.partial_preserve_rng <- function(expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(expr)
}

#' Return the deduplicated partial-confounding grid for one master seed
#'
#' Strength zero is represented once because all three axes have the same
#' random, batch-balanced allocation at that point.  The result therefore has
#' 16 scenarios: one zero condition plus 3 axes x 5 positive strengths.
partial_confounding_scenario_grid <- function(seed_id) {
  if (length(seed_id) != 1L || is.na(seed_id) || !grepl("^CONF[0-9]{2}$", seed_id)) {
    stop("seed_id must have form CONF01, ..., CONF99.", call. = FALSE)
  }
  positive <- expand.grid(
    confounding_axis = .partial_confounding_axes,
    nominal_strength = .partial_confounding_positive_strengths,
    stringsAsFactors = FALSE
  )
  positive <- positive[order(
    match(positive$nominal_strength, .partial_confounding_positive_strengths),
    match(positive$confounding_axis, .partial_confounding_axes)
  ), , drop = FALSE]
  grid <- rbind(
    data.frame(
      confounding_axis = "none",
      nominal_strength = 0,
      stringsAsFactors = FALSE
    ),
    positive
  )
  strength_tag <- gsub("[.]", "p", sprintf("%.2f", grid$nominal_strength))
  grid$seed_id <- seed_id
  grid$scenario_id <- sprintf(
    "%s_%s_%s", seed_id, grid$confounding_axis, strength_tag
  )
  grid$scenario_index <- seq_len(nrow(grid))
  grid <- grid[, c(
    "seed_id", "scenario_id", "scenario_index", "confounding_axis",
    "nominal_strength"
  )]
  rownames(grid) <- NULL
  stopifnot(nrow(grid) == 16L, sum(grid$nominal_strength == 0) == 1L)
  grid
}

#' Expand a prespecified CONF master seed into independent component seeds
#'
#' Fixed offsets make the expansion transparent and invariant to scenario
#' ordering.  All scenarios within a master seed reuse every component seed.
expand_partial_confounding_seed <- function(master_seed) {
  master_seed <- as.integer(master_seed)
  if (length(master_seed) != 1L || is.na(master_seed) || master_seed <= 0L) {
    stop("master_seed must be one positive integer.", call. = FALSE)
  }
  offsets <- c(
    clean_truth = 101L,
    batch_assignment = 211L,
    batch_magnitude = 307L,
    drift_assignment = 401L,
    drift_morphology = 503L,
    drift_shape = 601L,
    responsive_features = 701L,
    allocation_priority = 809L,
    method_Raw = 1201L,
    method_ComBat = 1301L,
    method_QC_RLSC = 1401L,
    method_QC_RFSC = 1501L,
    method_TIGER = 1601L,
    method_SERRF = 1701L,
    method_WINN_auto_QC = 1801L,
    method_WINN_auto_batch_QC = 1901L,
    method_WINN_fixed_no_QC = 2001L
  )
  seeds <- master_seed + offsets
  if (any(seeds > .Machine$integer.max)) {
    stop("Expanded seed exceeds the R integer range.", call. = FALSE)
  }
  stats::setNames(as.integer(seeds), names(offsets))
}

#' Select fixed phenotype-responsive features and effect directions
partial_confounding_feature_truth <- function(
  feature_ids,
  seed,
  responsive_fraction = 0.10,
  effect_magnitude_log = 0.35
) {
  feature_ids <- as.character(feature_ids)
  seed <- as.integer(seed)
  if (!length(feature_ids) || anyNA(feature_ids) || anyDuplicated(feature_ids)) {
    stop("feature_ids must be non-missing and unique.", call. = FALSE)
  }
  n_responsive <- as.integer(round(length(feature_ids) * responsive_fraction))
  if (n_responsive < 2L || n_responsive %% 2L != 0L) {
    stop("The responsive-feature count must be positive and even.", call. = FALSE)
  }
  if (!is.finite(effect_magnitude_log) || effect_magnitude_log <= 0) {
    stop("effect_magnitude_log must be positive and finite.", call. = FALSE)
  }

  .partial_preserve_rng({
    set.seed(seed)
    responsive_ids <- sample(feature_ids, n_responsive, replace = FALSE)
    directions <- sample(rep(c(-1L, 1L), each = n_responsive / 2L))
    truth <- data.frame(
      metabolite = feature_ids,
      responsive = feature_ids %in% responsive_ids,
      effect_direction = 0L,
      phenotype_effect_log = 0,
      stringsAsFactors = FALSE
    )
    match_index <- match(responsive_ids, truth$metabolite)
    truth$effect_direction[match_index] <- directions
    truth$phenotype_effect_log[match_index] <-
      directions * effect_magnitude_log
    stopifnot(
      sum(truth$responsive) == n_responsive,
      sum(truth$effect_direction == 1L) == n_responsive / 2L,
      sum(truth$effect_direction == -1L) == n_responsive / 2L
    )
    truth
  })
}

#' Add study-order ranks, scaled positions, and quartiles to metadata
partial_confounding_order_metadata <- function(metadata) {
  missing_fields <- setdiff(.partial_confounding_required_metadata, names(metadata))
  if (length(missing_fields)) {
    stop(
      "Metadata lacks field(s): ", paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(metadata$sample_id)) {
    stop("sample_id must be unique.", call. = FALSE)
  }
  metadata <- metadata[order(metadata$run_order), , drop = FALSE]
  metadata$study_order_rank <- NA_integer_
  metadata$within_plate_position_scaled <- NA_real_
  metadata$run_order_quartile <- NA_character_
  plate_levels <- unique(metadata$plate)
  for (plate in plate_levels) {
    index <- which(metadata$plate == plate & metadata$sample_type == "sample")
    index <- index[order(metadata$order_in_plate[index])]
    n <- length(index)
    metadata$study_order_rank[index] <- seq_len(n)
    metadata$within_plate_position_scaled[index] <-
      if (n > 1L) (seq_len(n) - 1) / (n - 1) else 0
    # Four deterministic, nearly equal groups; for n=90 their sizes are
    # 23, 22, 22, and 23.
    metadata$run_order_quartile[index] <- paste0(
      "Q", pmin(4L, floor((seq_len(n) - 1L) * 4L / n) + 1L)
    )
  }
  rownames(metadata) <- NULL
  metadata
}

#' Exact batch case counts at one nominal alignment strength
#'
#' The maximum plate case-fraction pattern is (1, 1, .5, 0, 0).  Counts are
#' interpolated from the balanced vector and rounded under the constraints that
#' paired extreme plates remain symmetric and the total remains exactly 225.
partial_confounding_batch_case_counts <- function(
  strength,
  study_counts = rep(90L, 5L),
  maximum_case_fractions = c(1, 1, 0.5, 0, 0)
) {
  strength <- as.numeric(strength)
  study_counts <- as.integer(study_counts)
  if (length(strength) != 1L || is.na(strength) || strength < 0 || strength > 1) {
    stop("strength must be a scalar in [0, 1].", call. = FALSE)
  }
  if (
    length(study_counts) != 5L || any(study_counts != 90L) ||
      !identical(as.numeric(maximum_case_fractions), c(1, 1, 0.5, 0, 0))
  ) {
    stop(
      "This prespecified design requires five plates with 90 study samples and ",
      "maximum fractions (1, 1, .5, 0, 0).",
      call. = FALSE
    )
  }
  balanced_count <- study_counts[[1L]] / 2L
  ideal_high <- balanced_count * (1 + strength)
  high <- as.integer(round(ideal_high))
  low <- as.integer(study_counts[[1L]] - high)
  counts <- c(high, high, as.integer(balanced_count), low, low)
  ideal <- (1 - strength) * (study_counts * 0.5) +
    strength * (study_counts * maximum_case_fractions)
  attr(counts, "ideal") <- ideal
  stopifnot(sum(counts) == 225L, all(counts >= 0L), all(counts <= study_counts))
  counts
}

.partial_progress_to_early_cases <- function(labels, order_rank, strength) {
  labels <- as.integer(labels)
  order_rank <- as.integer(order_rank)
  n_cases <- sum(labels == 1L)
  if (!n_cases || n_cases == length(labels) || strength <= 0) return(labels)
  target <- integer(length(labels))
  target[order(order_rank)[seq_len(n_cases)]] <- 1L
  missing_early <- which(target == 1L & labels == 0L)
  extra_late <- which(target == 0L & labels == 1L)
  stopifnot(length(missing_early) == length(extra_late))
  n_swap <- as.integer(round(strength * length(missing_early)))
  if (!n_swap) return(labels)
  missing_early <- missing_early[order(order_rank[missing_early])]
  extra_late <- extra_late[order(order_rank[extra_late], decreasing = TRUE)]
  labels[missing_early[seq_len(n_swap)]] <- 1L
  labels[extra_late[seq_len(n_swap)]] <- 0L
  stopifnot(sum(labels) == n_cases)
  labels
}

#' Generate fixed random allocation priorities for one master seed
partial_confounding_allocation_priority <- function(metadata, seed) {
  metadata <- partial_confounding_order_metadata(metadata)
  study_index <- which(metadata$sample_type == "sample")
  .partial_preserve_rng({
    set.seed(as.integer(seed))
    priority <- rep(NA_real_, nrow(metadata))
    priority[study_index] <- stats::runif(length(study_index))
    stats::setNames(priority, metadata$sample_id)
  })
}

#' Allocate one paired phenotype scenario
allocate_partial_confounding_phenotype <- function(
  metadata,
  allocation_priority,
  confounding_axis,
  nominal_strength
) {
  metadata <- partial_confounding_order_metadata(metadata)
  confounding_axis <- match.arg(
    confounding_axis,
    choices = c("none", .partial_confounding_axes)
  )
  nominal_strength <- as.numeric(nominal_strength)
  if (
    length(nominal_strength) != 1L || is.na(nominal_strength) ||
      nominal_strength < 0 || nominal_strength > 1
  ) {
    stop("nominal_strength must be a scalar in [0, 1].", call. = FALSE)
  }
  if (confounding_axis == "none" && nominal_strength != 0) {
    stop("The 'none' axis is valid only at strength zero.", call. = FALSE)
  }
  allocation_priority <- allocation_priority[metadata$sample_id]
  if (
    length(allocation_priority) != nrow(metadata) ||
      anyNA(allocation_priority[metadata$sample_type == "sample"])
  ) {
    stop("allocation_priority does not cover every study injection.", call. = FALSE)
  }

  plate_levels <- unique(metadata$plate)
  study_counts <- vapply(plate_levels, function(plate) {
    sum(metadata$plate == plate & metadata$sample_type == "sample")
  }, integer(1))
  use_batch_alignment <- confounding_axis %in% c("batch", "joint")
  target_case_counts <- if (use_batch_alignment) {
    partial_confounding_batch_case_counts(nominal_strength, study_counts)
  } else {
    as.integer(study_counts / 2L)
  }

  phenotype_numeric <- rep(NA_integer_, nrow(metadata))
  for (plate_index in seq_along(plate_levels)) {
    plate <- plate_levels[[plate_index]]
    index <- which(metadata$plate == plate & metadata$sample_type == "sample")
    k <- target_case_counts[[plate_index]]
    labels <- integer(length(index))
    if (k > 0L) {
      labels[order(allocation_priority[index])[seq_len(k)]] <- 1L
    }
    if (confounding_axis %in% c("run_order", "joint")) {
      labels <- .partial_progress_to_early_cases(
        labels,
        metadata$study_order_rank[index],
        nominal_strength
      )
    }
    phenotype_numeric[index] <- labels
  }

  phenotype <- ifelse(
    metadata$sample_type == "control", "QC",
    ifelse(phenotype_numeric == 1L, "case", "control")
  )
  allocation <- data.frame(
    sample_id = metadata$sample_id,
    plate = metadata$plate,
    sample_type = metadata$sample_type,
    run_order = metadata$run_order,
    order_in_plate = metadata$order_in_plate,
    study_order_rank = metadata$study_order_rank,
    within_plate_position_scaled = metadata$within_plate_position_scaled,
    run_order_quartile = metadata$run_order_quartile,
    allocation_priority = unname(allocation_priority),
    phenotype = phenotype,
    phenotype_numeric = phenotype_numeric,
    stringsAsFactors = FALSE
  )
  observed <- table(
    factor(allocation$plate[allocation$sample_type == "sample"], levels = plate_levels),
    factor(
      allocation$phenotype_numeric[allocation$sample_type == "sample"],
      levels = c(0L, 1L)
    )
  )
  stopifnot(
    sum(allocation$phenotype_numeric == 1L, na.rm = TRUE) == 225L,
    sum(allocation$phenotype_numeric == 0L, na.rm = TRUE) == 225L,
    all(is.na(allocation$phenotype_numeric[allocation$sample_type == "control"])),
    identical(as.integer(observed[, 2L]), as.integer(target_case_counts))
  )
  attr(allocation, "target_case_counts") <- stats::setNames(
    target_case_counts, plate_levels
  )
  allocation
}

.partial_cramers_v <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  expected <- outer(rowSums(tab), colSums(tab)) / n
  statistic <- sum((tab - expected)^2 / expected)
  denominator <- min(nrow(tab) - 1L, ncol(tab) - 1L)
  if (!is.finite(denominator) || denominator <= 0L) return(NA_real_)
  sqrt(statistic / (n * denominator))
}

.partial_safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || stats::sd(x[keep]) == 0 || stats::sd(y[keep]) == 0) {
    return(NA_real_)
  }
  stats::cor(x[keep], y[keep], method = "pearson")
}

#' Compute realized confounding and estimability diagnostics
partial_confounding_diagnostics <- function(allocation) {
  study <- allocation[allocation$sample_type == "sample", , drop = FALSE]
  phenotype_factor <- factor(study$phenotype, levels = c("control", "case"))
  plate_factor <- factor(study$plate, levels = unique(allocation$plate))
  pooled_cor <- .partial_safe_cor(
    study$phenotype_numeric, study$within_plate_position_scaled
  )
  per_plate_cor <- vapply(levels(plate_factor), function(plate) {
    index <- study$plate == plate
    .partial_safe_cor(
      study$phenotype_numeric[index],
      study$within_plate_position_scaled[index]
    )
  }, numeric(1))

  design <- stats::model.matrix(
    ~ phenotype_factor + plate_factor + within_plate_position_scaled,
    data = study
  )
  singular_values <- svd(design, nu = 0L, nv = 0L)$d
  design_rank <- qr(design)$rank
  condition_number <- if (min(singular_values) <= 0) {
    Inf
  } else {
    max(singular_values) / min(singular_values)
  }

  summary <- data.frame(
    n_study = nrow(study),
    n_case = sum(study$phenotype_numeric == 1L),
    n_control = sum(study$phenotype_numeric == 0L),
    cramers_v_phenotype_batch = .partial_cramers_v(
      phenotype_factor, plate_factor
    ),
    abs_cor_phenotype_within_batch_order = abs(pooled_cor),
    weighted_mean_abs_within_plate_cor = stats::weighted.mean(
      abs(per_plate_cor),
      w = as.numeric(table(plate_factor)),
      na.rm = TRUE
    ),
    max_abs_within_plate_cor = if (all(is.na(per_plate_cor))) {
      NA_real_
    } else {
      max(abs(per_plate_cor), na.rm = TRUE)
    },
    design_formula = "~ phenotype + plate + within_plate_position_scaled",
    design_nrow = nrow(design),
    design_ncol = ncol(design),
    design_rank = design_rank,
    design_full_rank = design_rank == ncol(design),
    design_condition_number = condition_number,
    design_min_singular_value = min(singular_values),
    stringsAsFactors = FALSE
  )

  plate_correlations <- data.frame(
    plate = names(per_plate_cor),
    pearson_cor_phenotype_order = unname(per_plate_cor),
    abs_pearson_cor_phenotype_order = abs(unname(per_plate_cor)),
    stringsAsFactors = FALSE
  )

  batch_grid <- expand.grid(
    plate = unique(allocation$plate),
    phenotype = c("control", "case"),
    stringsAsFactors = FALSE
  )
  batch_observed <- as.data.frame(table(study$plate, study$phenotype),
                                  stringsAsFactors = FALSE)
  names(batch_observed) <- c("plate", "phenotype", "n")
  batch_counts <- merge(
    batch_grid, batch_observed,
    by = c("plate", "phenotype"), all.x = TRUE, sort = FALSE
  )
  batch_counts$n[is.na(batch_counts$n)] <- 0L
  batch_counts <- batch_counts[order(
    match(batch_counts$plate, unique(allocation$plate)),
    match(batch_counts$phenotype, c("control", "case"))
  ), , drop = FALSE]

  quantile_grid <- expand.grid(
    plate = unique(allocation$plate),
    run_order_quartile = paste0("Q", 1:4),
    phenotype = c("control", "case"),
    stringsAsFactors = FALSE
  )
  quantile_observed <- as.data.frame(table(
    study$plate,
    factor(study$run_order_quartile, levels = paste0("Q", 1:4)),
    factor(study$phenotype, levels = c("control", "case"))
  ), stringsAsFactors = FALSE)
  names(quantile_observed) <- c(
    "plate", "run_order_quartile", "phenotype", "n"
  )
  quantile_counts <- merge(
    quantile_grid, quantile_observed,
    by = c("plate", "run_order_quartile", "phenotype"),
    all.x = TRUE, sort = FALSE
  )
  quantile_counts$n[is.na(quantile_counts$n)] <- 0L
  quantile_counts <- quantile_counts[order(
    match(quantile_counts$plate, unique(allocation$plate)),
    match(quantile_counts$run_order_quartile, paste0("Q", 1:4)),
    match(quantile_counts$phenotype, c("control", "case"))
  ), , drop = FALSE]

  list(
    summary = summary,
    within_plate_correlations = plate_correlations,
    batch_group_counts = batch_counts,
    run_order_quantile_counts = quantile_counts
  )
}

#' Construct clean and observed matrices for one phenotype allocation
#'
#' The case-control contrast is exactly +/-0.35 on the log scale: centered
#' phenotype coding (-0.5, +0.5) gives a contrast of one coefficient while
#' keeping the balanced-study grand mean unchanged.  Pooled QCs receive zero
#' biological effect.  `correction_metadata` deliberately excludes phenotype
#' and response truth.
construct_partial_confounding_matrices <- function(
  canonical_base,
  allocation,
  feature_truth,
  include_intensities = TRUE
) {
  required_artifacts <- c(
    "clean_log", "batch_shift_by_injection", "drift_by_injection"
  )
  if (
    is.null(canonical_base$artifact_matrices) ||
      length(setdiff(required_artifacts, names(canonical_base$artifact_matrices)))
  ) {
    stop(
      "canonical_base must be generated with include_artifact_matrices = TRUE.",
      call. = FALSE
    )
  }
  metadata <- canonical_base$metadata
  if (!identical(metadata$sample_id, allocation$sample_id)) {
    stop("Allocation order does not match canonical metadata.", call. = FALSE)
  }
  clean_log_baseline <- canonical_base$artifact_matrices$clean_log
  beta <- feature_truth$phenotype_effect_log[
    match(rownames(clean_log_baseline), feature_truth$metabolite)
  ]
  if (anyNA(beta)) stop("Feature truth does not cover the canonical matrix.", call. = FALSE)
  centered_phenotype <- allocation$phenotype_numeric - 0.5
  centered_phenotype[allocation$sample_type == "control"] <- 0
  biological_effect_log <- tcrossprod(beta, centered_phenotype)
  dimnames(biological_effect_log) <- dimnames(clean_log_baseline)
  technical_artifact_log <-
    canonical_base$artifact_matrices$batch_shift_by_injection +
    canonical_base$artifact_matrices$drift_by_injection
  clean_log_with_biology <- clean_log_baseline + biological_effect_log
  raw_log_with_biology <- clean_log_with_biology + technical_artifact_log

  correction_metadata <- metadata[, .partial_confounding_correction_fields,
                                  drop = FALSE]
  forbidden_present <- intersect(
    names(correction_metadata), .partial_confounding_forbidden_correction_fields
  )
  stopifnot(
    !length(forbidden_present),
    all(biological_effect_log[, allocation$sample_type == "control", drop = FALSE] == 0)
  )

  output <- list(
    correction_metadata = correction_metadata,
    evaluation_truth = list(
      phenotype = allocation[, c(
        "sample_id", "sample_type", "phenotype", "phenotype_numeric"
      )],
      feature_truth = feature_truth
    ),
    biological_effect_log = biological_effect_log,
    technical_artifact_log = technical_artifact_log,
    clean_log_with_biology = clean_log_with_biology,
    raw_log_with_biology = raw_log_with_biology
  )
  if (isTRUE(include_intensities)) {
    output$clean_ground_truth <- exp(clean_log_with_biology)
    output$raw_intensity <- exp(raw_log_with_biology)
  }
  output
}

#' Validate invariants for one constructed scenario
validate_partial_confounding_scenario <- function(
  canonical_base,
  allocation,
  feature_truth,
  constructed,
  tolerance = 1e-12
) {
  study <- allocation$sample_type == "sample"
  qc <- allocation$sample_type == "control"
  expected_artifact <-
    canonical_base$artifact_matrices$batch_shift_by_injection +
    canonical_base$artifact_matrices$drift_by_injection
  artifact_difference <- max(abs(
    (constructed$raw_log_with_biology - constructed$clean_log_with_biology) -
      expected_artifact
  ))
  checks <- data.frame(
    check = c(
      "study_balance_225_225",
      "qcs_phenotype_neutral",
      "responsive_features_exactly_100",
      "effect_directions_balanced",
      "effect_magnitude_abs_0p35",
      "technical_artifact_identical",
      "phenotype_excluded_from_correction_metadata",
      "correction_metadata_order_aligned",
      "all_constructed_logs_finite"
    ),
    observed = c(
      sprintf(
        "%d cases; %d controls",
        sum(allocation$phenotype_numeric[study] == 1L),
        sum(allocation$phenotype_numeric[study] == 0L)
      ),
      sprintf("%d non-neutral QC labels/effects", sum(
        !is.na(allocation$phenotype_numeric[qc]) |
          colSums(abs(constructed$biological_effect_log[, qc, drop = FALSE])) > 0
      )),
      as.character(sum(feature_truth$responsive)),
      sprintf(
        "%d positive; %d negative",
        sum(feature_truth$effect_direction == 1L),
        sum(feature_truth$effect_direction == -1L)
      ),
      paste(sort(unique(abs(
        feature_truth$phenotype_effect_log[feature_truth$responsive]
      ))), collapse = ","),
      format(artifact_difference, scientific = TRUE),
      paste(names(constructed$correction_metadata), collapse = ";"),
      as.character(identical(
        constructed$correction_metadata$sample_id,
        canonical_base$metadata$sample_id
      )),
      as.character(
        all(is.finite(constructed$clean_log_with_biology)) &&
          all(is.finite(constructed$raw_log_with_biology))
      )
    ),
    passed = c(
      sum(allocation$phenotype_numeric[study] == 1L) == 225L &&
        sum(allocation$phenotype_numeric[study] == 0L) == 225L,
      all(is.na(allocation$phenotype_numeric[qc])) &&
        all(constructed$biological_effect_log[, qc, drop = FALSE] == 0),
      sum(feature_truth$responsive) == 100L,
      sum(feature_truth$effect_direction == 1L) == 50L &&
        sum(feature_truth$effect_direction == -1L) == 50L,
      all(abs(feature_truth$phenotype_effect_log[feature_truth$responsive]) == 0.35),
      is.finite(artifact_difference) && artifact_difference <= tolerance,
      !length(intersect(
        names(constructed$correction_metadata),
        .partial_confounding_forbidden_correction_fields
      )),
      identical(
        constructed$correction_metadata$sample_id,
        canonical_base$metadata$sample_id
      ),
      all(is.finite(constructed$clean_log_with_biology)) &&
        all(is.finite(constructed$raw_log_with_biology))
    ),
    stringsAsFactors = FALSE
  )
  checks
}
