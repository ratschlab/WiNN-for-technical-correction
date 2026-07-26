# Pure generator for the canonical WiNN robustness simulation.
#
# This file deliberately performs no file I/O.  The companion
# generate_simulation_bundles.R script is responsible for provenance, hashes,
# and persistence. Keeping generation separate makes each realization easy to
# test without modifying its input configuration.

.canonical_required_seed_components <- c(
  "clean_truth",
  "batch_assignment",
  "batch_magnitude",
  "drift_assignment",
  "drift_morphology",
  "drift_shape"
)

.canonical_preserve_rng <- function(expr) {
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

.canonical_scale_to_amplitude <- function(x, amplitude) {
  x <- x - mean(x)
  largest <- max(abs(x))
  if (!is.finite(largest) || largest <= 1e-12) {
    return(rep(0, length(x)))
  }
  amplitude * (x / largest)
}

.canonical_simulate_plate_drift <- function(time, drift_type) {
  if (identical(drift_type, "linear")) {
    slope <- stats::rnorm(1L, mean = 0, sd = 1)
    shape <- slope * (time - mean(time))
    return(shape - mean(shape))
  }

  if (identical(drift_type, "nonlinear")) {
    phase <- stats::runif(1L, 0, 2 * pi)
    frequency <- sample(c(1, 2), size = 1L, prob = c(0.75, 0.25))
    sinusoid <- sin(2 * pi * frequency * time + phase)
    bump_center <- stats::runif(1L, 0.2, 0.8)
    bump_width <- stats::runif(1L, 0.08, 0.22)
    bump <- exp(-((time - bump_center)^2) / (2 * bump_width^2))
    shape <- sinusoid + stats::runif(1L, 0.5, 0.9) * bump
    return(shape - mean(shape))
  }

  if (identical(drift_type, "ar1")) {
    ar_coefficient <- stats::runif(1L, 0.45, 0.90)
    shape <- as.numeric(stats::arima.sim(
      model = list(ar = ar_coefficient), n = length(time)
    ))
    return(shape - mean(shape))
  }

  if (identical(drift_type, "mixed")) {
    slope <- stats::rnorm(1L, mean = 0, sd = 0.8)
    linear_component <- slope * (time - mean(time))

    phase <- stats::runif(1L, 0, 2 * pi)
    frequency <- sample(c(1, 2), size = 1L, prob = c(0.8, 0.2))
    nonlinear_component <- sin(2 * pi * frequency * time + phase)
    bump_center <- stats::runif(1L, 0.25, 0.75)
    bump_width <- stats::runif(1L, 0.08, 0.20)
    bump <- exp(-((time - bump_center)^2) / (2 * bump_width^2))
    nonlinear_component <- nonlinear_component + 0.7 * bump

    ar_coefficient <- stats::runif(1L, 0.40, 0.85)
    ar_component <- as.numeric(stats::arima.sim(
      model = list(ar = ar_coefficient), n = length(time)
    ))

    shape <- 0.45 * linear_component +
      0.75 * nonlinear_component +
      0.35 * ar_component
    return(shape - mean(shape))
  }

  stop("Unknown drift type: ", drift_type, call. = FALSE)
}

#' Generate one canonical metabolomics simulation realization
#'
#' This is the robustness-only, side-effect-free implementation of the design
#' in simulate_metabolomics.R.  The six stochastic components have independent
#' seeds so that every realization is fully specified by its seed ledger.
#'
#' @param seeds Named integer vector/list containing the six required component
#'   seeds listed in `.canonical_required_seed_components`.
#' @param include_artifact_matrices If TRUE, retain the full batch- and
#'   drift-effect matrices in the returned object.  They are useful for deep
#'   diagnostics but are not needed in persisted canonical bundles.
#' @return A named list containing raw and clean intensity matrices, metadata,
#'   feature truth, feature-by-plate truth, design parameters, and seeds.
generate_canonical_simulation <- function(
  seeds,
  include_artifact_matrices = FALSE
) {
  seeds <- unlist(seeds, use.names = TRUE)
  missing_components <- setdiff(.canonical_required_seed_components, names(seeds))
  if (length(missing_components)) {
    stop(
      "Missing canonical simulation seed component(s): ",
      paste(missing_components, collapse = ", "),
      call. = FALSE
    )
  }
  seeds <- stats::setNames(
    as.integer(seeds[.canonical_required_seed_components]),
    .canonical_required_seed_components
  )
  if (anyNA(seeds) || any(seeds <= 0L)) {
    stop("All canonical simulation seeds must be positive integers.", call. = FALSE)
  }

  .canonical_preserve_rng({
    n_features <- 1000L
    n_injections <- 500L
    n_plates <- 5L
    plate_size <- 100L
    qc_interval <- 10L

    prop_batch_affected <- 0.70
    sigma_plate_shift <- 0.60
    prop_batch_strong_within_affected <- 0.10
    sigma_plate_shift_strong <- 1.35

    prop_drift_feature_affected <- 0.60
    prop_drift_active_plate <- 0.55
    drift_type_probabilities <- c(
      linear = 0.35,
      nonlinear = 0.30,
      ar1 = 0.15,
      mixed = 0.20
    )

    stopifnot(n_plates * plate_size == n_injections)

    feature_ids <- paste0("M", seq_len(n_features))
    sample_ids <- sprintf("S%03d", seq_len(n_injections))
    plate_levels <- paste0("P", seq_len(n_plates))

    metadata <- data.frame(
      sample_id = sample_ids,
      plate = rep(plate_levels, each = plate_size),
      order_in_plate = rep(seq_len(plate_size), times = n_plates),
      run_order = seq_len(n_injections),
      stringsAsFactors = FALSE
    )
    metadata$sample_type <- ifelse(
      (metadata$order_in_plate - 1L) %% qc_interval == 0L,
      "control",
      "sample"
    )

    control_indices <- which(metadata$sample_type == "control")
    plate_ids <- match(metadata$plate, plate_levels)
    time_in_plate <- (metadata$order_in_plate - 1) / (plate_size - 1)
    plate_indices <- lapply(seq_len(n_plates), function(plate) {
      which(plate_ids == plate)
    })

    set.seed(seeds[["clean_truth"]])
    clean_log <- matrix(
      stats::rnorm(n_features * n_injections, mean = 0, sd = 0.35),
      nrow = n_features,
      ncol = n_injections,
      dimnames = list(feature_ids, sample_ids)
    )
    if (length(control_indices)) {
      anchor <- control_indices[[1L]]
      clean_log[, control_indices] <- matrix(
        clean_log[, anchor],
        nrow = n_features,
        ncol = length(control_indices)
      )
    }
    clean_intensity <- exp(clean_log)

    set.seed(seeds[["batch_assignment"]])
    batch_affected <- as.logical(stats::rbinom(
      n_features, size = 1, prob = prop_batch_affected
    ))
    batch_strong <- rep(FALSE, n_features)
    if (any(batch_affected)) {
      batch_strong[batch_affected] <- as.logical(stats::rbinom(
        sum(batch_affected),
        size = 1,
        prob = prop_batch_strong_within_affected
      ))
    }
    batch_moderate <- batch_affected & !batch_strong

    batch_shift_by_feature_plate <- matrix(
      0,
      nrow = n_features,
      ncol = n_plates,
      dimnames = list(feature_ids, plate_levels)
    )

    set.seed(seeds[["batch_magnitude"]])
    if (any(batch_moderate)) {
      moderate_shifts <- matrix(
        stats::rnorm(
          sum(batch_moderate) * n_plates,
          mean = 0,
          sd = sigma_plate_shift
        ),
        nrow = sum(batch_moderate),
        ncol = n_plates
      )
      moderate_shifts <- sweep(
        moderate_shifts, 1L, rowMeans(moderate_shifts), FUN = "-"
      )
      batch_shift_by_feature_plate[batch_moderate, ] <- moderate_shifts
    }
    if (any(batch_strong)) {
      strong_shifts <- matrix(
        stats::rnorm(
          sum(batch_strong) * n_plates,
          mean = 0,
          sd = sigma_plate_shift_strong
        ),
        nrow = sum(batch_strong),
        ncol = n_plates
      )
      strong_shifts <- sweep(
        strong_shifts, 1L, rowMeans(strong_shifts), FUN = "-"
      )
      batch_shift_by_feature_plate[batch_strong, ] <- strong_shifts
    }
    batch_shift_by_injection <- batch_shift_by_feature_plate[
      , match(metadata$plate, plate_levels), drop = FALSE
    ]

    set.seed(seeds[["drift_assignment"]])
    drift_feature_affected <- as.logical(stats::rbinom(
      n_features, size = 1, prob = prop_drift_feature_affected
    ))
    drift_active <- matrix(
      FALSE,
      nrow = n_features,
      ncol = n_plates,
      dimnames = list(feature_ids, plate_levels)
    )
    if (any(drift_feature_affected)) {
      active_draw <- matrix(
        as.logical(stats::rbinom(
          sum(drift_feature_affected) * n_plates,
          size = 1,
          prob = prop_drift_active_plate
        )),
        nrow = sum(drift_feature_affected),
        ncol = n_plates
      )
      zero_rows <- which(rowSums(active_draw) == 0L)
      if (length(zero_rows)) {
        forced_plate <- sample(
          seq_len(n_plates), size = length(zero_rows), replace = TRUE
        )
        active_draw[cbind(zero_rows, forced_plate)] <- TRUE
      }
      drift_active[drift_feature_affected, ] <- active_draw
    }

    drift_type <- matrix(
      "none",
      nrow = n_features,
      ncol = n_plates,
      dimnames = list(feature_ids, plate_levels)
    )
    set.seed(seeds[["drift_morphology"]])
    if (sum(drift_active) > 0L) {
      drift_type[drift_active] <- sample(
        names(drift_type_probabilities),
        size = sum(drift_active),
        replace = TRUE,
        prob = as.numeric(drift_type_probabilities)
      )
    }

    drift_amplitude <- matrix(
      0,
      nrow = n_features,
      ncol = n_plates,
      dimnames = list(feature_ids, plate_levels)
    )
    set.seed(seeds[["drift_shape"]])
    if (sum(drift_active) > 0L) {
      drift_amplitude[drift_active] <- stats::runif(
        sum(drift_active), min = 0.20, max = 0.85
      )
    }

    drift_by_injection <- matrix(
      0,
      nrow = n_features,
      ncol = n_injections,
      dimnames = list(feature_ids, sample_ids)
    )
    for (feature in seq_len(n_features)) {
      for (plate in seq_len(n_plates)) {
        if (!drift_active[feature, plate]) next
        indices <- plate_indices[[plate]]
        time <- time_in_plate[indices]
        raw_shape <- .canonical_simulate_plate_drift(
          time, drift_type[feature, plate]
        )
        drift_by_injection[feature, indices] <- .canonical_scale_to_amplitude(
          raw_shape, drift_amplitude[feature, plate]
        )
      }
    }

    raw_log <- clean_log + batch_shift_by_injection + drift_by_injection
    raw_intensity <- exp(raw_log)

    feature_truth <- data.frame(
      metabolite = feature_ids,
      batch_effect_applied = batch_affected,
      batch_effect_strength = ifelse(
        !batch_affected, "none", ifelse(batch_strong, "strong", "moderate")
      ),
      drift_effect_applied_any_plate = rowSums(drift_active) > 0L,
      n_plates_with_drift = as.integer(rowSums(drift_active)),
      drift_type_linear = as.integer(rowSums(drift_type == "linear")),
      drift_type_nonlinear = as.integer(rowSums(drift_type == "nonlinear")),
      drift_type_ar1 = as.integer(rowSums(drift_type == "ar1")),
      drift_type_mixed = as.integer(rowSums(drift_type == "mixed")),
      stringsAsFactors = FALSE
    )
    rownames(feature_truth) <- NULL

    feature_plate_truth <- data.frame(
      metabolite = rep(feature_ids, times = n_plates),
      plate = rep(plate_levels, each = n_features),
      batch_shift_log = as.vector(batch_shift_by_feature_plate),
      drift_active = as.vector(drift_active),
      drift_type = as.vector(drift_type),
      drift_amplitude_log = as.vector(drift_amplitude),
      stringsAsFactors = FALSE
    )

    design_parameters <- list(
      n_features = n_features,
      n_injections = n_injections,
      n_plates = n_plates,
      plate_size = plate_size,
      qc_interval = qc_interval,
      n_controls = length(control_indices),
      n_study_injections = n_injections - length(control_indices),
      clean_log_sd = 0.35,
      prop_batch_affected = prop_batch_affected,
      sigma_plate_shift = sigma_plate_shift,
      prop_batch_strong_within_affected = prop_batch_strong_within_affected,
      sigma_plate_shift_strong = sigma_plate_shift_strong,
      prop_drift_feature_affected = prop_drift_feature_affected,
      prop_drift_active_plate = prop_drift_active_plate,
      drift_type_probabilities = as.list(drift_type_probabilities),
      drift_amplitude_log_min = 0.20,
      drift_amplitude_log_max = 0.85
    )

    stopifnot(
      identical(dim(raw_intensity), c(n_features, n_injections)),
      identical(dim(clean_intensity), c(n_features, n_injections)),
      identical(dimnames(raw_intensity), dimnames(clean_intensity)),
      nrow(metadata) == n_injections,
      nrow(feature_truth) == n_features,
      nrow(feature_plate_truth) == n_features * n_plates,
      sum(metadata$sample_type == "control") == 50L,
      all(is.finite(raw_intensity)),
      all(is.finite(clean_intensity))
    )

    output <- list(
      raw_intensity = raw_intensity,
      clean_ground_truth = clean_intensity,
      metadata = metadata,
      feature_truth = feature_truth,
      feature_plate_truth = feature_plate_truth,
      design_parameters = design_parameters,
      effective_generation_seeds = seeds
    )
    if (isTRUE(include_artifact_matrices)) {
      output$artifact_matrices <- list(
        clean_log = clean_log,
        raw_log = raw_log,
        batch_shift_by_feature_plate = batch_shift_by_feature_plate,
        batch_shift_by_injection = batch_shift_by_injection,
        drift_active = drift_active,
        drift_type = drift_type,
        drift_amplitude = drift_amplitude,
        drift_by_injection = drift_by_injection
      )
    }
    output
  })
}

#' Select the prespecified hidden references for one simulation realization
#'
#' Exactly two of the ten controls on every plate are selected.  This function
#' varies only the assignment and leaves the generated data untouched.
select_canonical_hidden_references <- function(
  metadata,
  seed,
  hidden_per_plate = 2L
) {
  required_columns <- c("sample_id", "plate", "run_order", "sample_type")
  missing_columns <- setdiff(required_columns, names(metadata))
  if (length(missing_columns)) {
    stop(
      "Metadata lacks hidden-reference field(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  seed <- as.integer(seed)
  hidden_per_plate <- as.integer(hidden_per_plate)
  if (length(seed) != 1L || is.na(seed) || seed <= 0L) {
    stop("The hidden-reference seed must be one positive integer.", call. = FALSE)
  }
  if (length(hidden_per_plate) != 1L || is.na(hidden_per_plate) || hidden_per_plate < 1L) {
    stop("hidden_per_plate must be a positive integer.", call. = FALSE)
  }

  .canonical_preserve_rng({
    set.seed(seed)
    controls <- metadata[metadata$sample_type == "control", , drop = FALSE]
    plate_levels <- unique(metadata$plate)
    chosen <- unlist(lapply(plate_levels, function(plate) {
      ids <- controls$sample_id[controls$plate == plate]
      if (length(ids) < hidden_per_plate) {
        stop("Plate ", plate, " has too few controls for the requested split.", call. = FALSE)
      }
      sample(ids, size = hidden_per_plate, replace = FALSE)
    }), use.names = FALSE)

    assignment <- controls[, c("sample_id", "plate", "run_order"), drop = FALSE]
    assignment$is_hidden_reference <- assignment$sample_id %in% chosen
    assignment$reference_role <- ifelse(
      assignment$is_hidden_reference, "hidden_test", "training_control"
    )
    assignment <- assignment[order(assignment$run_order), , drop = FALSE]
    rownames(assignment) <- NULL

    stopifnot(
      sum(assignment$is_hidden_reference) == hidden_per_plate * length(plate_levels),
      all(table(assignment$plate[assignment$is_hidden_reference]) == hidden_per_plate)
    )
    assignment
  })
}
