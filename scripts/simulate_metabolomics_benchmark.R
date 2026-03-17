# Manuscript-aligned simulation benchmark generator for the WiNN paper.
# This script implements the synthetic data design described in the
# "Simulation benchmark" section and Supplementary Methods S6.

scale_to_amp <- function(x, amp) {
  x <- x - mean(x)
  max_abs <- max(abs(x))

  if (!is.finite(max_abs) || max_abs <= 1e-12) {
    return(rep(0, length(x)))
  }

  amp * (x / max_abs)
}

simulate_plate_drift_shape <- function(tt, drift_type) {
  if (drift_type == "linear") {
    slope <- stats::rnorm(1, mean = 0, sd = 1)
    shape <- slope * (tt - mean(tt))
    return(shape - mean(shape))
  }

  if (drift_type == "nonlinear") {
    phase <- stats::runif(1, 0, 2 * pi)
    freq <- sample(c(1, 2), size = 1, prob = c(0.75, 0.25))
    sinusoid <- sin(2 * pi * freq * tt + phase)
    bump_center <- stats::runif(1, 0.2, 0.8)
    bump_width <- stats::runif(1, 0.08, 0.22)
    bump <- exp(-((tt - bump_center) ^ 2) / (2 * bump_width ^ 2))
    shape <- sinusoid + stats::runif(1, 0.5, 0.9) * bump
    return(shape - mean(shape))
  }

  if (drift_type == "ar1") {
    ar_coef <- stats::runif(1, 0.45, 0.90)
    shape <- as.numeric(stats::arima.sim(model = list(ar = ar_coef), n = length(tt)))
    return(shape - mean(shape))
  }

  if (drift_type == "mixed") {
    slope <- stats::rnorm(1, mean = 0, sd = 0.8)
    linear_component <- slope * (tt - mean(tt))

    phase <- stats::runif(1, 0, 2 * pi)
    freq <- sample(c(1, 2), size = 1, prob = c(0.8, 0.2))
    nonlinear_component <- sin(2 * pi * freq * tt + phase)
    bump_center <- stats::runif(1, 0.25, 0.75)
    bump_width <- stats::runif(1, 0.08, 0.20)
    bump <- exp(-((tt - bump_center) ^ 2) / (2 * bump_width ^ 2))
    nonlinear_component <- nonlinear_component + 0.7 * bump

    ar_coef <- stats::runif(1, 0.40, 0.85)
    ar_component <- as.numeric(stats::arima.sim(model = list(ar = ar_coef), n = length(tt)))

    shape <- 0.45 * linear_component + 0.75 * nonlinear_component + 0.35 * ar_component
    return(shape - mean(shape))
  }

  stop("Unknown drift type: ", drift_type)
}

simulate_metabolomics_benchmark <- function(
  output_dir = file.path("data", "simulated"),
  verbose = interactive()
) {
  required_pkgs <- c("dplyr", "readr", "tibble")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0L) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }

  seeds <- list(
    gt = 20260311L,
    batch_assignment = 20260312L,
    batch_shift = 20260313L,
    drift_assignment = 20260314L,
    drift_type = 20260315L,
    drift_shape = 20260316L
  )

  n_metabolites <- 1000L
  n_injections <- 500L
  n_plates <- 5L
  plate_size <- 100L
  qc_interval <- 10L

  stopifnot(n_plates * plate_size == n_injections)

  prop_batch_affected <- 0.70
  sigma_plate_shift <- 0.60
  prop_batch_strong_within_affected <- 0.10
  sigma_plate_shift_strong <- 1.35

  prop_drift_met_affected <- 0.60
  prop_drift_active_plate <- 0.55
  drift_type_probs <- c(
    linear = 0.35,
    nonlinear = 0.30,
    ar1 = 0.15,
    mixed = 0.20
  )

  meta <- tibble::tibble(
    sample_id = sprintf("S%03d", seq_len(n_injections)),
    plate = rep(paste0("P", seq_len(n_plates)), each = plate_size),
    order_in_plate = rep(seq_len(plate_size), times = n_plates),
    run_order = seq_len(n_injections)
  ) |>
    dplyr::mutate(
      sample_type = dplyr::if_else((order_in_plate - 1L) %% qc_interval == 0L, "control", "sample")
    )

  control_idx <- which(meta$sample_type == "control")
  plate_levels <- unique(meta$plate)
  plate_id <- as.integer(factor(meta$plate, levels = plate_levels))
  t_in_plate <- (meta$order_in_plate - 1) / (plate_size - 1)
  plate_indices <- lapply(seq_len(n_plates), function(idx) which(plate_id == idx))

  set.seed(seeds$gt)
  sd_gt <- 0.35
  ground_truth_log <- matrix(
    stats::rnorm(n_metabolites * n_injections, mean = 0, sd = sd_gt),
    nrow = n_metabolites,
    ncol = n_injections,
    dimnames = list(paste0("M", seq_len(n_metabolites)), meta$sample_id)
  )

  if (length(control_idx) > 0L) {
    anchor <- control_idx[1]
    ground_truth_log[, control_idx] <- matrix(
      ground_truth_log[, anchor],
      nrow = n_metabolites,
      ncol = length(control_idx)
    )
  }

  ground_truth <- exp(ground_truth_log)

  set.seed(seeds$batch_assignment)
  batch_affected <- as.logical(stats::rbinom(n_metabolites, size = 1, prob = prop_batch_affected))
  batch_strong <- rep(FALSE, n_metabolites)
  if (any(batch_affected)) {
    batch_strong[batch_affected] <- as.logical(stats::rbinom(
      sum(batch_affected),
      size = 1,
      prob = prop_batch_strong_within_affected
    ))
  }
  batch_moderate <- batch_affected & !batch_strong

  shift_by_metabolite_plate <- matrix(
    0,
    nrow = n_metabolites,
    ncol = n_plates,
    dimnames = list(rownames(ground_truth_log), plate_levels)
  )

  set.seed(seeds$batch_shift)
  if (any(batch_moderate)) {
    raw_shifts <- matrix(
      stats::rnorm(sum(batch_moderate) * n_plates, mean = 0, sd = sigma_plate_shift),
      nrow = sum(batch_moderate),
      ncol = n_plates
    )
    raw_shifts <- sweep(raw_shifts, 1, rowMeans(raw_shifts), FUN = "-")
    shift_by_metabolite_plate[batch_moderate, ] <- raw_shifts
  }

  if (any(batch_strong)) {
    strong_shifts <- matrix(
      stats::rnorm(sum(batch_strong) * n_plates, mean = 0, sd = sigma_plate_shift_strong),
      nrow = sum(batch_strong),
      ncol = n_plates
    )
    strong_shifts <- sweep(strong_shifts, 1, rowMeans(strong_shifts), FUN = "-")
    shift_by_metabolite_plate[batch_strong, ] <- strong_shifts
  }

  plate_shift_by_sample <- shift_by_metabolite_plate[, match(meta$plate, plate_levels), drop = FALSE]

  set.seed(seeds$drift_assignment)
  drift_met_affected <- as.logical(stats::rbinom(
    n_metabolites,
    size = 1,
    prob = prop_drift_met_affected
  ))

  drift_active <- matrix(
    FALSE,
    nrow = n_metabolites,
    ncol = n_plates,
    dimnames = list(rownames(ground_truth_log), plate_levels)
  )

  if (any(drift_met_affected)) {
    active_draw <- matrix(
      as.logical(stats::rbinom(
        sum(drift_met_affected) * n_plates,
        size = 1,
        prob = prop_drift_active_plate
      )),
      nrow = sum(drift_met_affected),
      ncol = n_plates
    )

    zero_rows <- which(rowSums(active_draw) == 0L)
    if (length(zero_rows) > 0L) {
      forced_plate <- sample(seq_len(n_plates), size = length(zero_rows), replace = TRUE)
      active_draw[cbind(zero_rows, forced_plate)] <- TRUE
    }

    drift_active[drift_met_affected, ] <- active_draw
  }

  drift_type_map <- matrix(
    "none",
    nrow = n_metabolites,
    ncol = n_plates,
    dimnames = list(rownames(ground_truth_log), plate_levels)
  )

  set.seed(seeds$drift_type)
  if (sum(drift_active) > 0L) {
    drift_type_map[drift_active] <- sample(
      names(drift_type_probs),
      size = sum(drift_active),
      replace = TRUE,
      prob = as.numeric(drift_type_probs)
    )
  }

  drift_amplitude <- matrix(
    0,
    nrow = n_metabolites,
    ncol = n_plates,
    dimnames = list(rownames(ground_truth_log), plate_levels)
  )

  set.seed(seeds$drift_shape)
  if (sum(drift_active) > 0L) {
    drift_amplitude[drift_active] <- stats::runif(sum(drift_active), min = 0.20, max = 0.85)
  }

  drift_matrix <- matrix(
    0,
    nrow = n_metabolites,
    ncol = n_injections,
    dimnames = list(rownames(ground_truth_log), meta$sample_id)
  )

  for (metabolite_idx in seq_len(n_metabolites)) {
    for (plate_idx in seq_len(n_plates)) {
      if (!drift_active[metabolite_idx, plate_idx]) {
        next
      }

      sample_idx <- plate_indices[[plate_idx]]
      tt <- t_in_plate[sample_idx]
      raw_shape <- simulate_plate_drift_shape(tt, drift_type_map[metabolite_idx, plate_idx])
      drift_matrix[metabolite_idx, sample_idx] <- scale_to_amp(
        raw_shape,
        amp = drift_amplitude[metabolite_idx, plate_idx]
      )
    }
  }

  shifted_profiles_log <- ground_truth_log + plate_shift_by_sample + drift_matrix
  shifted_profiles <- exp(shifted_profiles_log)

  design_table <- tibble::tibble(
    metabolite = rownames(ground_truth_log),
    batch_effect_applied = batch_affected,
    batch_effect_strength = dplyr::case_when(
      !batch_affected ~ "none",
      batch_strong ~ "strong",
      TRUE ~ "moderate"
    ),
    drift_effect_applied_any_plate = rowSums(drift_active) > 0,
    n_plates_with_drift = rowSums(drift_active),
    drift_type_linear = rowSums(drift_type_map == "linear"),
    drift_type_nonlinear = rowSums(drift_type_map == "nonlinear"),
    drift_type_ar1 = rowSums(drift_type_map == "ar1"),
    drift_type_mixed = rowSums(drift_type_map == "mixed")
  )

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ground_truth_path <- file.path(output_dir, "ground_truth.tsv")
  shifted_profiles_path <- file.path(output_dir, "shifted_profiles.tsv")
  metadata_path <- file.path(output_dir, "sample_metadata.tsv")
  design_table_path <- file.path(output_dir, "simulation_design_table.tsv")

  readr::write_tsv(
    tibble::as_tibble(ground_truth, rownames = "metabolite"),
    ground_truth_path
  )
  readr::write_tsv(
    tibble::as_tibble(shifted_profiles, rownames = "metabolite"),
    shifted_profiles_path
  )
  readr::write_tsv(meta, metadata_path)
  readr::write_tsv(design_table, design_table_path)

  strong_subset_pct <- if (any(batch_affected)) {
    100 * mean(batch_strong[batch_affected])
  } else {
    0
  }

  summary <- list(
    n_metabolites = n_metabolites,
    n_injections = n_injections,
    n_plates = n_plates,
    controls = length(control_idx),
    study_samples = n_injections - length(control_idx),
    qc_interval = qc_interval,
    batch_affected_pct = 100 * mean(batch_affected),
    strong_batch_subset_pct = strong_subset_pct,
    drift_active_pair_pct = 100 * mean(as.numeric(drift_active))
  )

  if (isTRUE(verbose)) {
    message("Saved simulation benchmark files to: ", normalizePath(output_dir, winslash = "/"))
    message("Controls inserted every ", qc_interval, " injections per plate.")
    message(sprintf("Batch effects applied to %.1f%% of metabolites.", summary$batch_affected_pct))
    message(sprintf(
      "Strong batch shifts assigned to %.1f%% of batch-affected metabolites.",
      summary$strong_batch_subset_pct
    ))
    message(sprintf(
      "Drift-active metabolite-plate pairs: %.1f%% of all pairs.",
      summary$drift_active_pair_pct
    ))
  }

  invisible(list(
    ground_truth = ground_truth,
    shifted_profiles = shifted_profiles,
    metadata = meta,
    design_table = design_table,
    output_paths = c(
      ground_truth = ground_truth_path,
      shifted_profiles = shifted_profiles_path,
      sample_metadata = metadata_path,
      simulation_design_table = design_table_path
    ),
    summary = summary
  ))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output_dir <- if (length(args) >= 1L) args[[1]] else file.path("data", "simulated")
  simulate_metabolomics_benchmark(output_dir = output_dir, verbose = TRUE)
}
