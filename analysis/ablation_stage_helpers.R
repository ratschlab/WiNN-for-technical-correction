# Derive consecutive ablation-stage changes without mutating cumulative values.

ablation_stage_spec <- function() {
  data.frame(
    stage_index = 1:4,
    stage = c("Outlier_shrinkage", "Selective_drift", "Selective_batch", "PQN"),
    stage_label = c("Outlier shrinkage", "Selective drift", "Selective batch", "PQN"),
    previous_column = c("Raw", "Outlier_shrinkage", "Selective_drift", "Selective_batch"),
    current_column = c("Outlier_shrinkage", "Selective_drift", "Selective_batch", "PQN"),
    stringsAsFactors = FALSE
  )
}

ablation_direction_sign <- function(direction) {
  result <- rep(NA_real_, length(direction))
  result[direction == "higher"] <- 1
  result[direction == "lower"] <- -1
  invalid <- !is.na(direction) & !direction %in% c("higher", "lower", "descriptive")
  if (any(invalid)) {
    stop("Unrecognised metric direction: ", paste(unique(direction[invalid]), collapse = ", "), call. = FALSE)
  }
  result
}

normalize_ablation_cumulative_metric_names <- function(cumulative) {
  cumulative <- as.data.frame(cumulative, stringsAsFactors = FALSE, check.names = FALSE)
  metric_by_table <- c(
    heldout_qc_cv = "heldout_qc_cv_mean",
    gam_deviance = "residual_run_order_gam_mean_deviance",
    batch_r2 = "batch_weighted_pc_r2_categorical",
    feature_repeatability_ratio = "feature_repeatability_ratio_median"
  )
  missing_metric <- is.na(cumulative$metric) | !nzchar(cumulative$metric)
  can_fill <- missing_metric & cumulative$table_id %in% names(metric_by_table)
  cumulative$metric[can_fill] <- unname(metric_by_table[cumulative$table_id[can_fill]])
  if (any(is.na(cumulative$metric) | !nzchar(cumulative$metric))) {
    stop("A cumulative ablation row has no explicit metric identifier.", call. = FALSE)
  }
  cumulative
}

derive_ablation_stage_impacts <- function(cumulative) {
  cumulative <- normalize_ablation_cumulative_metric_names(cumulative)
  id_columns <- c("dataset_key", "dataset", "table_id", "metric", "direction")
  cumulative_columns <- c("Raw", "Outlier_shrinkage", "Selective_drift", "Selective_batch", "PQN")
  required <- c(id_columns, cumulative_columns)
  missing <- setdiff(required, names(cumulative))
  if (length(missing)) {
    stop("Cumulative ablation table is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  duplicate_key <- duplicated(cumulative[id_columns]) | duplicated(cumulative[id_columns], fromLast = TRUE)
  if (any(duplicate_key)) {
    stop("Cumulative ablation rows are duplicated for one or more dataset/metric keys.", call. = FALSE)
  }

  # These copies are deliberately created before any differences are evaluated.
  immutable <- cumulative[id_columns]
  immutable$cum_raw <- as.numeric(cumulative$Raw)
  immutable$cum_outlier <- as.numeric(cumulative$Outlier_shrinkage)
  immutable$cum_drift <- as.numeric(cumulative$Selective_drift)
  immutable$cum_batch <- as.numeric(cumulative$Selective_batch)
  immutable$cum_pqn <- as.numeric(cumulative$PQN)

  spec <- ablation_stage_spec()
  cumulative_map <- c(
    Raw = "cum_raw", Outlier_shrinkage = "cum_outlier",
    Selective_drift = "cum_drift", Selective_batch = "cum_batch", PQN = "cum_pqn"
  )
  result <- do.call(rbind, lapply(seq_len(nrow(spec)), function(index) {
    previous <- immutable[[cumulative_map[[spec$previous_column[index]]]]]
    current <- immutable[[cumulative_map[[spec$current_column[index]]]]]
    raw_change <- current - previous
    direction_sign <- ablation_direction_sign(immutable$direction)
    data.frame(
      immutable,
      stage_index = spec$stage_index[index],
      stage = spec$stage[index],
      stage_label = spec$stage_label[index],
      previous_stage = spec$previous_column[index],
      cumulative_previous = previous,
      cumulative_current = current,
      raw_change = raw_change,
      direction_sign = direction_sign,
      direction_normalized_benefit = direction_sign * raw_change,
      change_definition = "raw_change = cumulative_current - cumulative_previous",
      benefit_definition = "direction_normalized_benefit = direction_sign * raw_change; positive values favor adding the stage",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  rownames(result) <- NULL
  result
}

validate_ablation_stage_impacts <- function(cumulative, impacts, tolerance = 1e-12) {
  cumulative <- as.data.frame(cumulative, stringsAsFactors = FALSE, check.names = FALSE)
  impacts <- as.data.frame(impacts, stringsAsFactors = FALSE, check.names = FALSE)
  id_columns <- c("dataset_key", "dataset", "table_id", "metric", "direction")
  expected <- derive_ablation_stage_impacts(cumulative)
  comparison_columns <- c(
    id_columns, "stage_index", "stage", "stage_label", "previous_stage",
    "cum_raw", "cum_outlier", "cum_drift", "cum_batch", "cum_pqn",
    "cumulative_previous", "cumulative_current", "raw_change", "direction_sign",
    "direction_normalized_benefit"
  )
  missing <- setdiff(comparison_columns, names(impacts))
  if (length(missing)) stop("Impact table is missing: ", paste(missing, collapse = ", "), call. = FALSE)

  order_key <- function(value) {
    do.call(order, c(value[c(id_columns, "stage_index")], list(na.last = TRUE)))
  }
  expected <- expected[order_key(expected), comparison_columns, drop = FALSE]
  observed <- impacts[order_key(impacts), comparison_columns, drop = FALSE]
  if (nrow(expected) != nrow(observed)) stop("Ablation impact row count does not match the cumulative source.", call. = FALSE)

  character_columns <- comparison_columns[!vapply(expected[comparison_columns], is.numeric, logical(1))]
  for (column in character_columns) {
    if (!identical(as.character(expected[[column]]), as.character(observed[[column]]))) {
      stop("Ablation impact identity mismatch in column: ", column, call. = FALSE)
    }
  }
  numeric_columns <- setdiff(comparison_columns, character_columns)
  for (column in numeric_columns) {
    difference <- abs(expected[[column]] - observed[[column]])
    compatible_na <- is.na(expected[[column]]) & is.na(observed[[column]])
    if (any(!compatible_na & (!is.finite(difference) | difference > tolerance))) {
      stop("Ablation impact numerical mismatch in column: ", column, call. = FALSE)
    }
  }

  terminal_columns <- c("cum_outlier", "cum_drift", "cum_batch", "cum_pqn")
  for (index in seq_along(terminal_columns)) {
    rows <- impacts$stage_index <= index
    reconstructed <- impacts$cum_raw[rows] + ave(
      impacts$raw_change[rows],
      interaction(impacts[rows, id_columns], drop = TRUE),
      FUN = cumsum
    )
    target <- impacts[[terminal_columns[index]]][rows]
    stage_rows <- impacts$stage_index[rows] == index
    compatible_na <- is.na(reconstructed[stage_rows]) & is.na(target[stage_rows])
    difference <- abs(reconstructed[stage_rows] - target[stage_rows])
    if (any(!compatible_na & (!is.finite(difference) | difference > tolerance))) {
      stop("A cumulative trajectory does not telescope at stage ", index, ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}
