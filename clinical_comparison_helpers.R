# Generic adjacent-pair handling for the restricted FIA-MS benchmark. This file
# defines the correction-unit transform only; it contains no cohort records,
# identifiers, feature names, or private input paths.

align_meta_to_mat <- function(x_mat, meta_df) {
  required <- c(
    "sample_id", "rep_unit", "dip", "run_order", "batch", "role",
    "sample", "replicate"
  )
  missing <- setdiff(required, names(meta_df))
  if (length(missing)) {
    stop("Pair metadata lack column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  metadata <- meta_df[match(colnames(x_mat), meta_df$sample_id), , drop = FALSE]
  if (anyNA(metadata$sample_id) || !identical(colnames(x_mat), as.character(metadata$sample_id))) {
    stop("Matrix columns do not align exactly with pair metadata.", call. = FALSE)
  }
  metadata
}

prepare_clinical_dip_pair_input <- function(x_mat, meta_df) {
  x_mat <- as.matrix(x_mat)
  storage.mode(x_mat) <- "double"
  metadata <- align_meta_to_mat(x_mat, meta_df)
  if (any(!is.finite(x_mat)) || any(x_mat < 0)) {
    stop("Pair correction requires finite, nonnegative intensities.", call. = FALSE)
  }

  metadata <- metadata[order(metadata$run_order), , drop = FALSE]
  unit_ids <- unique(as.character(metadata$rep_unit))
  split_index <- lapply(unit_ids, function(unit_id) which(metadata$rep_unit == unit_id))
  names(split_index) <- unit_ids
  valid <- vapply(split_index, function(index) {
    length(index) == 2L &&
      identical(sort(as.integer(metadata$dip[index])), c(1L, 2L)) &&
      diff(sort(as.numeric(metadata$run_order[index]))) == 1 &&
      length(unique(metadata$batch[index])) == 1L &&
      length(unique(metadata$role[index])) == 1L &&
      length(unique(metadata$sample[index])) == 1L &&
      length(unique(metadata$replicate[index])) == 1L
  }, logical(1))
  if (!all(valid)) {
    stop(sum(!valid), " correction unit(s) are not adjacent two-injection pairs.", call. = FALSE)
  }

  unit_log <- vapply(
    split_index,
    function(index) rowMeans(log1p(x_mat[, metadata$sample_id[index], drop = FALSE])),
    numeric(nrow(x_mat))
  )
  if (is.null(dim(unit_log))) unit_log <- matrix(unit_log, nrow = nrow(x_mat), ncol = 1L)
  rownames(unit_log) <- rownames(x_mat)
  colnames(unit_log) <- unit_ids
  unit_matrix <- expm1(unit_log)

  first_index <- vapply(split_index, `[[`, integer(1), 1L)
  unit_metadata <- metadata[first_index, , drop = FALSE]
  unit_metadata$sample_id <- unit_ids
  unit_metadata$dip <- NA_integer_
  unit_metadata$run_order <- vapply(
    split_index, function(index) mean(metadata$run_order[index]), numeric(1)
  )
  unit_metadata$n_dips <- 2L
  unit_metadata$within_batch_order <- ave(
    unit_metadata$run_order, unit_metadata$batch, FUN = rank
  )
  unit_metadata <- unit_metadata[match(unit_ids, unit_metadata$sample_id), , drop = FALSE]

  mapping <- do.call(rbind, lapply(seq_along(split_index), function(index) {
    injection_index <- split_index[[index]]
    data.frame(
      sample_id = metadata$sample_id[injection_index],
      correction_unit_id = unit_ids[index],
      dip = metadata$dip[injection_index],
      original_run_order = metadata$run_order[injection_index],
      correction_run_order = mean(metadata$run_order[injection_index]),
      batch = metadata$batch[injection_index],
      role = metadata$role[injection_index],
      stringsAsFactors = FALSE
    )
  }))
  mapping <- mapping[match(colnames(x_mat), mapping$sample_id), , drop = FALSE]

  list(
    data = unit_matrix, log_data = unit_log, meta = unit_metadata,
    mapping = mapping,
    control_ids = unit_metadata$sample_id[unit_metadata$role == "control"],
    nist_ids = unit_metadata$sample_id[unit_metadata$role == "NIST1950"],
    clinical_ids = unit_metadata$sample_id[unit_metadata$role == "clinical"]
  )
}

expand_clinical_dip_pair_correction <- function(corrected_unit_mat, pair_input, original_mat) {
  unit_raw <- pair_input$data
  corrected_unit_mat <- as.matrix(corrected_unit_mat)
  if (!all(rownames(unit_raw) %in% rownames(corrected_unit_mat)) ||
      !all(colnames(unit_raw) %in% colnames(corrected_unit_mat))) {
    stop("Corrected matrix lost required pair units.", call. = FALSE)
  }
  corrected_unit_mat <- corrected_unit_mat[
    rownames(unit_raw), colnames(unit_raw), drop = FALSE
  ]
  if (any(!is.finite(corrected_unit_mat)) || any(corrected_unit_mat < 0)) {
    stop("Corrected pair matrix contains invalid intensities.", call. = FALSE)
  }

  log_adjustment <- log1p(corrected_unit_mat) - log1p(unit_raw)
  unit_index <- match(pair_input$mapping$correction_unit_id, colnames(unit_raw))
  original <- original_mat[
    rownames(unit_raw), pair_input$mapping$sample_id, drop = FALSE
  ]
  expanded <- expm1(log1p(original) + log_adjustment[, unit_index, drop = FALSE])
  expanded[expanded < 0] <- 0
  dimnames(expanded) <- dimnames(original_mat)
  expanded
}
