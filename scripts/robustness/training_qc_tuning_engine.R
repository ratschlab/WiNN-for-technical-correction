# Shared data preparation and training-QC candidate fitting for the benchmark.
# This file contains no endpoint evaluation and accepts only the technical
# metadata and reference identities required by each correction method.

tq_dataset_map <- c(
  mtbls79 = "MTBLS79",
  batchcorr_set1 = "BatchCorr Set 1",
  clinical_fiams = "Clinical FIA-MS (pair-aware)",
  sacurine = "Sacurine",
  waveica = "WaveICA adenocarcinoma"
)

tq_registry_path <- function(repo_root) {
  file.path(repo_root, "analysis", "config", "training_qc_candidate_registry.csv")
}

tq_read_registry <- function(repo_root, dataset, method = NULL) {
  path <- tq_registry_path(repo_root)
  if (!file.exists(path)) stop("Missing candidate registry: ", path, call. = FALSE)
  registry <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("candidate_id", "params", "dataset", "method", "selection_class")
  missing <- setdiff(required, names(registry))
  if (length(missing)) {
    stop("Candidate registry lacks column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  display <- unname(tq_dataset_map[[dataset]])
  if (is.null(display)) stop("Unsupported dataset: ", dataset, call. = FALSE)
  registry <- registry[registry$dataset == display, , drop = FALSE]
  if (!is.null(method)) registry <- registry[registry$method == method, , drop = FALSE]
  method_order <- c("ComBat", "QC-RLSC", "QC-RFSC", "SERRF")
  registry <- registry[
    order(match(registry$method, method_order), registry$candidate_id),
    ,
    drop = FALSE
  ]
  rownames(registry) <- NULL
  if (!nrow(registry)) {
    stop(
      "No candidate definitions for ", display,
      if (is.null(method)) "." else paste0(" / ", method, "."),
      call. = FALSE
    )
  }
  registry
}

tq_prepare_dataset <- function(repo_root, dataset) {
  loaded <- load_winn_ablation_dataset(repo_root, dataset)
  x <- as.matrix(loaded$x)
  storage.mode(x) <- "double"
  meta <- loaded$meta
  correction_unit <- "injection"

  if (dataset == "clinical_fiams") {
    pair <- prepare_clinical_dip_pair_input(x, meta)
    x <- pair$data
    meta <- pair$meta
    correction_unit <- "adjacent_dip_pair"
    meta$is_qc <- meta$role != "clinical"
    meta$is_study <- meta$role == "clinical"
    training_ids <- as.character(meta$sample_id[meta$role == "control"])
    external_ids <- as.character(meta$sample_id[meta$role == "NIST1950"])
    hidden_ids <- character()
  } else {
    hidden_ids <- as.character(loaded$hidden_ids)
    if (!"is_qc" %in% names(meta)) stop(dataset, ": metadata lacks is_qc.", call. = FALSE)
    training_ids <- as.character(
      meta$sample_id[as.logical(meta$is_qc) & !meta$sample_id %in% hidden_ids]
    )
    external_ids <- character()
  }

  meta <- meta[match(colnames(x), meta$sample_id), , drop = FALSE]
  if (anyNA(meta$sample_id) || !identical(colnames(x), as.character(meta$sample_id))) {
    stop(dataset, ": correction matrix and metadata are not aligned.", call. = FALSE)
  }
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop(dataset, ": duplicate matrix identifiers.", call. = FALSE)
  }
  if (any(!is.finite(x)) || any(x < 0)) {
    stop(dataset, ": correction matrix must be finite and nonnegative.", call. = FALSE)
  }
  if (length(intersect(training_ids, hidden_ids)) ||
      length(intersect(training_ids, external_ids)) ||
      length(intersect(hidden_ids, external_ids))) {
    stop(dataset, ": reference roles overlap.", call. = FALSE)
  }
  if (!all(c(training_ids, hidden_ids, external_ids) %in% colnames(x))) {
    stop(dataset, ": one or more reference IDs are absent from the correction matrix.", call. = FALSE)
  }

  list(
    x = x,
    meta = meta,
    training_ids = training_ids,
    hidden_ids = hidden_ids,
    external_ids = external_ids,
    correction_unit = correction_unit
  )
}

tq_method_metadata <- function(meta, fit_qc_ids) {
  fields <- intersect(
    c("sample_id", "batch", "run_order", "within_batch_order"),
    names(meta)
  )
  out <- meta[, fields, drop = FALSE]
  out$class <- "Sample"
  out$sample_type <- "Study"
  out$is_qc <- FALSE
  out$is_study <- TRUE
  out$role <- "ordinary"
  selected <- out$sample_id %in% fit_qc_ids
  out$class[selected] <- "QC"
  out$sample_type[selected] <- "QC"
  out$is_qc[selected] <- TRUE
  out$is_study[selected] <- FALSE
  out$role[selected] <- "training_reference"
  out
}

tq_technical_metadata <- function(meta) {
  fields <- intersect(
    c("sample_id", "batch", "run_order", "within_batch_order"),
    names(meta)
  )
  meta[, fields, drop = FALSE]
}

tq_candidate_parameter_fields <- c(
  "par_prior", "span", "ntree", "coCV", "Frule", "use_injection",
  "jitter_eps", "test", "spline_method", "acorr_fdr", "anova_fdr",
  "normalization", "scale_by_batch", "pelt_penalty"
)

tq_candidate_parameters <- function(candidate_row) {
  fields <- intersect(tq_candidate_parameter_fields, names(candidate_row))
  result <- as.list(candidate_row[1L, fields, drop = FALSE])
  if (identical(as.character(candidate_row$method[1L]), "ComBat") &&
      (is.null(result$par_prior) || is.na(result$par_prior))) {
    result$par_prior <- TRUE
  }
  result
}

tq_qc_rlsc_min_qcs_per_batch <- function(dataset) {
  # The MTBLS79 holdout leaves three training QCs in six of eight batches;
  # the direct subtractive qcrlscR fit supports this design. Other datasets
  # retain the four-QC safeguard.
  if (identical(dataset, "mtbls79")) 3L else 4L
}

tq_run_candidate <- function(
  method, x, fit_qc_ids, meta_method, meta_technical, candidate_row,
  dataset = NULL
) {
  parameters <- tq_candidate_parameters(candidate_row)
  switch(
    method,
    ComBat = run_combat(
      x, meta_technical, par_prior = as.logical(parameters$par_prior)
    ),
    `QC-RLSC` = run_qc_rlsc_id_safe(
      x, fit_qc_ids, meta_method, span = as.numeric(parameters$span),
      degree = 1L, shift_batches = TRUE,
      min_qcs_per_batch = tq_qc_rlsc_min_qcs_per_batch(dataset)
    ),
    `QC-RFSC` = run_qc_rfsc_with_controls(
      x, fit_qc_ids, meta_method, ntree = as.integer(parameters$ntree),
      coCV = as.numeric(parameters$coCV), Frule = as.numeric(parameters$Frule)
    ),
    SERRF = run_serrf_all_corrected(
      x, fit_qc_ids, meta_method, jitter_eps = as.numeric(parameters$jitter_eps)
    ),
    stop("Unsupported candidate method: ", method, call. = FALSE)
  )
}
