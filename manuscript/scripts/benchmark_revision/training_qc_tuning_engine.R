# Shared execution primitives for the frozen training-QC-only tuning
# sensitivity analysis.  This file deliberately contains no top-level data
# loading or endpoint evaluation: candidate selection scripts can source it
# without making held-out/external outcomes available to the selection phase.

tq_dataset_map <- c(
  mtbls79 = "MTBLS79",
  batchcorr_set1 = "BatchCorr Set 1",
  clinical_fiams = "Clinical FIA-MS (pair-aware)",
  sacurine = "Sacurine",
  waveica = "WaveICA adenocarcinoma"
)

tq_method_order <- c(
  "Raw", "ComBat", "QC-RLSC", "QC-RFSC", "TIGER", "SERRF",
  "WINN auto (QC)", "WINN auto-batch (QC)", "WINN default (no QC)"
)

tq_method_slug <- function(method) {
  value <- gsub("[^a-z0-9]+", "_", tolower(method))
  gsub("^_+|_+$", "", value)
}

tq_sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

tq_sha256_object <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

tq_write_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
}

tq_relative_path <- function(path, repo_root) {
  normalized <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(normalizePath(repo_root, mustWork = TRUE), .Platform$file.sep)
  if (startsWith(normalized, prefix)) {
    substring(normalized, nchar(prefix) + 1L)
  } else {
    normalized
  }
}

tq_atomic_save_rds <- function(value, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(basename(path), ".tmp_"), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(value, temporary, compress = compress)
  if (!file.rename(temporary, path)) stop("Could not atomically install cache: ", path)
  invisible(path)
}

# A bundle manifest deliberately excludes only itself: hashing the manifest
# from within the manifest would be circular.  Every other regular file below
# the bundle root, including cached matrices and human-readable logs, is listed
# with an exact byte count and SHA-256 digest.
tq_artifact_manifest_name <- "artifact_manifest.csv"

tq_list_bundle_files <- function(directory, manifest_name = tq_artifact_manifest_name) {
  if (!dir.exists(directory)) return(character())
  paths <- list.files(
    directory, recursive = TRUE, full.names = TRUE, all.files = TRUE,
    no.. = TRUE, include.dirs = FALSE
  )
  info <- file.info(paths)
  paths <- paths[!is.na(info$isdir) & !info$isdir]
  manifest_path <- normalizePath(
    file.path(directory, manifest_name), mustWork = FALSE
  )
  paths[normalizePath(paths, mustWork = FALSE) != manifest_path]
}

tq_build_artifact_manifest <- function(directory, manifest_name = tq_artifact_manifest_name) {
  directory <- normalizePath(directory, mustWork = TRUE)
  paths <- tq_list_bundle_files(directory, manifest_name)
  relative <- if (length(paths)) {
    substring(normalizePath(paths, mustWork = TRUE), nchar(directory) + 2L)
  } else {
    character()
  }
  info <- file.info(paths)
  data.frame(
    path = relative,
    size_bytes = as.numeric(info$size),
    sha256 = vapply(paths, tq_sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
}

tq_write_artifact_manifest <- function(directory, manifest_name = tq_artifact_manifest_name) {
  manifest <- tq_build_artifact_manifest(directory, manifest_name)
  path <- file.path(directory, manifest_name)
  temporary <- tempfile(paste0(manifest_name, ".tmp_"), tmpdir = directory)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  write.csv(manifest, temporary, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(temporary, path)) stop("Could not atomically install artifact manifest: ", path)
  invisible(manifest)
}

tq_validate_artifact_manifest <- function(directory, manifest_name = tq_artifact_manifest_name) {
  path <- file.path(directory, manifest_name)
  fail <- function(reason, manifest = data.frame()) {
    list(valid = FALSE, reason = reason, manifest = manifest)
  }
  if (!file.exists(path)) return(fail(paste0("missing ", manifest_name)))
  manifest <- tryCatch(
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) e
  )
  if (inherits(manifest, "error")) return(fail(conditionMessage(manifest)))
  required <- c("path", "size_bytes", "sha256")
  if (!all(required %in% names(manifest))) return(fail("artifact manifest lacks required columns", manifest))
  if (anyDuplicated(manifest$path)) return(fail("artifact manifest contains duplicate paths", manifest))
  unsafe <- !nzchar(manifest$path) | grepl("(^/|(^|/)\\.\\.(/|$))", manifest$path)
  if (any(unsafe)) return(fail("artifact manifest contains an unsafe path", manifest))

  actual_paths <- tq_list_bundle_files(directory, manifest_name)
  root <- normalizePath(directory, mustWork = TRUE)
  actual_relative <- if (length(actual_paths)) {
    substring(normalizePath(actual_paths, mustWork = TRUE), nchar(root) + 2L)
  } else {
    character()
  }
  if (!identical(sort(as.character(manifest$path)), sort(actual_relative))) {
    return(fail("artifact manifest is incomplete or lists absent files", manifest))
  }
  resolved <- file.path(directory, manifest$path)
  info <- file.info(resolved)
  size_ok <- !anyNA(info$size) && identical(
    as.numeric(info$size), as.numeric(manifest$size_bytes)
  )
  hash_ok <- all(vapply(seq_along(resolved), function(index) {
    identical(tq_sha256_file(resolved[index]), as.character(manifest$sha256[index]))
  }, logical(1)))
  if (!size_ok) return(fail("one or more artifact byte counts changed", manifest))
  if (!hash_ok) return(fail("one or more artifact SHA-256 digests changed", manifest))
  list(valid = TRUE, reason = "", manifest = manifest)
}

tq_capture_call <- function(fn) {
  warnings <- character()
  messages <- character()
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(
      fn(),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) e
  )
  list(
    value = value,
    elapsed_sec = proc.time()[["elapsed"]] - started,
    warnings = unique(warnings),
    messages = unique(messages),
    error = if (inherits(value, "error")) conditionMessage(value) else ""
  )
}

tq_parse_bool <- function(value) {
  if (is.logical(value)) return(value)
  normalized <- tolower(trimws(as.character(value)))
  result <- rep(NA, length(normalized))
  result[normalized %in% c("true", "t", "1")] <- TRUE
  result[normalized %in% c("false", "f", "0")] <- FALSE
  result
}

tq_registry_path <- function(repo_root) {
  file.path(
    repo_root, "results", "robustness", "07_training_qc_tuning_sensitivity",
    "candidate_registry.csv"
  )
}

tq_fold_manifest_path <- function(repo_root) {
  file.path(
    repo_root, "results", "robustness", "07_training_qc_tuning_sensitivity",
    "internal_qc_fold_manifest.csv"
  )
}

tq_output_root <- function(repo_root) {
  file.path(
    repo_root, "results", "robustness", "07_training_qc_tuning_sensitivity"
  )
}

tq_read_registry <- function(repo_root, dataset, method = NULL) {
  path <- tq_registry_path(repo_root)
  if (!file.exists(path)) stop("Missing frozen candidate registry: ", path)
  registry <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  display <- unname(tq_dataset_map[[dataset]])
  if (is.null(display)) stop("Unsupported tuning-sensitivity dataset: ", dataset)
  registry <- registry[registry$dataset == display, , drop = FALSE]
  if (!is.null(method)) registry <- registry[registry$method == method, , drop = FALSE]
  registry <- registry[order(match(registry$method, tq_method_order), registry$candidate_id), , drop = FALSE]
  rownames(registry) <- NULL
  if (!nrow(registry)) stop("No frozen registry rows for ", display, " / ", method %||% "all methods", ".")
  registry
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

tq_tuning_feature_spec <- function(dataset) {
  switch(
    dataset,
    mtbls79 = list(seed = 20260306L, n_features = 220L),
    batchcorr_set1 = list(seed = 20260722L, n_features = 120L),
    clinical_fiams = list(seed = 20260312L, n_features = 250L),
    sacurine = list(seed = 20260812L, n_features = Inf),
    waveica = list(seed = 20260812L, n_features = 220L),
    stop("No frozen tuning-feature specification for ", dataset, ".")
  )
}

tq_tuning_features <- function(x, dataset) {
  spec <- tq_tuning_feature_spec(dataset)
  n_features <- if (is.finite(spec$n_features)) {
    min(nrow(x), as.integer(spec$n_features))
  } else {
    nrow(x)
  }
  if (n_features == nrow(x)) return(rownames(x))
  set.seed(spec$seed)
  sample(rownames(x), n_features, replace = FALSE)
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
    if (!"is_qc" %in% names(meta)) stop(dataset, ": metadata lacks is_qc.")
    training_ids <- as.character(meta$sample_id[as.logical(meta$is_qc) & !meta$sample_id %in% hidden_ids])
    external_ids <- character()
  }

  meta <- meta[match(colnames(x), meta$sample_id), , drop = FALSE]
  if (anyNA(meta$sample_id) || !identical(colnames(x), as.character(meta$sample_id))) {
    stop(dataset, ": correction matrix and metadata are not aligned.")
  }
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop(dataset, ": duplicate matrix identifiers.")
  }
  if (any(!is.finite(x)) || any(x < 0)) {
    stop(dataset, ": correction matrix must be finite and nonnegative.")
  }
  if (length(intersect(training_ids, hidden_ids)) ||
      length(intersect(training_ids, external_ids)) ||
      length(intersect(hidden_ids, external_ids))) {
    stop(dataset, ": reference roles overlap.")
  }
  if (!all(c(training_ids, hidden_ids, external_ids) %in% colnames(x))) {
    stop(dataset, ": one or more reference IDs are absent from the correction matrix.")
  }

  list(
    x = x,
    meta = meta,
    training_ids = training_ids,
    hidden_ids = hidden_ids,
    external_ids = external_ids,
    correction_unit = correction_unit,
    archived_fixed_path = loaded$fixed_path
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

tq_fold_definitions <- function(repo_root, dataset, training_ids, meta) {
  display <- unname(tq_dataset_map[[dataset]])
  feasibility <- read.csv(
    file.path(tq_output_root(repo_root), "qc_fold_feasibility.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  row <- feasibility[feasibility$dataset == display, , drop = FALSE]
  if (nrow(row) != 1L) stop("Missing unique fold-feasibility row for ", display, ".")
  mode <- as.character(row$tuning_mode)

  if (identical(mode, "in_sample_training_qc_sensitivity")) {
    return(list(
      mode = mode,
      rows = list(data.frame(
        repeat_id = 0L,
        repeat_seed = NA_integer_,
        validation_fold = 0L,
        fold_id = "full_training_qcs_in_sample",
        fit_qc_ids = I(list(sort(training_ids))),
        validation_qc_ids = I(list(sort(training_ids))),
        stringsAsFactors = FALSE
      )),
      ledger = data.frame(
        repeat_id = 0L, repeat_seed = NA_integer_, batch = "all",
        qc_id = sort(training_ids), run_order = meta$run_order[match(sort(training_ids), meta$sample_id)],
        protected_training_endpoint = TRUE, validation_fold = 0L,
        dataset = display, stringsAsFactors = FALSE
      )
    ))
  }

  manifest_path <- tq_fold_manifest_path(repo_root)
  manifest <- read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
  ledger <- manifest[manifest$dataset == display, , drop = FALSE]
  if (!nrow(ledger)) stop("No internal QC folds for ", display, ".")
  manifest_ids <- sort(unique(as.character(ledger$qc_id)))
  if (!identical(manifest_ids, sort(as.character(training_ids)))) {
    stop(display, ": frozen fold manifest does not exactly match the training-QC panel.")
  }

  combinations <- unique(ledger[, c("repeat_id", "repeat_seed", "validation_fold"), drop = FALSE])
  combinations <- combinations[!is.na(combinations$validation_fold), , drop = FALSE]
  combinations <- combinations[order(combinations$repeat_id, combinations$validation_fold), , drop = FALSE]
  rows <- lapply(seq_len(nrow(combinations)), function(index) {
    key <- combinations[index, , drop = FALSE]
    in_repeat <- ledger$repeat_id == key$repeat_id
    validation <- as.character(ledger$qc_id[
      in_repeat & ledger$validation_fold == key$validation_fold & !is.na(ledger$validation_fold)
    ])
    fit <- setdiff(as.character(training_ids), validation)
    data.frame(
      repeat_id = as.integer(key$repeat_id),
      repeat_seed = as.integer(key$repeat_seed),
      validation_fold = as.integer(key$validation_fold),
      fold_id = sprintf("repeat%02d_fold%02d", key$repeat_id, key$validation_fold),
      fit_qc_ids = I(list(sort(fit))),
      validation_qc_ids = I(list(sort(validation))),
      stringsAsFactors = FALSE
    )
  })
  list(mode = mode, rows = rows, ledger = ledger)
}

tq_candidate_seed <- function(dataset, method, candidate_id, fold_id, phase = "candidate") {
  key <- paste(dataset, method, as.integer(candidate_id), fold_id, phase, sep = "::")
  # A seven-hex-digit prefix is deterministic and safely inside R's
  # signed-integer range.  The nonzero offset avoids the sentinel seed 0.
  hex <- substr(digest::digest(key, algo = "sha256", serialize = FALSE), 1L, 7L)
  as.integer(strtoi(hex, base = 16L) + 1L)
}

tq_candidate_parameter_fields <- c(
  "par_prior", "span", "ntree", "coCV", "Frule", "use_injection",
  "jitter_eps", "test", "spline_method", "acorr_fdr", "anova_fdr",
  "normalization", "scale_by_batch", "pelt_penalty"
)

tq_candidate_parameters <- function(candidate_row) {
  fields <- intersect(tq_candidate_parameter_fields, names(candidate_row))
  result <- as.list(candidate_row[1L, fields, drop = FALSE])
  # Fixed ComBat rows in Sacurine/WaveICA encode the published/default setting
  # in `params`, rather than duplicating it in the par_prior column.
  if (identical(as.character(candidate_row$method[1L]), "ComBat") &&
      (is.null(result$par_prior) || is.na(result$par_prior))) {
    result$par_prior <- TRUE
  }
  result
}

tq_qc_rlsc_min_qcs_per_batch <- function(dataset) {
  # The final MTBLS79 12-QC holdout archive was produced with three training
  # QCs in six of eight batches (the qcrlscR direct subtractive fit supports
  # that design). Four is retained everywhere else and remains the threshold
  # that makes an additional internal MTBLS79 QC fold infeasible.  These are
  # distinct constraints: do not turn the fold-feasibility safeguard into a
  # full-data method-fit prohibition.
  if (identical(dataset, "mtbls79")) 3L else 4L
}

tq_run_candidate <- function(method, x, fit_qc_ids, meta_method, meta_technical, candidate_row,
                             dataset = NULL) {
  p <- tq_candidate_parameters(candidate_row)
  switch(
    method,
    Raw = x,
    ComBat = run_combat(x, meta_technical, par_prior = as.logical(p$par_prior)),
    `QC-RLSC` = run_qc_rlsc_id_safe(
      x, fit_qc_ids, meta_method, span = as.numeric(p$span),
      degree = 1L, shift_batches = TRUE,
      min_qcs_per_batch = tq_qc_rlsc_min_qcs_per_batch(dataset)
    ),
    `QC-RFSC` = run_qc_rfsc_with_controls(
      x, fit_qc_ids, meta_method, ntree = as.integer(p$ntree),
      coCV = as.numeric(p$coCV), Frule = as.numeric(p$Frule)
    ),
    TIGER = with_no_cluster(run_tiger_all_corrected(
      x, fit_qc_ids, meta_method,
      use_injection = as.logical(p$use_injection),
      mtry_percent = 0.4, nodesize_percent = 0.4,
      ntree = 5L, parallel_cores = 1L
    )),
    SERRF = run_serrf_all_corrected(
      x, fit_qc_ids, meta_method, jitter_eps = as.numeric(p$jitter_eps)
    ),
    `WINN auto (QC)` = run_winn_mode1_auto_qc(
      x, fit_qc_ids, meta_method, params = candidate_row
    ),
    `WINN auto-batch (QC)` = run_winn_mode2_auto_batch_qc(
      x, fit_qc_ids, meta_method, params = candidate_row
    ),
    `WINN default (no QC)` = run_winn_mode3_default_no_qc(x, meta_technical),
    stop("Unsupported candidate method: ", method)
  )
}

tq_sanitize_candidate <- function(value, x, method) {
  if (is.list(value) && !is.null(value$data)) value <- value$data
  value <- as.matrix(value)
  storage.mode(value) <- "double"
  if (is.null(rownames(value)) && nrow(value) == nrow(x)) rownames(value) <- rownames(x)
  if (is.null(colnames(value)) && ncol(value) == ncol(x)) colnames(value) <- colnames(x)
  if (!identical(colnames(value), colnames(x)) || anyDuplicated(rownames(value)) ||
      !all(rownames(value) %in% rownames(x))) {
    stop(method, " candidate returned misaligned identifiers.")
  }
  value <- value[rownames(x)[rownames(x) %in% rownames(value)], , drop = FALSE]
  finite <- apply(value, 1L, function(row) all(is.finite(row)))
  if (any(!finite)) {
    if (!identical(method, "TIGER")) stop(method, " candidate returned non-finite values.")
    value <- value[finite, , drop = FALSE]
  }
  if (!nrow(value)) stop(method, " candidate retained no finite features.")
  pmax(value, 0)
}

tq_reference_profile_metrics <- function(mat, ids) {
  ids <- as.character(ids)
  ids <- ids[ids %in% colnames(mat)]
  score <- training_qc_profile_score(mat, ids)
  q <- mat[, ids, drop = FALSE]
  cvs <- if (length(ids) >= 2L) apply(q, 1L, function(v) {
    mu <- mean(v, na.rm = TRUE)
    if (!is.finite(mu) || abs(mu) <= .Machine$double.eps) NA_real_ else stats::sd(v, na.rm = TRUE) / abs(mu)
  }) else rep(NA_real_, nrow(mat))
  correlations <- if (length(ids) >= 2L) suppressWarnings(stats::cor(q, use = "pairwise.complete.obs")) else matrix(NA_real_, 0L, 0L)
  correlation_values <- if (length(correlations)) correlations[lower.tri(correlations)] else NA_real_
  cbind(
    score,
    median_reference_cv_fraction = stats::median(cvs, na.rm = TRUE),
    sd_reference_cv_fraction = stats::sd(cvs, na.rm = TRUE),
    median_reference_profile_correlation = stats::median(correlation_values, na.rm = TRUE)
  )
}

# The archived-comparator map is deliberately explicit.  Searching by filename
# prefix is unsafe here: the MTBLS79 directory contains mutually incompatible
# historical holdout/split caches, and the clinical archive is stored at the
# injection rather than adjacent-dip-pair level.  Each entry therefore states
# exactly which final archived artifact is compared with the QC-only rerun.
#
# `full_same_units` has an aligned matrix at the current correction unit.
# `clinical_injection_to_pair` is an archived injection-level matrix which is
# collapsed with the documented adjacent-dip-pair transform before comparison.
# `mtbls79_endpoint_composite` joins a study-only historical matrix with the
# matching final held-out-QC output.  This is sufficient for every requested
# external/preservation endpoint (study samples plus the fixed 12 references),
# but it is deliberately labelled endpoint-only rather than full-matrix
# equivalent because archived training-QC columns were not retained.
tq_archived_matrix_spec <- function(repo_root, dataset, method) {
  path <- function(...) file.path(repo_root, ...)
  mtbls_holdout <- function(filename) path(
    "mtbls79_plots_publication", "method_cache_qc_holdout", filename
  )
  mtbls_cache <- function(filename) path(
    "mtbls79_plots_publication", "method_cache", filename
  )
  method_paths <- list(
    mtbls79 = list(
      `Raw` = list(representation = "full_same_units", paths = c(
        matrix = mtbls_cache("raw.rds")
      )),
      `ComBat` = list(representation = "full_same_units", paths = c(
        matrix = mtbls_cache("combat.rds")
      )),
      `QC-RLSC` = list(representation = "full_same_units", paths = c(
        matrix = mtbls_holdout("qc_rlsc_holdout12_seed20260327_12761_fullrun_v5_qcrlscsubtract_v1_span_0.500.rds")
      )),
      `QC-RFSC` = list(representation = "full_same_units", paths = c(
        matrix = mtbls_holdout("qc_rfsc_holdout12_seed20260327_12761_fullrun_v5_nt200_cv20.00_fr0.80.rds")
      )),
      `TIGER` = list(representation = "mtbls79_endpoint_composite", paths = c(
        study = mtbls_cache("tiger_tigerr_default_v1.rds"),
        heldout = mtbls_holdout("tiger_qccv_holdout_holdout12_seed20260327_12761_fullrun_v5_inj_1_mtry_0.4_nodesize_0.4_ntree_5_fullrun.rds")
      )),
      `SERRF` = list(representation = "mtbls79_endpoint_composite", paths = c(
        study = mtbls_cache("serrf.rds"),
        heldout = mtbls_holdout("serrf_qccv_holdout_holdout12_seed20260327_12761_fullrun_v5_jit_0.000001_fullrun.rds")
      )),
      `WINN auto (QC)` = list(representation = "full_same_units", paths = c(
        matrix = mtbls_holdout("winn_auto_qc__holdout12_seed20260327_12761_fullrun_v5_winnauto_v3.rds")
      )),
      `WINN auto-batch (QC)` = list(representation = "full_same_units", paths = c(
        matrix = mtbls_holdout("winn_auto_batch_qc__holdout12_seed20260327_12761_fullrun_v5_winnautobatch_v3.rds")
      )),
      `WINN default (no QC)` = list(representation = "full_same_units", paths = c(
        matrix = mtbls_cache("winn_default_no_qc_.rds")
      ))
    ),
    batchcorr_set1 = list(), clinical_fiams = list(), sacurine = list(), waveica = list()
  )
  cache_specs <- list(
    batchcorr_set1 = c(
      Raw = "raw_batchcorr_set1_freshqc10_qcrlsc_corrected_v2.rds",
      ComBat = "combat_batchcorr_set1_freshqc10_qcrlsc_corrected_v2.rds",
      `QC-RLSC` = "qc_rlsc_batchcorr_set1_freshqc10_qcrlsc_corrected_v2.rds",
      `QC-RFSC` = "qc_rfsc_batchcorr_set1_freshqc10_qcrlsc_corrected_v2.rds",
      TIGER = "tiger_batchcorr_set1_freshqc10_qcrlsc_corrected_v2.rds",
      SERRF = "serrf_batchcorr_set1_freshqc10_qcrlsc_corrected_v2.rds",
      `WINN auto (QC)` = "winn_auto_qc__batchcorr_set1_freshqc10_qcrlsc_corrected_v2.rds",
      `WINN auto-batch (QC)` = "winn_auto_batch_qc__batchcorr_set1_freshqc10_qcrlsc_corrected_v2.rds",
      `WINN default (no QC)` = "winn_default_no_qc__batchcorr_set1_freshqc10_qcrlsc_corrected_v2.rds"
    ),
    clinical_fiams = c(
      Raw = "raw_clinical_dippair_v1.rds", ComBat = "combat_clinical_dippair_v1.rds",
      `QC-RLSC` = "qc_rlsc_clinical_dippair_v1.rds", `QC-RFSC` = "qc_rfsc_clinical_dippair_v1.rds",
      TIGER = "tiger_clinical_dippair_v1.rds", SERRF = "serrf_clinical_dippair_v1.rds",
      `WINN auto (QC)` = "winn_auto_qc__clinical_dippair_v1.rds",
      `WINN auto-batch (QC)` = "winn_auto_batch_qc__clinical_dippair_v1.rds",
      `WINN default (no QC)` = "winn_default_no_qc__clinical_dippair_v1.rds"
    ),
    sacurine = c(
      Raw = "raw_sacurine_hiddenqc_20260810_v1.rds", ComBat = "combat_sacurine_hiddenqc_20260810_v1.rds",
      `QC-RLSC` = "qc_rlsc_sacurine_hiddenqc_20260810_v1.rds", `QC-RFSC` = "qc_rfsc_sacurine_hiddenqc_20260810_v1.rds",
      TIGER = "tiger_sacurine_hiddenqc_20260810_v1.rds", SERRF = "serrf_sacurine_hiddenqc_20260810_v1.rds",
      `WINN auto (QC)` = "winn_auto_qc__sacurine_hiddenqc_20260810_v1.rds",
      `WINN auto-batch (QC)` = "winn_auto_batch_qc__sacurine_hiddenqc_20260810_v1.rds",
      `WINN default (no QC)` = "winn_default_no_qc__sacurine_hiddenqc_20260810_v1.rds"
    ),
    waveica = c(
      Raw = "raw_waveica_adenocarcinoma_hiddenqc_20260811_v1.rds", ComBat = "combat_waveica_adenocarcinoma_hiddenqc_20260811_v1.rds",
      `QC-RLSC` = "qc_rlsc_waveica_adenocarcinoma_hiddenqc_20260811_v1.rds", `QC-RFSC` = "qc_rfsc_waveica_adenocarcinoma_hiddenqc_20260811_v1.rds",
      TIGER = "tiger_waveica_adenocarcinoma_hiddenqc_20260811_v1.rds", SERRF = "serrf_waveica_adenocarcinoma_hiddenqc_20260811_v1.rds",
      `WINN auto (QC)` = "winn_auto_qc__waveica_adenocarcinoma_hiddenqc_20260811_v1.rds",
      `WINN auto-batch (QC)` = "winn_auto_batch_qc__waveica_adenocarcinoma_hiddenqc_20260811_v1.rds",
      `WINN default (no QC)` = "winn_default_no_qc__waveica_adenocarcinoma_hiddenqc_20260811_v1.rds"
    )
  )
  if (!dataset %in% names(method_paths) || !method %in% tq_method_order) {
    stop("Unknown archived comparator request: ", dataset, " / ", method, ".")
  }
  if (!is.null(method_paths[[dataset]][[method]])) return(method_paths[[dataset]][[method]])
  filename <- cache_specs[[dataset]][[method]]
  if (is.null(filename) || is.na(filename)) stop("No explicit archived comparator for ", dataset, " / ", method, ".")
  directory <- switch(
    dataset,
    batchcorr_set1 = path("results", "batchcorr_set1", "method_cache"),
    clinical_fiams = path("clinical_comparison_outputs", "method_cache"),
    sacurine = path("results", "sacurine", "method_cache"),
    waveica = path("results", "waveica_adenocarcinoma", "method_cache")
  )
  list(
    representation = if (dataset == "clinical_fiams") "clinical_injection_to_pair" else "full_same_units",
    paths = c(matrix = file.path(directory, filename))
  )
}

tq_load_archived_endpoint_matrix <- function(repo_root, dataset, method, prepared) {
  spec <- tq_archived_matrix_spec(repo_root, dataset, method)
  if (!all(file.exists(unname(spec$paths)))) {
    missing <- unname(spec$paths)[!file.exists(unname(spec$paths))]
    stop("Archived comparator artifact missing for ", dataset, " / ", method, ": ",
         paste(tq_relative_path(missing, repo_root), collapse = "; "))
  }
  matrix <- switch(
    spec$representation,
    full_same_units = tq_extract_matrix(spec$paths[["matrix"]]),
    clinical_injection_to_pair = {
      original <- load_winn_ablation_dataset(repo_root, "clinical_fiams")
      injected <- tq_extract_matrix(spec$paths[["matrix"]])
      prepare_clinical_dip_pair_input(injected, original$meta)$data
    },
    mtbls79_endpoint_composite = {
      study <- tq_extract_matrix(spec$paths[["study"]])
      heldout <- tq_extract_matrix(spec$paths[["heldout"]])
      study_ids <- as.character(prepared$meta$sample_id[as.logical(prepared$meta$is_study)])
      reference_ids <- as.character(prepared$hidden_ids)
      if (!all(study_ids %in% colnames(study)) || !all(reference_ids %in% colnames(heldout))) {
        stop("MTBLS79 archived endpoint composite lacks required study or held-out QC IDs for ", method, ".")
      }
      common_features <- intersect(rownames(study), rownames(heldout))
      if (!length(common_features)) stop("MTBLS79 archived endpoint composite has no shared features for ", method, ".")
      cbind(
        study[common_features, study_ids, drop = FALSE],
        heldout[common_features, reference_ids, drop = FALSE]
      )
    },
    stop("Unsupported archived comparator representation: ", spec$representation)
  )
  storage.mode(matrix) <- "double"
  required_ids <- unique(c(
    as.character(prepared$meta$sample_id[as.logical(prepared$meta$is_study)]),
    as.character(if (length(prepared$external_ids)) prepared$external_ids else prepared$hidden_ids)
  ))
  if (!all(required_ids %in% colnames(matrix))) {
    stop("Archived comparator lacks one or more required endpoint IDs for ", dataset, " / ", method, ".")
  }
  if (anyDuplicated(rownames(matrix)) || anyDuplicated(colnames(matrix)) ||
      any(!is.finite(matrix)) || any(matrix < 0)) {
    stop("Archived comparator has invalid identifiers or values for ", dataset, " / ", method, ".")
  }
  list(
    matrix = matrix,
    representation = spec$representation,
    source_paths = unname(spec$paths),
    source_sha256 = vapply(unname(spec$paths), tq_sha256_file, character(1)),
    endpoint_panel_complete = identical(sort(colnames(matrix)), sort(colnames(prepared$x))),
    endpoint_required_ids = required_ids
  )
}

# Kept as a compatibility shim for read-only callers.  Evaluation must use the
# representation-aware loader above, not this path alone.
tq_archived_matrix_path <- function(repo_root, dataset, method) {
  spec <- tq_archived_matrix_spec(repo_root, dataset, method)
  if (identical(spec$representation, "full_same_units")) spec$paths[["matrix"]] else NA_character_
}

tq_extract_matrix <- function(path) {
  if (!file.exists(path)) stop("Missing matrix cache: ", path)
  value <- readRDS(path)
  if (is.list(value) && !is.null(value$data)) value <- value$data
  value <- as.matrix(value)
  storage.mode(value) <- "double"
  value
}

tq_code_manifest <- function(repo_root, runner_path, winn_source = "local") {
  winn_paths <- if (identical(winn_source, "local")) {
    list.files(file.path(repo_root, "winn", "R"), pattern = "\\.R$", full.names = TRUE)
  } else {
    installed <- find.package("winn")
    paths <- c(
      file.path(installed, "DESCRIPTION"), file.path(installed, "NAMESPACE"),
      file.path(installed, "R", "winn.rdb"), file.path(installed, "R", "winn.rdx")
    )
    paths[file.exists(paths)]
  }
  paths <- unique(c(
    runner_path,
    file.path(repo_root, "scripts", "robustness", "training_qc_tuning_engine.R"),
    file.path(repo_root, "scripts", "robustness", "training_qc_tuning_policy.R"),
    file.path(repo_root, "scripts", "public_benchmark_method_helpers.R"),
    file.path(repo_root, "scripts", "winn_ablation_helpers.R"),
    file.path(repo_root, "clinical_comparison_helpers.R"),
    file.path(repo_root, "scripts", "weighted_pc_r2.R"),
    winn_paths
  ))
  paths <- paths[file.exists(paths)]
  data.frame(
    path = vapply(paths, tq_relative_path, character(1), repo_root = repo_root),
    sha256 = vapply(paths, tq_sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
}
