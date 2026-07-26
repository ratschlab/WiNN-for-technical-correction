release_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop(name, " is not set.", call. = FALSE)
  normalizePath(value, mustWork = TRUE)
}

release_root <- release_env("WINN_RELEASE_ROOT")
source_root <- release_env("WINN_CANONICAL_SOURCE_ROOT")

source(file.path(release_root, "package", "assert_frozen_winn.R"), local = FALSE)
frozen_winn <- assert_frozen_winn(release_root)

required_packages <- c(
  "digest", "jsonlite", "dplyr", "tibble", "sva", "qcrlscR",
  "statTarget", "TIGERr", "malbacR", "pmartR", "winn"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
source(file.path(source_root, "clinical_comparison_helpers.R"), local = FALSE)
source(file.path(release_root, "scripts", "benchmark_helpers.R"), local = FALSE)
source(file.path(source_root, "scripts", "public_benchmark_method_helpers.R"), local = FALSE)
source(file.path(source_root, "scripts", "weighted_pc_r2.R"), local = FALSE)
source(file.path(source_root, "scripts", "winn_ablation_helpers.R"), local = FALSE)
source(file.path(source_root, "scripts", "robustness", "training_qc_tuning_engine.R"), local = FALSE)

release_sha256_object <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

release_sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256")
}

release_atomic_save_rds <- function(value, path, compress = "gzip") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(basename(path), ".tmp_"), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(value, temporary, version = 3, compress = compress)
  if (!file.rename(temporary, path)) stop("Could not atomically save ", path, call. = FALSE)
  invisible(path)
}

release_capture <- function(fn) {
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
    runtime_sec = proc.time()[["elapsed"]] - started,
    warnings = unique(warnings),
    messages = unique(messages),
    error = if (inherits(value, "error")) conditionMessage(value) else ""
  )
}

release_load_simulation <- function(seed_id = "SIM01") {
  bundle <- file.path(source_root, "data", "simulated", "canonical", seed_id)
  x <- readRDS(file.path(bundle, "raw_intensity.rds"))
  truth <- readRDS(file.path(bundle, "clean_ground_truth.rds"))
  meta <- read.csv(
    file.path(bundle, "sample_metadata.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  hidden_ids <- as.character(meta$sample_id[as.logical(meta$is_hidden_reference)])
  training_ids <- as.character(meta$sample_id[as.logical(meta$is_qc) & !as.logical(meta$is_hidden_reference)])
  meta$is_study <- !as.logical(meta$is_qc)
  meta$biological_label <- ifelse(meta$is_study, "study", "QC")
  list(
    x = x, meta = meta, truth = truth,
    training_ids = training_ids, hidden_ids = hidden_ids,
    external_ids = character(), correction_unit = "injection",
    bundle = bundle
  )
}

release_load_dataset <- function(dataset_key) {
  if (identical(dataset_key, "simulation")) {
    prepared <- release_load_simulation("SIM01")
  } else {
    prepared <- tq_prepare_dataset(source_root, dataset_key)
  }
  prepared$x <- as.matrix(prepared$x)
  storage.mode(prepared$x) <- "double"
  prepared$meta <- prepared$meta[
    match(colnames(prepared$x), prepared$meta$sample_id), , drop = FALSE
  ]
  if (anyNA(prepared$meta$sample_id) ||
      !identical(colnames(prepared$x), as.character(prepared$meta$sample_id))) {
    stop(dataset_key, ": matrix and metadata are not aligned.", call. = FALSE)
  }
  if (any(!is.finite(prepared$x)) || any(prepared$x < 0) ||
      anyDuplicated(rownames(prepared$x)) || anyDuplicated(colnames(prepared$x))) {
    stop(dataset_key, ": canonical input validation failed.", call. = FALSE)
  }
  prepared$meta_method <- tq_method_metadata(prepared$meta, prepared$training_ids)
  prepared$meta_technical <- tq_technical_metadata(prepared$meta)
  prepared
}

release_expected_input_hash <- function(dataset_key) {
  manifest <- read.csv(
    file.path(release_root, "analysis", "config", "dataset_manifest.csv"),
    stringsAsFactors = FALSE
  )
  value <- manifest$matrix_sha256[match(dataset_key, manifest$dataset_key)]
  if (length(value) != 1L || is.na(value)) stop("No frozen input hash for ", dataset_key)
  value
}

release_reference_assignments <- function(dataset_key) {
  path <- if (identical(dataset_key, "clinical_fiams")) {
    file.path(
      source_root, "data", "private", "clinical_fiams",
      "reference_split_assignments.csv"
    )
  } else {
    file.path(release_root, "analysis", "config", "reference_split_assignments.csv")
  }
  if (!file.exists(path)) {
    stop("Missing reference-split assignment file for ", dataset_key, ": ", path, call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

release_validate_input <- function(prepared, dataset_key) {
  observed <- release_sha256_object(prepared$x)
  expected <- release_expected_input_hash(dataset_key)
  if (!identical(observed, expected)) {
    stop(
      dataset_key, ": canonical input hash mismatch; expected ", expected,
      ", observed ", observed, ".", call. = FALSE
    )
  }
  invisible(observed)
}

release_selected_parameters <- function(dataset_key, method_id) {
  manifest <- read.csv(
    file.path(release_root, "analysis", "config", "endpoint_free_selection_manifest.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  row <- manifest[manifest$dataset_key == dataset_key & manifest$method_id == method_id, , drop = FALSE]
  if (nrow(row) != 1L) stop("Missing unique selected configuration for ", dataset_key, " / ", method_id)
  row
}

release_reuse_decision <- function(dataset_key, method_id) {
  "rerun_required"
}

release_run_winn <- function(prepared, auto_batch, parameters) {
  controls <- if (identical(parameters, "auto")) {
    which(prepared$meta$sample_id %in% prepared$training_ids)
  } else {
    NULL
  }
  winn::winn(
    data = prepared$x,
    batch = if (isTRUE(auto_batch)) NULL else prepared$meta$batch,
    run_order = prepared$meta$run_order,
    control_samples = controls,
    parameters = parameters,
    ljung_box_fitdf = 0L,
    return_details = TRUE
  )
}

release_run_method <- function(dataset_key, method_id, prepared) {
  cfg <- release_selected_parameters(dataset_key, method_id)
  value <- switch(
    method_id,
    combat = run_combat(
      prepared$x, prepared$meta_technical,
      par_prior = grepl("par_prior=TRUE", cfg$selected_parameters, fixed = TRUE)
    ),
    qc_rlsc = run_qc_rlsc_id_safe(
      prepared$x, prepared$training_ids, prepared$meta_method,
      span = as.numeric(sub(".*span=([0-9.]+).*", "\\1", cfg$selected_parameters)),
      degree = 1L, shift_batches = TRUE,
      min_qcs_per_batch = if (identical(dataset_key, "mtbls79")) 3L else 4L
    ),
    qc_rfsc = {
      ntree <- as.integer(sub(".*ntree=([0-9]+).*", "\\1", cfg$selected_parameters))
      cocv <- as.numeric(sub(".*coCV=([0-9.]+).*", "\\1", cfg$selected_parameters))
      frule <- as.numeric(sub(".*Frule=([0-9.]+).*", "\\1", cfg$selected_parameters))
      run_qc_rfsc_with_controls(
        prepared$x, prepared$training_ids, prepared$meta_method,
        ntree = ntree, coCV = cocv, Frule = frule
      )
    },
    tiger = with_no_cluster(run_tiger_all_corrected(
      prepared$x, prepared$training_ids, prepared$meta_method,
      use_injection = TRUE,
      mtry_percent = seq(0.2, 0.8, by = 0.2),
      nodesize_percent = seq(0.2, 0.8, by = 0.2),
      ntree = 500L,
      parallel_cores = as.integer(Sys.getenv("WINN_TIGER_CORES", unset = "1"))
    )),
    serrf = {
      jitter <- as.numeric(sub(".*jitter_eps=([0-9.eE+-]+).*", "\\1", cfg$selected_parameters))
      run_serrf_all_corrected(
        prepared$x, prepared$training_ids, prepared$meta_method,
        jitter_eps = jitter
      )
    },
    winn_auto_qc = release_run_winn(prepared, auto_batch = FALSE, parameters = "auto"),
    winn_auto_batch_qc = release_run_winn(prepared, auto_batch = TRUE, parameters = "auto"),
    winn_fixed_no_qc = release_run_winn(prepared, auto_batch = FALSE, parameters = "fixed"),
    raw = prepared$x,
    stop("Unsupported method_id: ", method_id, call. = FALSE)
  )
  if (inherits(value, "winn_result")) {
    list(data = value$data, details = unclass(value)[setdiff(names(value), "data")])
  } else {
    list(data = value, details = list(reuse_decision = "rerun_required"))
  }
}

release_validate_output <- function(value, prepared, method_id) {
  tiger_fallback <- attr(value, "tiger_nonfinite_fallback", exact = TRUE)
  value <- as.matrix(value)
  storage.mode(value) <- "double"
  if (!identical(dim(value), dim(prepared$x)) ||
      !identical(rownames(value), rownames(prepared$x)) ||
      !identical(colnames(value), colnames(prepared$x))) {
    stop(method_id, ": output dimensions or identifiers changed.", call. = FALSE)
  }
  if (any(!is.finite(value))) {
    stop(method_id, ": output contains non-finite values.", call. = FALSE)
  }
  negative <- value < 0
  if (any(negative)) {
    floor_record <- list(
      applied = TRUE,
      rule = "Negative corrected intensities were floored at zero before evaluation.",
      n_values = sum(negative),
      n_features = sum(rowSums(negative) > 0L),
      n_samples = sum(colSums(negative) > 0L),
      minimum_before_floor = min(value)
    )
    value[negative] <- 0
    attr(value, "intensity_floor") <- floor_record
  }
  if (!is.null(tiger_fallback)) {
    attr(value, "tiger_nonfinite_fallback") <- tiger_fallback
  }
  value
}
