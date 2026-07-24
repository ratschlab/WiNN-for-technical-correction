# Reduced correction engine for the reviewer-requested partial-confounding study.
#
# This deliberately runs only the uncorrected Raw baseline and the manuscript's
# fixed, QC-free WiNN mode.  The general 18-unit engine remains unchanged for
# the canonical benchmark and ablation studies.

run_partial_confounding_winn_only_methods <- function(
  x,
  meta_blind,
  seed_ledger,
  method_parameters,
  run_dir,
  global_context,
  force = FALSE,
  log_line = function(...) invisible(NULL)
) {
  expected_methods <- c("Raw", "WINN default (no QC)")
  if (!identical(method_parameters$method, expected_methods)) {
    stop("WiNN-only method ledger must contain Raw then WINN default (no QC).",
         call. = FALSE)
  }
  for (subdir in c("cache", "logs", "matrices", "diagnostics")) {
    dir.create(file.path(run_dir, subdir), recursive = TRUE,
               showWarnings = FALSE)
  }

  method_specs <- list(
    "Raw" = function() x,
    "WINN default (no QC)" = function() run_winn_mode3_default_no_qc(
      x, meta_df = meta_blind
    )
  )
  method_matrices <- list()
  status_rows <- list()
  error_rows <- list()
  runtime_rows <- list()
  cache_rows <- list()
  log_rows <- list()

  for (method in expected_methods) {
    parameter <- canonical_seed_parameter_row(method_parameters, method)
    rng_seed <- canonical_seed_rng(
      seed_ledger, parameter$rng_component[[1L]]
    )
    context <- c(global_context, list(
      analysis_unit = method,
      parameters = as.list(parameter),
      rng_seed = rng_seed
    ))
    cache_path <- file.path(
      run_dir, "cache", paste0(canonical_seed_safe_name(method), ".rds")
    )
    log_line("Running ", method, " (seed ", rng_seed, ").")
    captured <- canonical_seed_capture(function() {
      set.seed(rng_seed)
      canonical_cached_call(
        cache_path, context, method_specs[[method]], force = force
      )
    })

    if (inherits(captured$value, "canonical_seed_method_error")) {
      error <- captured$value
      status_rows[[method]] <- data.frame(
        method = method, status = "failed", error = error$message,
        attempted = TRUE, stringsAsFactors = FALSE
      )
      error_rows[[method]] <- data.frame(
        method = method, error = error$message, call = error$call,
        error_class = paste(error$class, collapse = ";"),
        stringsAsFactors = FALSE
      )
      runtime_rows[[method]] <- data.frame(
        method = method, runtime_sec = NA_real_, cache_hit = FALSE,
        stringsAsFactors = FALSE
      )
      log_line(method, " failed: ", error$message)
    } else {
      cached <- captured$value
      validated <- canonical_seed_capture(function() {
        canonical_seed_validate_matrix(cached$value, method, x)
      })
      if (inherits(validated$value, "canonical_seed_method_error")) {
        error <- validated$value
        status_rows[[method]] <- data.frame(
          method = method, status = "failed_validation", error = error$message,
          attempted = TRUE, stringsAsFactors = FALSE
        )
        error_rows[[method]] <- data.frame(
          method = method, error = error$message, call = error$call,
          error_class = paste(error$class, collapse = ";"),
          stringsAsFactors = FALSE
        )
      } else {
        method_matrices[[method]] <- validated$value
        status_rows[[method]] <- data.frame(
          method = method,
          status = if (cached$cache_hit) "completed_cached" else "completed",
          error = "", attempted = TRUE, stringsAsFactors = FALSE
        )
        cache_rows[[method]] <- data.frame(
          method = method, cache_key = cached$cache_key,
          cache_hit = cached$cache_hit, cache_reason = cached$cache_reason,
          stringsAsFactors = FALSE
        )
      }
      runtime_rows[[method]] <- data.frame(
        method = method, runtime_sec = cached$runtime_sec,
        cache_hit = cached$cache_hit, stringsAsFactors = FALSE
      )
      if (!inherits(validated$value, "canonical_seed_method_error")) {
        log_line(
          method, " completed in ", sprintf("%.3f", cached$runtime_sec),
          " seconds; cache_hit=", cached$cache_hit
        )
      }
    }
    log_rows[[method]] <- data.frame(
      method = method,
      warnings = paste(captured$warnings, collapse = "\n"),
      messages = paste(captured$messages, collapse = "\n"),
      stringsAsFactors = FALSE
    )
    writeLines(
      c("WARNINGS", captured$warnings, "", "MESSAGES", captured$messages),
      file.path(
        run_dir, "logs", paste0(canonical_seed_safe_name(method), ".txt")
      )
    )
  }

  ordered_bind <- function(rows) {
    bound <- dplyr::bind_rows(rows)
    bound[match(expected_methods, bound$method), , drop = FALSE]
  }
  status <- ordered_bind(status_rows)
  runtime <- ordered_bind(runtime_rows)
  errors <- if (length(error_rows)) {
    dplyr::bind_rows(error_rows)
  } else {
    data.frame(
      method = character(), error = character(), call = character(),
      error_class = character(), stringsAsFactors = FALSE
    )
  }
  cache_manifest <- dplyr::bind_rows(cache_rows)
  warnings_messages <- ordered_bind(log_rows)

  seed_id <- as.character(global_context$seed_id)
  add_seed_id <- function(value) {
    value$seed_id <- rep(seed_id, nrow(value))
    value[, c("seed_id", setdiff(names(value), "seed_id")), drop = FALSE]
  }
  status <- add_seed_id(status)
  runtime <- add_seed_id(runtime)
  errors <- add_seed_id(errors)
  cache_manifest <- add_seed_id(cache_manifest)
  warnings_messages <- add_seed_id(warnings_messages)

  matrix_rows <- list()
  for (method in names(method_matrices)) {
    path <- file.path(
      run_dir, "matrices",
      paste0(canonical_seed_safe_name(method), ".rds")
    )
    saveRDS(method_matrices[[method]], path, version = 3, compress = "gzip")
    matrix_rows[[method]] <- data.frame(
      seed_id = seed_id, method = method,
      relative_path = file.path("matrices", basename(path)),
      n_features = nrow(method_matrices[[method]]),
      n_samples = ncol(method_matrices[[method]]),
      bytes = file.info(path)$size,
      sha256 = canonical_sha256_file(path),
      stringsAsFactors = FALSE
    )
  }
  matrix_manifest <- dplyr::bind_rows(matrix_rows)

  list(
    matrices = method_matrices,
    status = status,
    runtime = runtime,
    errors = errors,
    cache_manifest = cache_manifest,
    warnings_messages = warnings_messages,
    matrix_manifest = matrix_manifest,
    ablation_diagnostics = list(),
    ablation_configurations = list()
  )
}
