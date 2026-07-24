# Method execution engine for one canonical simulation seed.
#
# The caller must load the local WiNN source and source the shared benchmark,
# public-method, ablation, and content-cache helpers before calling this file.

canonical_seed_safe_name <- function(value) {
  cleaned <- gsub("[^A-Za-z0-9]+", "_", value)
  gsub("^_|_$", "", cleaned)
}

# jsonlite::toJSON() returns an object with class "json".  Keep serialized
# parameters as an ordinary character scalar so rows originating from frozen
# CSV parameters and rows originating from auto-selection can always be bound
# into one audit table without a vctrs class conflict.
canonical_seed_plain_json <- function(...) {
  as.character(jsonlite::toJSON(...))
}

canonical_seed_capture <- function(fn) {
  warnings <- character()
  messages <- character()
  error_call <- ""
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
    error = function(e) {
      error_call <<- paste(deparse(conditionCall(e)), collapse = " ")
      structure(
        list(message = conditionMessage(e), call = error_call, class = class(e)),
        class = "canonical_seed_method_error"
      )
    }
  )
  list(value = value, warnings = unique(warnings), messages = unique(messages))
}

canonical_seed_validate_matrix <- function(value, method, template) {
  if (!is.matrix(value) || !is.numeric(value)) stop(method, " did not return a numeric matrix.", call. = FALSE)
  if (!identical(dim(value), dim(template))) stop(method, " output dimensions differ from the canonical input.", call. = FALSE)
  if (!identical(rownames(value), rownames(template))) stop(method, " changed feature identity/order.", call. = FALSE)
  if (!identical(colnames(value), colnames(template))) stop(method, " changed injection identity/order.", call. = FALSE)
  if (any(!is.finite(value))) stop(method, " returned nonfinite values.", call. = FALSE)
  if (any(value <= -1)) stop(method, " returned values <= -1, invalid for log1p evaluation.", call. = FALSE)
  value
}

canonical_seed_parameter_row <- function(method_parameters, method) {
  row <- method_parameters[method_parameters$method == method, , drop = FALSE]
  if (nrow(row) != 1L) stop("Missing or duplicate frozen parameters for ", method, call. = FALSE)
  row
}

canonical_seed_rng <- function(seed_ledger, component) {
  row <- seed_ledger[seed_ledger$component == component, , drop = FALSE]
  if (nrow(row) != 1L) stop("Missing or duplicate seed component: ", component, call. = FALSE)
  as.integer(row$effective_seed[[1L]])
}

run_canonical_seed_methods <- function(
  x,
  meta_blind,
  training_ids,
  seed_ledger,
  method_parameters,
  run_dir,
  global_context,
  force = FALSE,
  log_line = function(...) invisible(NULL)
) {
  for (subdir in c("cache", "logs", "matrices", "diagnostics")) {
    dir.create(file.path(run_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  }

  method_specs <- list(
    "Raw" = function() x,
    "ComBat" = function() run_combat(x, meta_df = meta_blind, par_prior = TRUE),
    "QC-RLSC" = function() run_qc_rlsc_id_safe(
      x, control_ids = training_ids, meta_df = meta_blind,
      span = 0.75, degree = 1L, shift_batches = TRUE
    ),
    "QC-RFSC" = function() run_qc_rfsc_with_controls(
      x, control_ids = training_ids, meta_df = meta_blind,
      ntree = 500L, coCV = 30, Frule = 0.8
    ),
    "TIGER" = function() with_no_cluster(run_tiger_all_corrected(
      x, control_ids = training_ids, meta_df = meta_blind,
      use_injection = TRUE, mtry_percent = 0.4, nodesize_percent = 0.4,
      ntree = 5L, parallel_cores = 1L
    )),
    "SERRF" = function() run_serrf_all_corrected(
      x, control_ids = training_ids, meta_df = meta_blind, jitter_eps = 0
    ),
    "WINN auto (QC)" = function() run_winn_mode1_auto_qc(
      x, control_ids = training_ids, meta_df = meta_blind, params = NULL
    ),
    "WINN auto-batch (QC)" = function() run_winn_mode2_auto_batch_qc(
      x, control_ids = training_ids, meta_df = meta_blind, params = NULL
    ),
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

  for (method in names(method_specs)) {
    parameter <- canonical_seed_parameter_row(method_parameters, method)
    rng_seed <- canonical_seed_rng(seed_ledger, parameter$rng_component[[1L]])
    context <- c(global_context, list(
      analysis_unit = method,
      parameters = as.list(parameter),
      rng_seed = rng_seed
    ))
    cache_path <- file.path(run_dir, "cache", paste0(canonical_seed_safe_name(method), ".rds"))
    log_line("Running ", method, " (seed ", rng_seed, ").")
    captured <- canonical_seed_capture(function() {
      set.seed(rng_seed)
      canonical_cached_call(cache_path, context, method_specs[[method]], force = force)
    })

    if (inherits(captured$value, "canonical_seed_method_error")) {
      error <- captured$value
      status_rows[[method]] <- data.frame(
        method = method, status = "failed", error = error$message,
        attempted = TRUE, stringsAsFactors = FALSE
      )
      error_rows[[method]] <- data.frame(
        method = method, error = error$message, call = error$call,
        error_class = paste(error$class, collapse = ";"), stringsAsFactors = FALSE
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
          error_class = paste(error$class, collapse = ";"), stringsAsFactors = FALSE
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
        log_line(method, " completed in ", sprintf("%.3f", cached$runtime_sec),
                 " seconds; cache_hit=", cached$cache_hit)
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
      file.path(run_dir, "logs", paste0(canonical_seed_safe_name(method), ".txt"))
    )
  }

  ablation_context <- c(global_context, list(
    analysis_unit = "WiNN fixed ablation bundle",
    parameters = as.list(method_parameters[method_parameters$analysis_family != "cross_method", , drop = FALSE]),
    rng_seed = canonical_seed_rng(seed_ledger, "method_rng::WINN fixed (QC identities withheld)")
  ))
  ablation_cache_path <- file.path(run_dir, "cache", "WiNN_ablation_bundle.rds")
  log_line("Running WiNN fixed ablation bundle.")
  ablation_capture <- canonical_seed_capture(function() {
    set.seed(canonical_seed_rng(seed_ledger, "method_rng::WINN fixed (QC identities withheld)"))
    canonical_cached_call(
      ablation_cache_path,
      ablation_context,
      function() {
        base <- winn::winn_ablation(
          x, batch = meta_blind$batch, run_order = meta_blind$run_order,
          control_samples = NULL, parameters = "fixed",
          use_outlier_shrinkage = TRUE, drift_gate = "selective",
          batch_gate = "selective", pqn_mode = "shrink",
          fdr_threshold = 0.05, test = "Ljung-Box", lag = NULL,
          spline_method = "conservative", return_intermediates = TRUE,
          return_diagnostics = TRUE
        )
        all_drift <- winn::winn_ablation(
          x, batch = meta_blind$batch, run_order = meta_blind$run_order,
          control_samples = NULL, parameters = "fixed",
          use_outlier_shrinkage = TRUE, drift_gate = "all",
          batch_gate = "selective", pqn_mode = "shrink",
          fdr_threshold = 0.05, test = "Ljung-Box", lag = NULL,
          spline_method = "conservative", return_intermediates = TRUE,
          return_diagnostics = TRUE
        )
        selective_drift_all_batch <- derive_forced_batch_variant(
          base, meta_blind$batch, control_samples = NULL
        )
        all_drift_all_batch <- derive_forced_batch_variant(
          all_drift, meta_blind$batch, control_samples = NULL
        )
        variants <- list(
          C0_RAW = empty_raw_result(x),
          C1_OUTLIER = subset_cumulative_result(base, "outlier"),
          C2_SELECTIVE_DRIFT = subset_cumulative_result(base, "drift"),
          C3_SELECTIVE_BATCH = subset_cumulative_result(base, "batch"),
          C4_FULL_FIXED = base,
          G_SS = base,
          G_AS = all_drift,
          G_SA = selective_drift_all_batch,
          G_AA = all_drift_all_batch
        )
        list(
          matrices = lapply(variants, `[[`, "data"),
          runtime_sec = vapply(variants, `[[`, numeric(1), "total_runtime_sec"),
          configurations = lapply(variants, `[[`, "configuration"),
          diagnostics = lapply(variants, `[[`, "diagnostics")
        )
      },
      force = force
    )
  })

  ablation_methods <- method_parameters$method[
    method_parameters$analysis_family != "cross_method"
  ]
  ablation_diagnostics <- list()
  ablation_configurations <- list()
  if (inherits(ablation_capture$value, "canonical_seed_method_error")) {
    error <- ablation_capture$value
    for (method in ablation_methods) {
      status_rows[[method]] <- data.frame(
        method = method, status = "failed", error = error$message,
        attempted = TRUE, stringsAsFactors = FALSE
      )
      error_rows[[method]] <- data.frame(
        method = method, error = error$message, call = error$call,
        error_class = paste(error$class, collapse = ";"), stringsAsFactors = FALSE
      )
      runtime_rows[[method]] <- data.frame(
        method = method, runtime_sec = NA_real_, cache_hit = FALSE,
        stringsAsFactors = FALSE
      )
      log_rows[[method]] <- data.frame(
        method = method,
        warnings = paste(ablation_capture$warnings, collapse = "\n"),
        messages = paste(ablation_capture$messages, collapse = "\n"),
        stringsAsFactors = FALSE
      )
      writeLines(
        c("WARNINGS", ablation_capture$warnings, "", "MESSAGES", ablation_capture$messages),
        file.path(run_dir, "logs", paste0(canonical_seed_safe_name(method), ".txt"))
      )
    }
    log_line("WiNN ablation bundle failed: ", error$message)
  } else {
    cached <- ablation_capture$value
    ablation <- cached$value
    ablation_diagnostics <- ablation$diagnostics
    ablation_configurations <- ablation$configurations
    saveRDS(
      ablation_diagnostics,
      file.path(run_dir, "diagnostics", "winn_ablation_diagnostics.rds"),
      version = 3, compress = "gzip"
    )
    for (method in names(ablation$matrices)) {
      validated <- canonical_seed_capture(function() {
        canonical_seed_validate_matrix(ablation$matrices[[method]], method, x)
      })
      if (inherits(validated$value, "canonical_seed_method_error")) {
        error <- validated$value
        status_rows[[method]] <- data.frame(
          method = method, status = "failed_validation", error = error$message,
          attempted = TRUE, stringsAsFactors = FALSE
        )
        error_rows[[method]] <- data.frame(
          method = method, error = error$message, call = error$call,
          error_class = paste(error$class, collapse = ";"), stringsAsFactors = FALSE
        )
      } else {
        method_matrices[[method]] <- validated$value
        status_rows[[method]] <- data.frame(
          method = method,
          status = if (cached$cache_hit) "completed_cached" else "completed",
          error = "", attempted = TRUE, stringsAsFactors = FALSE
        )
      }
      runtime_rows[[method]] <- data.frame(
        method = method, runtime_sec = unname(ablation$runtime_sec[[method]]),
        cache_hit = cached$cache_hit, stringsAsFactors = FALSE
      )
      cache_rows[[method]] <- data.frame(
        method = method, cache_key = cached$cache_key,
        cache_hit = cached$cache_hit, cache_reason = cached$cache_reason,
        stringsAsFactors = FALSE
      )
      log_rows[[method]] <- data.frame(
        method = method,
        warnings = paste(ablation_capture$warnings, collapse = "\n"),
        messages = paste(ablation_capture$messages, collapse = "\n"),
        stringsAsFactors = FALSE
      )
      writeLines(
        c("WARNINGS", ablation_capture$warnings, "", "MESSAGES", ablation_capture$messages),
        file.path(run_dir, "logs", paste0(canonical_seed_safe_name(method), ".txt"))
      )
    }
    log_line("WiNN ablation bundle completed in ", sprintf("%.3f", cached$runtime_sec),
             " seconds; cache_hit=", cached$cache_hit)
  }

  ordered_bind <- function(rows) {
    bound <- dplyr::bind_rows(rows)
    bound[match(method_parameters$method, bound$method), , drop = FALSE]
  }
  status <- ordered_bind(status_rows)
  runtime <- ordered_bind(runtime_rows)
  errors <- if (length(error_rows)) {
    dplyr::bind_rows(error_rows)
  } else {
    data.frame(method = character(), error = character(), call = character(), error_class = character())
  }
  cache_manifest <- dplyr::bind_rows(cache_rows)
  warnings_messages <- dplyr::bind_rows(log_rows)

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

  write.csv(status, file.path(run_dir, "method_status.csv"), row.names = FALSE, quote = TRUE)
  write.csv(runtime, file.path(run_dir, "runtime.csv"), row.names = FALSE, quote = TRUE)
  write.csv(errors, file.path(run_dir, "errors.csv"), row.names = FALSE, quote = TRUE)
  write.csv(cache_manifest, file.path(run_dir, "cache_manifest.csv"), row.names = FALSE, quote = TRUE)
  write.csv(warnings_messages, file.path(run_dir, "warnings_messages_errors.csv"), row.names = FALSE, quote = TRUE)

  matrix_rows <- list()
  for (method in names(method_matrices)) {
    path <- file.path(run_dir, "matrices", paste0(canonical_seed_safe_name(method), ".rds"))
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
  write.csv(matrix_manifest, file.path(run_dir, "matrix_checksums.csv"), row.names = FALSE, quote = TRUE)

  list(
    matrices = method_matrices,
    status = status,
    runtime = runtime,
    errors = errors,
    cache_manifest = cache_manifest,
    warnings_messages = warnings_messages,
    matrix_manifest = matrix_manifest,
    ablation_diagnostics = ablation_diagnostics,
    ablation_configurations = ablation_configurations
  )
}
