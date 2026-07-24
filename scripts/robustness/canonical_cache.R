# Content-addressed cache helpers for the robustness simulation workstream.
#
# Cache identity is deliberately derived from serialized input content rather
# than a human-readable method label.  A cached value is reused only when the
# complete expected key matches the key stored in its envelope.

canonical_sha256_file <- function(path) {
  if (!file.exists(path)) stop("Cannot hash missing file: ", path, call. = FALSE)
  digest::digest(path, algo = "sha256", file = TRUE)
}

canonical_sha256_object <- function(object) {
  digest::digest(object, algo = "sha256", serialize = TRUE)
}

canonical_cache_key <- function(context) {
  canonical_sha256_object(list(cache_schema = "canonical-simulation-v1", context = context))
}

canonical_cache_read <- function(path, expected_key) {
  if (!file.exists(path)) {
    return(list(hit = FALSE, reason = "missing", value = NULL, envelope = NULL))
  }
  envelope <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(envelope, "error")) {
    return(list(
      hit = FALSE,
      reason = paste0("unreadable: ", conditionMessage(envelope)),
      value = NULL,
      envelope = NULL
    ))
  }
  required <- c("cache_schema", "key", "created_at", "value")
  if (!is.list(envelope) || !all(required %in% names(envelope))) {
    return(list(hit = FALSE, reason = "invalid_envelope", value = NULL, envelope = envelope))
  }
  if (!identical(envelope$cache_schema, "canonical-simulation-v1")) {
    return(list(hit = FALSE, reason = "schema_mismatch", value = NULL, envelope = envelope))
  }
  if (!identical(envelope$key, expected_key)) {
    return(list(hit = FALSE, reason = "key_mismatch", value = NULL, envelope = envelope))
  }
  list(hit = TRUE, reason = "exact_key_match", value = envelope$value, envelope = envelope)
}

canonical_cache_write <- function(path, key, value, context, runtime_sec) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary, force = TRUE), add = TRUE)
  envelope <- list(
    cache_schema = "canonical-simulation-v1",
    key = key,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    context_sha256 = canonical_sha256_object(context),
    runtime_sec = as.numeric(runtime_sec),
    value = value
  )
  saveRDS(envelope, temporary, version = 3, compress = "gzip")
  if (!file.rename(temporary, path)) stop("Failed to commit cache file: ", path, call. = FALSE)
  invisible(envelope)
}

canonical_cached_call <- function(path, context, fn, force = FALSE) {
  key <- canonical_cache_key(context)
  cached <- if (isTRUE(force)) {
    list(hit = FALSE, reason = "force_requested", value = NULL, envelope = NULL)
  } else {
    canonical_cache_read(path, key)
  }
  if (isTRUE(cached$hit)) {
    return(list(
      value = cached$value,
      cache_hit = TRUE,
      cache_reason = cached$reason,
      cache_key = key,
      runtime_sec = as.numeric(cached$envelope$runtime_sec)
    ))
  }

  started <- proc.time()[["elapsed"]]
  value <- fn()
  elapsed <- proc.time()[["elapsed"]] - started
  canonical_cache_write(path, key, value, context, elapsed)
  list(
    value = value,
    cache_hit = FALSE,
    cache_reason = cached$reason,
    cache_key = key,
    runtime_sec = elapsed
  )
}

run_canonical_cache_negative_tests <- function(output_path) {
  base_context <- list(input_sha256 = "input-a", method = "Raw", parameters = list(mode = "raw"))
  key_a <- canonical_cache_key(base_context)
  key_a_repeat <- canonical_cache_key(base_context)
  parameter_context <- base_context
  parameter_context$parameters$mode <- "changed"
  input_context <- base_context
  input_context$input_sha256 <- "input-b"

  temporary <- tempfile("canonical_cache_negative_", fileext = ".rds")
  on.exit(if (file.exists(temporary)) unlink(temporary, force = TRUE), add = TRUE)
  canonical_cache_write(temporary, key_a, matrix(1, 1, 1), base_context, 0)
  stale_read <- canonical_cache_read(temporary, canonical_cache_key(parameter_context))
  valid_read <- canonical_cache_read(temporary, key_a)

  tests <- data.frame(
    check = c(
      "same_context_same_key",
      "parameter_change_changes_key",
      "input_change_changes_key",
      "stale_envelope_is_rejected",
      "exact_envelope_is_accepted"
    ),
    passed = c(
      identical(key_a, key_a_repeat),
      !identical(key_a, canonical_cache_key(parameter_context)),
      !identical(key_a, canonical_cache_key(input_context)),
      !isTRUE(stale_read$hit) && identical(stale_read$reason, "key_mismatch"),
      isTRUE(valid_read$hit)
    ),
    observed = c(
      key_a,
      canonical_cache_key(parameter_context),
      canonical_cache_key(input_context),
      stale_read$reason,
      valid_read$reason
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(tests, output_path, row.names = FALSE, quote = TRUE)
  if (!all(tests$passed)) stop("Content-addressed cache negative tests failed.", call. = FALSE)
  invisible(tests)
}
