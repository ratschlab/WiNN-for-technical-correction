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
