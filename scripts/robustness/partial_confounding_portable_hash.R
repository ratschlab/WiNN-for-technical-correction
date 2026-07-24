# Cross-platform reconstruction identity for the partial-confounding study.
#
# Exact R-object hashes remain useful provenance diagnostics, but floating-point
# elementary functions can differ by a few ulps across operating systems and R
# patch releases. Exponentiation can amplify those harmless differences onto a
# decimal rounding boundary. Scientific identity is therefore gated before
# exponentiation, using the complete raw/clean log matrices after decimal
# quantization. Intensity object and intensity-round9 hashes are diagnostics.
# Matrix dimensions and dimnames remain part of the serialized digest.

.partial_confounding_portable_hash_digits <- 9L
.partial_confounding_portable_max_abs_tolerance <-
  10^(-.partial_confounding_portable_hash_digits)
.partial_confounding_portable_hash_schema <-
  "r_serialized_pre_exp_log_matrices_rounded_9_decimals_v3"

partial_confounding_quantized_matrix <- function(
  matrix,
  digits = .partial_confounding_portable_hash_digits
) {
  if (!is.matrix(matrix) || !is.numeric(matrix)) {
    stop("Portable hashing requires a numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(matrix))) {
    stop("Portable hashing requires finite matrix values.", call. = FALSE)
  }
  if (length(digits) != 1L || is.na(digits) ||
      as.integer(digits) != digits || digits < 0L || digits > 15L) {
    stop("digits must be one integer from 0 through 15.", call. = FALSE)
  }
  round(matrix, digits = as.integer(digits))
}

partial_confounding_portable_matrix_sha256 <- function(
  matrix,
  digits = .partial_confounding_portable_hash_digits
) {
  if (!exists("canonical_sha256_object", mode = "function")) {
    stop("Source canonical_cache.R before portable hashing.", call. = FALSE)
  }
  canonical_sha256_object(partial_confounding_quantized_matrix(
    matrix, digits = digits
  ))
}

partial_confounding_matrix_identity <- function(
  matrix,
  expected_authoritative_quantized_sha256,
  reference_object_sha256,
  digits = .partial_confounding_portable_hash_digits
) {
  valid_sha256 <- function(value) {
    length(value) == 1L && !is.na(value) &&
      grepl("^[0-9a-f]{64}$", value)
  }
  if (!valid_sha256(expected_authoritative_quantized_sha256) ||
      !valid_sha256(reference_object_sha256)) {
    stop("Expected portable and reference hashes must be SHA-256 strings.",
         call. = FALSE)
  }
  observed_object <- canonical_sha256_object(matrix)
  observed_quantized <- partial_confounding_portable_matrix_sha256(
    matrix, digits = digits
  )
  list(
    observed_object_sha256 = observed_object,
    observed_authoritative_quantized_sha256 = observed_quantized,
    authoritative_hash_match = identical(
      observed_quantized, expected_authoritative_quantized_sha256
    ),
    exact_object_hash_match_diagnostic = identical(
      observed_object, reference_object_sha256
    ),
    exact_object_hash_is_gate = FALSE
  )
}

partial_confounding_v3_reconstruction_gate <- function(reconstruction_hashes) {
  required <- c(
    "object", "hash_role", "portable_hash_schema",
    "expected_quantized_sha256", "observed_quantized_sha256",
    "quantized_hash_match_diagnostic", "reference_platform_object_sha256",
    "observed_platform_object_sha256", "is_cross_platform_gate"
  )
  if (!is.data.frame(reconstruction_hashes) ||
      !all(required %in% names(reconstruction_hashes)) ||
      nrow(reconstruction_hashes) != 4L) {
    return(FALSE)
  }
  expected_objects <- c(
    "raw_log_with_biology", "clean_log_with_biology",
    "raw_intensity", "clean_ground_truth"
  )
  expected_roles <- c(
    "authoritative_cross_platform_gate",
    "authoritative_cross_platform_gate",
    "diagnostic_only", "diagnostic_only"
  )
  gate <- as.logical(reconstruction_hashes$is_cross_platform_gate)
  valid_sha <- function(value) {
    !is.na(value) & grepl("^[0-9a-f]{64}$", value)
  }
  identical(as.character(reconstruction_hashes$object), expected_objects) &&
    identical(as.character(reconstruction_hashes$hash_role), expected_roles) &&
    identical(gate, c(TRUE, TRUE, FALSE, FALSE)) &&
    all(reconstruction_hashes$portable_hash_schema ==
          .partial_confounding_portable_hash_schema) &&
    all(valid_sha(as.character(
      reconstruction_hashes$expected_quantized_sha256
    ))) &&
    all(valid_sha(as.character(
      reconstruction_hashes$observed_quantized_sha256
    ))) &&
    all(valid_sha(as.character(
      reconstruction_hashes$reference_platform_object_sha256
    ))) &&
    all(valid_sha(as.character(
      reconstruction_hashes$observed_platform_object_sha256
    ))) &&
    all(reconstruction_hashes$quantized_hash_match_diagnostic[gate]) &&
    identical(
      as.character(reconstruction_hashes$expected_quantized_sha256[gate]),
      as.character(reconstruction_hashes$observed_quantized_sha256[gate])
    )
}

partial_confounding_portable_hash_policy <- function() {
  list(
    schema = .partial_confounding_portable_hash_schema,
    digits = .partial_confounding_portable_hash_digits,
    quantization_unit =
      .partial_confounding_portable_max_abs_tolerance,
    equivalence_gate = paste(
      "SHA-256 identity of raw_log_with_biology and",
      "clean_log_with_biology after round(x, 9); dimensions and dimnames",
      "are included in R serialization"
    ),
    equivalence_tolerance = paste0(
      "one decimal-quantization cell (",
      format(.partial_confounding_portable_max_abs_tolerance,
             scientific = TRUE), ")"
    ),
    authoritative_objects = c(
      "raw_log_with_biology", "clean_log_with_biology"
    ),
    diagnostic_objects = c("raw_intensity", "clean_ground_truth"),
    intensity_hash_status = paste(
      "intensity object and intensity round9 hashes are recorded as",
      "platform diagnostics and are not reconstruction gates"
    ),
    exact_object_hash_status = "recorded diagnostic only"
  )
}
