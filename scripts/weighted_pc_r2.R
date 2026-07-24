# Explicit weighted-PC R-squared metric used by benchmark evaluations.
#
# The target type is deliberately mandatory. Numeric category codes must never
# be interpreted as a one-degree-of-freedom continuous trend.
weighted_pc_r2_explicit <- function(mat, target, target_type = NULL, n_pcs = 5L) {
  if (is.null(target_type) || length(target_type) != 1L) {
    stop("target_type must be declared as either 'categorical' or 'continuous'.")
  }
  target_type <- match.arg(target_type, c("categorical", "continuous"))

  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  if (length(target) != ncol(mat)) {
    stop("target length (", length(target), ") does not equal the number of matrix columns (", ncol(mat), ").")
  }
  if (length(n_pcs) != 1L || !is.finite(n_pcs) || n_pcs < 1L) {
    stop("n_pcs must be a positive finite scalar.")
  }

  if (target_type == "categorical") {
    observed <- !is.na(target)
    target <- droplevels(factor(target[observed]))
    if (nlevels(target) < 2L) return(NA_real_)
  } else {
    if (!is.numeric(target)) {
      stop("A continuous target must be numeric; explicit coercion is not performed.")
    }
    observed <- is.finite(target)
    target <- as.numeric(target[observed])
    if (length(target) < 2L || !is.finite(stats::sd(target)) || stats::sd(target) <= 0) {
      return(NA_real_)
    }
  }

  mat <- mat[, observed, drop = FALSE]
  if (ncol(mat) < 3L) return(NA_real_)
  if (any(mat < -1, na.rm = TRUE)) {
    stop("weighted_pc_r2_explicit requires intensities >= -1 for log1p transformation.")
  }

  z <- t(log1p(mat))
  keep <- apply(z, 2L, function(v) {
    all(is.finite(v)) && is.finite(stats::sd(v)) && stats::sd(v) > 0
  })
  z <- z[, keep, drop = FALSE]
  if (nrow(z) < 3L || ncol(z) < 2L) return(NA_real_)

  pca <- tryCatch(
    stats::prcomp(z, center = TRUE, scale. = TRUE),
    error = function(e) NULL
  )
  if (is.null(pca) || !length(pca$sdev) || sum(pca$sdev^2) <= 0) return(NA_real_)

  pcs <- min(as.integer(n_pcs), ncol(pca$x), length(pca$sdev))
  variance_weights <- (pca$sdev^2 / sum(pca$sdev^2))[seq_len(pcs)]
  r2 <- vapply(seq_len(pcs), function(i) {
    fit <- tryCatch(stats::lm(pca$x[, i] ~ target), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    value <- summary(fit)$r.squared
    if (length(value) == 1L && is.finite(value)) value else NA_real_
  }, numeric(1))

  valid <- is.finite(r2) & is.finite(variance_weights) & variance_weights > 0
  if (!any(valid)) return(NA_real_)
  sum(r2[valid] * variance_weights[valid]) / sum(variance_weights[valid])
}
