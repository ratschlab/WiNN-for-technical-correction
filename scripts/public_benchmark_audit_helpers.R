repo_root_from_script <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg) != 1L) stop("Run this script with Rscript or R --file.")
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)))
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))
}

sha256_file <- function(path) {
  unname(openssl::sha256(file(path)))
}

write_json_pretty <- function(value, path) {
  jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, na = "null")
}

cramers_v <- function(table_value) {
  table_value <- as.matrix(table_value)
  if (sum(table_value) == 0L || min(dim(table_value)) < 2L) return(NA_real_)
  test <- suppressWarnings(stats::chisq.test(table_value, correct = FALSE))
  denom <- sum(table_value) * min(nrow(table_value) - 1L, ncol(table_value) - 1L)
  if (denom <= 0) NA_real_ else sqrt(unname(test$statistic) / denom)
}

standardized_mean_difference <- function(x, group) {
  group <- factor(group)
  if (nlevels(group) != 2L) return(NA_real_)
  x1 <- x[group == levels(group)[1]]
  x2 <- x[group == levels(group)[2]]
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  if (length(x1) < 2L || length(x2) < 2L) return(NA_real_)
  pooled <- sqrt(((length(x1) - 1L) * stats::var(x1) + (length(x2) - 1L) * stats::var(x2)) /
                   (length(x1) + length(x2) - 2L))
  if (!is.finite(pooled) || pooled == 0) return(NA_real_)
  (mean(x2) - mean(x1)) / pooled
}

binary_runs_test <- function(labels) {
  labels <- as.character(labels)
  labels <- labels[!is.na(labels) & nzchar(labels)]
  unique_labels <- unique(labels)
  if (length(unique_labels) != 2L) {
    return(c(n = length(labels), runs = NA_real_, expected = NA_real_, z = NA_real_, p_value = NA_real_))
  }
  n1 <- sum(labels == unique_labels[1])
  n2 <- sum(labels == unique_labels[2])
  n <- n1 + n2
  runs <- 1L + sum(labels[-1] != labels[-length(labels)])
  expected <- 1 + 2 * n1 * n2 / n
  variance <- 2 * n1 * n2 * (2 * n1 * n2 - n) / (n^2 * (n - 1))
  if (!is.finite(variance) || variance <= 0) {
    return(c(n = n, runs = runs, expected = expected, z = NA_real_, p_value = NA_real_))
  }
  z <- (runs - expected) / sqrt(variance)
  c(n = n, runs = runs, expected = expected, z = z, p_value = 2 * stats::pnorm(-abs(z)))
}

rank_biserial_order <- function(binary_group, order_value) {
  binary_group <- as.character(binary_group)
  ok <- !is.na(binary_group) & is.finite(order_value)
  binary_group <- binary_group[ok]
  order_value <- order_value[ok]
  levels_found <- sort(unique(binary_group))
  if (length(levels_found) != 2L) return(NA_real_)
  x0 <- order_value[binary_group == levels_found[1]]
  x1 <- order_value[binary_group == levels_found[2]]
  if (!length(x0) || !length(x1)) return(NA_real_)
  u <- as.numeric(stats::wilcox.test(x1, x0, exact = FALSE)$statistic)
  2 * u / (length(x0) * length(x1)) - 1
}

intensity_audit <- function(mat) {
  values <- as.numeric(mat)
  finite_values <- values[is.finite(values)]
  quantile_probs <- c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 1)
  list(
    dimensions = c(features = nrow(mat), injections = ncol(mat)),
    minimum = min(finite_values),
    maximum = max(finite_values),
    median = stats::median(finite_values),
    quantiles = stats::setNames(as.list(stats::quantile(finite_values, quantile_probs, names = FALSE)), as.character(quantile_probs)),
    zero_count = sum(finite_values == 0),
    zero_fraction = mean(finite_values == 0),
    nonfinite_count = sum(!is.finite(values)),
    nonfinite_fraction = mean(!is.finite(values)),
    integer_like_fraction = mean(abs(finite_values - round(finite_values)) < 1e-8),
    zero_variance_features = sum(apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) == 0), na.rm = TRUE),
    maximum_feature_missingness = max(rowMeans(!is.finite(mat))),
    maximum_sample_missingness = max(colMeans(!is.finite(mat)))
  )
}

histogram_source <- function(mat, dataset) {
  values <- as.numeric(mat)
  values <- values[is.finite(values)]
  make_hist <- function(x, scale_label) {
    h <- graphics::hist(x, breaks = 80L, plot = FALSE)
    data.frame(
      dataset = dataset,
      scale = scale_label,
      bin_left = head(h$breaks, -1L),
      bin_right = tail(h$breaks, -1L),
      bin_mid = h$mids,
      count = h$counts,
      density = h$density,
      stringsAsFactors = FALSE
    )
  }
  rbind(make_hist(values, "intensity"), make_hist(log1p(values), "log1p_intensity"))
}

qc_spacing_by_batch <- function(metadata) {
  batches <- unique(as.character(metadata$batch))
  do.call(rbind, lapply(batches, function(batch_id) {
    d <- metadata[metadata$batch == batch_id, , drop = FALSE]
    d <- d[order(d$within_batch_order), , drop = FALSE]
    qc_positions <- d$within_batch_order[d$is_qc]
    gaps <- diff(c(0L, qc_positions, nrow(d) + 1L)) - 1L
    data.frame(
      batch = batch_id,
      total_injections = nrow(d),
      study_injections = sum(!d$is_qc),
      qc_injections = sum(d$is_qc),
      first_qc_position = min(qc_positions),
      last_qc_position = max(qc_positions),
      median_study_between_qcs = stats::median(gaps),
      maximum_study_between_qcs = max(gaps),
      stringsAsFactors = FALSE
    )
  }))
}

write_suitability_csv <- function(criteria, path) {
  stopifnot(all(c("criterion_id", "criterion", "status", "evidence") %in% names(criteria)))
  write.csv(criteria, path, row.names = FALSE, quote = TRUE, na = "")
}
