build_contiguous_run_segments <- function(meta_df, sample_id_col = "sample_id", order_col = "run_order", batch_col = "batch") {
  md <- meta_df[order(meta_df[[order_col]], na.last = TRUE), , drop = FALSE]
  if (nrow(md) == 0L) {
    md$segment_index <- integer(0)
    md$segment_id <- character(0)
    return(md)
  }

  batch_vals <- as.character(md[[batch_col]])
  segment_index <- cumsum(c(TRUE, batch_vals[-1] != batch_vals[-length(batch_vals)]))

  md$segment_index <- segment_index
  md$segment_id <- paste0("batch_", batch_vals, "_segment_", sprintf("%03d", segment_index))
  md
}

fit_run_order_gam_explained <- function(series, run_order, min_obs = 6L, k_max = 6L) {
  ok <- is.finite(series) & is.finite(run_order)
  n_obs <- sum(ok)

  if (n_obs < min_obs) {
    return(c(explained = NA_real_, n_obs = n_obs, k = NA_real_))
  }

  y <- as.numeric(series[ok])
  x <- as.numeric(run_order[ok])
  n_unique_x <- length(unique(x))

  if (n_unique_x < 4L) {
    return(c(explained = NA_real_, n_obs = n_obs, k = NA_real_))
  }

  if (stats::sd(y) <= .Machine$double.eps) {
    return(c(explained = 0, n_obs = n_obs, k = NA_real_))
  }

  k_use <- min(as.integer(k_max), n_unique_x - 1L, max(3L, floor(n_obs / 2L)))
  if (!is.finite(k_use) || k_use < 3L) {
    return(c(explained = NA_real_, n_obs = n_obs, k = NA_real_))
  }

  fit <- tryCatch(
    mgcv::gam(y ~ s(x, bs = "cr", k = k_use), method = "REML"),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(c(explained = NA_real_, n_obs = n_obs, k = as.numeric(k_use)))
  }

  fit_summary <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(fit_summary)) {
    return(c(explained = NA_real_, n_obs = n_obs, k = as.numeric(k_use)))
  }

  explained <- if (is.finite(fit_summary$dev.expl)) fit_summary$dev.expl else fit_summary$r.sq
  explained <- if (is.finite(explained)) min(max(as.numeric(explained), 0), 1) else NA_real_

  c(explained = explained, n_obs = n_obs, k = as.numeric(k_use))
}

compute_metabolite_segment_gam <- function(
  mat,
  meta_df,
  method_label,
  sample_id_col = "sample_id",
  order_col = "run_order",
  batch_col = "batch",
  transform_fun = log1p,
  min_obs = 6L,
  k_max = 6L
) {
  md <- meta_df[meta_df[[sample_id_col]] %in% colnames(mat), , drop = FALSE]
  md <- build_contiguous_run_segments(md, sample_id_col = sample_id_col, order_col = order_col, batch_col = batch_col)

  if (nrow(md) == 0L || nrow(mat) == 0L) {
    return(tibble::tibble(
      method = character(),
      metabolite = character(),
      batch = character(),
      segment_id = character(),
      n_obs = integer(),
      k = numeric(),
      explained = numeric()
    ))
  }

  seg_ids <- unique(md$segment_id)
  transformed <- transform_fun(mat[, md[[sample_id_col]], drop = FALSE])

  dplyr::bind_rows(lapply(seg_ids, function(seg_id) {
    seg_md <- md[md$segment_id == seg_id, , drop = FALSE]
    seg_sample_ids <- seg_md[[sample_id_col]]
    seg_mat <- transformed[, seg_sample_ids, drop = FALSE]

    gam_stats <- t(vapply(
      seq_len(nrow(seg_mat)),
      function(i) fit_run_order_gam_explained(seg_mat[i, ], seg_md[[order_col]], min_obs = min_obs, k_max = k_max),
      numeric(3)
    ))

    tibble::tibble(
      method = method_label,
      metabolite = rownames(seg_mat),
      batch = as.character(seg_md[[batch_col]][1]),
      segment_id = seg_id,
      n_obs = as.integer(gam_stats[, "n_obs"]),
      k = as.numeric(gam_stats[, "k"]),
      explained = as.numeric(gam_stats[, "explained"])
    )
  }))
}

summarise_metabolite_segment_gam <- function(df, digits = 4) {
  df |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      n_pairs = sum(is.finite(explained)),
      mean_explained = mean(explained, na.rm = TRUE),
      median_explained = stats::median(explained, na.rm = TRUE),
      sd_explained = stats::sd(explained, na.rm = TRUE),
      mean_sd = ifelse(
        is.finite(mean_explained),
        sprintf(paste0("%.", digits, "f \u00b1 %.", digits, "f"), mean_explained, ifelse(is.finite(sd_explained), sd_explained, NA_real_)),
        NA_character_
      ),
      .groups = "drop"
    )
}

ljung_box_series_metrics <- function(series, max_lag = 10L, min_obs = 8L) {
  x <- as.numeric(series)
  x <- x[is.finite(x)]
  n_obs <- length(x)
  lag_use <- min(as.integer(max_lag), max(1L, floor(n_obs / 3L)), n_obs - 1L)

  if (n_obs < min_obs || lag_use < 1L || isTRUE(stats::sd(x, na.rm = TRUE) < .Machine$double.eps)) {
    return(c(lb_stat = NA_real_, p_value = NA_real_, lag = as.numeric(lag_use), n_obs = as.numeric(n_obs)))
  }

  lb <- tryCatch(
    stats::Box.test(x, lag = lag_use, type = "Ljung-Box"),
    error = function(e) NULL
  )

  if (is.null(lb)) {
    return(c(lb_stat = NA_real_, p_value = NA_real_, lag = as.numeric(lag_use), n_obs = as.numeric(n_obs)))
  }

  c(
    lb_stat = as.numeric(lb$statistic),
    p_value = as.numeric(lb$p.value),
    lag = as.numeric(lag_use),
    n_obs = as.numeric(n_obs)
  )
}

compute_metabolite_segment_ljung_box <- function(
  mat,
  meta_df,
  method_label,
  sample_id_col = "sample_id",
  order_col = "run_order",
  batch_col = "batch",
  transform_fun = log1p,
  min_obs = 8L,
  max_lag = 10L
) {
  md <- meta_df[meta_df[[sample_id_col]] %in% colnames(mat), , drop = FALSE]
  md <- build_contiguous_run_segments(md, sample_id_col = sample_id_col, order_col = order_col, batch_col = batch_col)

  if (nrow(md) == 0L || nrow(mat) == 0L) {
    return(tibble::tibble(
      method = character(),
      metabolite = character(),
      batch = character(),
      segment_id = character(),
      n_obs = integer(),
      lag = numeric(),
      lb_stat = numeric(),
      p_value = numeric()
    ))
  }

  seg_ids <- unique(md$segment_id)
  transformed <- transform_fun(mat[, md[[sample_id_col]], drop = FALSE])

  dplyr::bind_rows(lapply(seg_ids, function(seg_id) {
    seg_md <- md[md$segment_id == seg_id, , drop = FALSE]
    seg_sample_ids <- seg_md[[sample_id_col]]
    seg_mat <- transformed[, seg_sample_ids, drop = FALSE]

    lb_stats <- t(vapply(
      seq_len(nrow(seg_mat)),
      function(i) ljung_box_series_metrics(seg_mat[i, ], max_lag = max_lag, min_obs = min_obs),
      numeric(4)
    ))

    tibble::tibble(
      method = method_label,
      metabolite = rownames(seg_mat),
      batch = as.character(seg_md[[batch_col]][1]),
      segment_id = seg_id,
      n_obs = as.integer(lb_stats[, "n_obs"]),
      lag = as.numeric(lb_stats[, "lag"]),
      lb_stat = as.numeric(lb_stats[, "lb_stat"]),
      p_value = as.numeric(lb_stats[, "p_value"])
    )
  }))
}

annotate_metabolite_segment_ljung_box <- function(df, alpha = 0.05, adjust_method = "BH") {
  df |>
    dplyr::group_by(method) |>
    dplyr::mutate(
      p_adj = stats::p.adjust(p_value, method = adjust_method),
      is_autocorrelated = is.finite(p_adj) & (p_adj < alpha)
    ) |>
    dplyr::ungroup()
}

summarise_metabolite_segment_ljung_box <- function(df) {
  df |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      n_tested = sum(is.finite(p_value)),
      n_autocorrelated = sum(is_autocorrelated, na.rm = TRUE),
      prop_autocorrelated = ifelse(n_tested > 0, n_autocorrelated / n_tested, NA_real_),
      median_lb_stat = stats::median(lb_stat, na.rm = TRUE),
      .groups = "drop"
    )
}
