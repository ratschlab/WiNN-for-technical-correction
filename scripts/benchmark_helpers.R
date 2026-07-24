# Shared plotting, export, correction, and metric helpers for the public
# benchmarks. Dataset-specific input adapters deliberately live elsewhere.

theme_publication <- function() {
  ggplot2::theme_minimal(base_size = 12, base_family = "Helvetica") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.4, linetype = "dotted"),
      axis.title = ggplot2::element_text(color = "#111827", face = "bold", size = 12),
      axis.text = ggplot2::element_text(color = "#374151", size = 10),
      axis.ticks = ggplot2::element_line(color = "#9CA3AF"),
      strip.text = ggplot2::element_text(face = "bold", color = "#111827", size = 11),
      strip.background = ggplot2::element_rect(fill = "#F3F4F6", color = NA),
      plot.title = ggplot2::element_text(face = "bold", size = 13, color = "#111827", margin = ggplot2::margin(b = 5), lineheight = 1.05),
      plot.subtitle = ggplot2::element_text(size = 10.2, color = "#4B5563", margin = ggplot2::margin(b = 10), lineheight = 1.05),
      legend.title = ggplot2::element_text(face = "bold", color = "#111827"),
      legend.text = ggplot2::element_text(color = "#374151"),
      legend.position = "top",
      text = ggplot2::element_text(color = "#111827")
    )
}

save_pdf <- function(plot, plot_dir, filename, width, height) {
  ggplot2::ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot,
    device = "pdf",
    width = width,
    height = height,
    units = "in"
  )
}

save_csv <- function(df, dir_path, filename) {
  write.csv(df, file = file.path(dir_path, filename), row.names = FALSE, quote = TRUE)
}

write_value_exports <- function(df, dir_path, full_filename, summary_filename, group_cols, value_col = "value", digits = 4, write_full = TRUE) {
  if (isTRUE(write_full)) {
    save_csv(df, dir_path, full_filename)
  }

  summary_df <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n = sum(is.finite(.data[[value_col]])),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = stats::sd(.data[[value_col]], na.rm = TRUE),
      mean_sd = ifelse(
        is.finite(mean),
        sprintf(paste0("%.", digits, "f \u00b1 %.", digits, "f"), mean, ifelse(is.finite(sd), sd, NA_real_)),
        NA_character_
      ),
      .groups = "drop"
    )

  save_csv(summary_df, dir_path, summary_filename)
  summary_df
}

format_value_label <- function(x, digits = 2, integer_threshold = 1e-10) {
  ifelse(abs(x - round(x)) < integer_threshold, sprintf("%.0f", round(x)), sprintf(paste0("%.", digits, "f"), x))
}

format_mean_label <- function(mean_x, digits = 2) {
  ifelse(
    is.finite(mean_x),
    sprintf(paste0("%.", digits, "f"), mean_x),
    NA_character_
  )
}

export_distribution_suite <- function(
  df,
  plot_dir,
  csv_dir,
  stem,
  title,
  subtitle,
  y_label,
  method_palette,
  facet_var = NULL,
  ecdf_subtitle = NULL,
  violin_width = 14,
  violin_height = 6,
  ecdf_width = 11,
  ecdf_height = 6,
  label_digits = 2
) {
  df <- df |>
    dplyr::filter(is.finite(value))

  save_csv(df, csv_dir, paste0(stem, "_full_values.csv"))
  write_value_exports(
    df,
    dir_path = csv_dir,
    full_filename = paste0(stem, "_long.csv"),
    summary_filename = paste0(stem, "_summary_mean_sd.csv"),
    group_cols = c("method", if (!is.null(facet_var)) facet_var else character(0)),
    value_col = "value"
  )

  summary_df <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c("method", if (!is.null(facet_var)) facet_var else character(0))))) |>
    dplyr::summarise(
      n = sum(is.finite(value)),
      mean_value = mean(value, na.rm = TRUE),
      sd_value = stats::sd(value, na.rm = TRUE),
      se_value = sd_value / sqrt(pmax(n, 1)),
      mean_label = format_mean_label(mean_value, digits = label_digits),
      .groups = "drop"
    )

  save_csv(summary_df, csv_dir, paste0(stem, "_bar_summary_mean_se.csv"))

  p_violin <- ggplot2::ggplot(df, ggplot2::aes(method, value, fill = method)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.28, color = NA) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.65, color = "#111827") +
    ggplot2::scale_fill_manual(values = method_palette, guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.01), expand = ggplot2::expansion(mult = c(0.02, 0.06))) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = y_label) +
    theme_publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 28, hjust = 1))

  if (!is.null(facet_var)) {
    facet_formula <- stats::as.formula(paste("~", facet_var))
    p_violin <- p_violin + ggplot2::facet_wrap(facet_formula, scales = "free_y")
  }

  save_pdf(p_violin, plot_dir, paste0(stem, "_distribution.pdf"), width = violin_width, height = violin_height)

  y_expand_mult <- if (is.null(facet_var)) c(0, 0.16) else c(0, 0.2)
  p_bar <- ggplot2::ggplot(summary_df, ggplot2::aes(method, mean_value, fill = method)) +
    ggplot2::geom_col(width = 0.72, color = "#1F2937", linewidth = 0.25) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
      width = 0.16,
      linewidth = 0.5,
      color = "#111827"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = mean_label, y = mean_value + se_value),
      vjust = -0.28,
      size = 2.7,
      color = "#111827"
    ) +
    ggplot2::scale_fill_manual(values = method_palette, guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.001), expand = ggplot2::expansion(mult = y_expand_mult)) +
    ggplot2::labs(
      title = paste0(title, " (mean with SE bars)"),
      subtitle = subtitle,
      x = NULL,
      y = paste0("Mean ", y_label)
    ) +
    theme_publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 28, hjust = 1))

  if (!is.null(facet_var)) {
    facet_formula <- stats::as.formula(paste("~", facet_var))
    p_bar <- p_bar + ggplot2::facet_wrap(facet_formula, scales = "free_y")
  }

  save_pdf(p_bar, plot_dir, paste0(stem, "_bar_mean_se.pdf"), width = violin_width, height = violin_height)

  p_ecdf <- ggplot2::ggplot(df, ggplot2::aes(value, color = method)) +
    ggplot2::stat_ecdf(linewidth = 0.95) +
    ggplot2::scale_color_manual(values = method_palette, name = "Method") +
    ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.01), expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(
      title = paste0(title, " ECDF"),
      subtitle = ifelse(
        is.null(ecdf_subtitle),
        "Right-shifted curves indicate stronger performance for agreement metrics and lower curves indicate stronger performance for CV metrics.",
        ecdf_subtitle
      ),
      x = y_label,
      y = "ECDF"
    ) +
    theme_publication()

  if (!is.null(facet_var)) {
    facet_formula <- stats::as.formula(paste("~", facet_var))
    p_ecdf <- p_ecdf + ggplot2::facet_wrap(facet_formula, scales = "free_x")
  }

  save_pdf(p_ecdf, plot_dir, paste0(stem, "_ecdf.pdf"), width = ecdf_width, height = ecdf_height)

  list(violin = p_violin, bar = p_bar, ecdf = p_ecdf, summary = summary_df)
}

compute_delta_vs_raw_long <- function(df, id_col = "feature", raw_method = "Raw") {
  stopifnot(all(c("method", id_col, "agreement_stat", "value") %in% names(df)))

  method_levels <- if (is.factor(df$method)) levels(df$method) else unique(as.character(df$method))
  agreement_levels <- if (is.factor(df$agreement_stat)) levels(df$agreement_stat) else unique(as.character(df$agreement_stat))

  raw_df <- df |>
    dplyr::mutate(method = as.character(method), agreement_stat = as.character(agreement_stat)) |>
    dplyr::filter(method == raw_method) |>
    dplyr::select(dplyr::all_of(c(id_col, "agreement_stat")), raw_value = value)

  if (nrow(raw_df) == 0L) {
    stop("Raw method rows were not found while computing delta vs raw.")
  }

  df |>
    dplyr::mutate(method = as.character(method), agreement_stat = as.character(agreement_stat)) |>
    dplyr::left_join(raw_df, by = c(id_col, "agreement_stat")) |>
    dplyr::filter(is.finite(value), is.finite(raw_value)) |>
    dplyr::mutate(
      corrected_value = value,
      value = corrected_value - raw_value,
      method = factor(method, levels = method_levels),
      agreement_stat = factor(agreement_stat, levels = agreement_levels)
    )
}

clip_nonnegative <- function(m) {
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  m[is.infinite(m)] <- NA_real_
  pmax(m, 0)
}

subset_to_features <- function(mat, feature_ids) {
  feature_ids <- intersect(feature_ids, rownames(mat))
  mat[feature_ids, , drop = FALSE]
}

safe_cor <- function(a, b, method = "pearson", min_n = 6L) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < min_n) {
    return(NA_real_)
  }
  if (stats::sd(a[ok]) <= 0 || stats::sd(b[ok]) <= 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(a[ok], b[ok], method = method))
}

calc_icc_a1_matrix <- function(ratings) {
  ratings <- as.matrix(ratings)
  ratings <- ratings[stats::complete.cases(ratings), , drop = FALSE]

  n <- nrow(ratings)
  k <- ncol(ratings)
  if (n < 2L || k < 2L) {
    return(NA_real_)
  }

  grand_mean <- mean(ratings)
  row_means <- rowMeans(ratings)
  col_means <- colMeans(ratings)

  ss_rows <- k * sum((row_means - grand_mean)^2)
  ss_cols <- n * sum((col_means - grand_mean)^2)
  ss_total <- sum((ratings - grand_mean)^2)
  ss_error <- ss_total - ss_rows - ss_cols

  ms_rows <- ss_rows / (n - 1L)
  ms_cols <- ss_cols / (k - 1L)
  ms_error <- ss_error / ((n - 1L) * (k - 1L))

  denom <- ms_rows + (k - 1L) * ms_error + (k * (ms_cols - ms_error) / n)
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }

  (ms_rows - ms_error) / denom
}

align_meta_to_mat <- function(x_mat, meta_df) {
  sample_ids <- colnames(x_mat)
  md <- meta_df[match(sample_ids, meta_df$sample_id), , drop = FALSE]

  if (anyNA(md$sample_id)) {
    stop("Sample IDs in matrix do not match metadata sample IDs.")
  }

  md
}

run_silently <- function(code) {
  suppressWarnings(suppressMessages(force(code)))
}

run_method <- function(label, fn, tuning_sec = 0) {
  process_sec <- system.time({
    out_raw <- fn()
  })[["elapsed"]]

  retained_attrs <- attributes(out_raw)[intersect(
    names(attributes(out_raw)),
    c(
      "dip_pair_protocol", "negative_values_clipped", "detected_batch",
      "detected_batch_units", "pelt_penalty",
      "qcrfsc_uncorrected_feature_ids", "qcrfsc_uncorrected_sample_ids",
      "qcrfsc_uncorrected_feature_count", "qcrfsc_uncorrected_sample_count"
    )
  )]
  out <- clip_nonnegative(out_raw)
  for (attr_name in names(retained_attrs)) attr(out, attr_name) <- retained_attrs[[attr_name]]

  tuning_sec <- as.numeric(tuning_sec)[1]
  if (!is.finite(tuning_sec)) {
    tuning_sec <- 0
  }

  list(
    label = label,
    data = out,
    tuning_sec = tuning_sec,
    process_sec = process_sec,
    elapsed_sec = tuning_sec + process_sec
  )
}

run_method_cached <- function(label, fn, plot_dir, tuning_sec = 0, force = FALSE, cache_version = "") {
  cache_dir <- file.path(plot_dir, "method_cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  safe_label <- gsub("[^A-Za-z0-9]+", "_", tolower(label))
  cache_suffix <- ifelse(nzchar(cache_version), paste0("_", cache_version), "")
  cache_file <- file.path(cache_dir, paste0(safe_label, cache_suffix, ".rds"))

  if (!isTRUE(force) && file.exists(cache_file)) {
    out <- readRDS(cache_file)
    out$tuning_sec <- as.numeric(tuning_sec)[1]
    if (!is.finite(out$tuning_sec)) out$tuning_sec <- 0
    if (is.null(out$process_sec) || !is.finite(out$process_sec)) {
      out$process_sec <- out$elapsed_sec
    }
    out$elapsed_sec <- out$process_sec + out$tuning_sec
    return(out)
  }

  out <- run_method(label = label, fn = fn, tuning_sec = tuning_sec)
  saveRDS(out, cache_file)
  out
}

run_matrix_cached <- function(label, fn, plot_dir, force = FALSE, cache_tag = "") {
  cache_dir <- file.path(plot_dir, "matrix_cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  safe_label <- gsub("[^A-Za-z0-9]+", "_", tolower(label))
  tag <- ifelse(nzchar(cache_tag), paste0("_", cache_tag), "")
  cache_file <- file.path(cache_dir, paste0(safe_label, tag, ".rds"))

  if (!isTRUE(force) && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  out <- clip_nonnegative(fn())
  saveRDS(out, cache_file)
  out
}

with_no_cluster <- function(code, patch_do_parallel = FALSE) {
  parallel_ns <- asNamespace("parallel")
  original_parallel <- list(
    makeCluster = get("makeCluster", envir = parallel_ns),
    clusterExport = get("clusterExport", envir = parallel_ns),
    stopCluster = get("stopCluster", envir = parallel_ns)
  )

  if (bindingIsLocked("makeCluster", parallel_ns)) unlockBinding("makeCluster", parallel_ns)
  if (bindingIsLocked("clusterExport", parallel_ns)) unlockBinding("clusterExport", parallel_ns)
  if (bindingIsLocked("stopCluster", parallel_ns)) unlockBinding("stopCluster", parallel_ns)

  assign("makeCluster", function(spec, ...) NULL, envir = parallel_ns)
  assign("clusterExport", function(cl, ...) invisible(NULL), envir = parallel_ns)
  assign("stopCluster", function(cl, ...) invisible(NULL), envir = parallel_ns)

  lockBinding("makeCluster", parallel_ns)
  lockBinding("clusterExport", parallel_ns)
  lockBinding("stopCluster", parallel_ns)

  do_parallel_ns <- NULL
  original_register <- NULL

  if (isTRUE(patch_do_parallel) && requireNamespace("doParallel", quietly = TRUE)) {
    do_parallel_ns <- asNamespace("doParallel")
    original_register <- get("registerDoParallel", envir = do_parallel_ns)
    if (bindingIsLocked("registerDoParallel", do_parallel_ns)) unlockBinding("registerDoParallel", do_parallel_ns)
    assign("registerDoParallel", function(cl, ...) foreach::registerDoSEQ(), envir = do_parallel_ns)
    lockBinding("registerDoParallel", do_parallel_ns)
  }

  on.exit({
    if (!is.null(do_parallel_ns) && !is.null(original_register)) {
      if (bindingIsLocked("registerDoParallel", do_parallel_ns)) unlockBinding("registerDoParallel", do_parallel_ns)
      assign("registerDoParallel", original_register, envir = do_parallel_ns)
      lockBinding("registerDoParallel", do_parallel_ns)
    }

    if (bindingIsLocked("makeCluster", parallel_ns)) unlockBinding("makeCluster", parallel_ns)
    if (bindingIsLocked("clusterExport", parallel_ns)) unlockBinding("clusterExport", parallel_ns)
    if (bindingIsLocked("stopCluster", parallel_ns)) unlockBinding("stopCluster", parallel_ns)

    assign("makeCluster", original_parallel$makeCluster, envir = parallel_ns)
    assign("clusterExport", original_parallel$clusterExport, envir = parallel_ns)
    assign("stopCluster", original_parallel$stopCluster, envir = parallel_ns)

    lockBinding("makeCluster", parallel_ns)
    lockBinding("clusterExport", parallel_ns)
    lockBinding("stopCluster", parallel_ns)
  }, add = TRUE)

  force(code)
}

make_malbac_object <- function(x_mat, meta_df, jitter_eps = 0) {
  x_mat <- as.matrix(x_mat)
  if (jitter_eps > 0) {
    x_mat <- x_mat + matrix(stats::runif(length(x_mat), min = 0, max = jitter_eps), nrow = nrow(x_mat))
  }

  safe_feature <- paste0("Feature_", seq_len(nrow(x_mat)))
  safe_sample <- paste0("Sample_", seq_len(ncol(x_mat)))
  feature_map <- stats::setNames(rownames(x_mat), safe_feature)
  sample_map <- stats::setNames(meta_df$sample_id, safe_sample)

  x_safe <- x_mat
  rownames(x_safe) <- safe_feature
  colnames(x_safe) <- safe_sample

  f_data <- data.frame(
    sampleID = safe_sample,
    sample_type = meta_df$class,
    batch = meta_df$batch,
    run_order = meta_df$run_order,
    stringsAsFactors = FALSE
  )
  rownames(f_data) <- f_data$sampleID

  e_data <- data.frame(
    feature = safe_feature,
    x_safe,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  obj <- pmartR::as.metabData(
    e_data = e_data,
    f_data = f_data,
    edata_cname = "feature",
    fdata_cname = "sampleID"
  )

  obj <- pmartR::group_designation(
    omicsData = obj,
    main_effects = "sample_type",
    batch_id = "batch"
  )

  list(obj = obj, sample_map = sample_map, feature_map = feature_map)
}

extract_malbac_matrix <- function(malbac_obj, sample_map, feature_map = NULL) {
  e_data <- malbac_obj$e_data
  m <- as.matrix(e_data[, -1, drop = FALSE])
  feature_ids <- as.character(e_data[, 1])
  if (is.null(feature_map)) {
    rownames(m) <- feature_ids
  } else {
    mapped_features <- unname(feature_map[feature_ids])
    mapped_features[is.na(mapped_features)] <- feature_ids[is.na(mapped_features)]
    rownames(m) <- mapped_features
  }

  mapped <- sample_map[colnames(m)]
  colnames(m) <- unname(mapped)
  m
}

run_combat <- function(x_mat, meta_df, par_prior = TRUE) {
  md <- align_meta_to_mat(x_mat, meta_df)
  x_log <- log1p(x_mat)
  y_log <- sva::ComBat(
    dat = x_log,
    batch = md$batch,
    mod = NULL,
    par.prior = par_prior,
    prior.plots = FALSE
  )
  expm1(y_log)
}

run_qc_rlsc_with_controls <- function(x_mat, control_ids, meta_df, span = 0.75) {
  md <- align_meta_to_mat(x_mat, meta_df) |>
    dplyr::arrange(run_order)
  x_ord <- x_mat[, md$sample_id, drop = FALSE]
  dat <- as.data.frame(t(x_ord), check.names = FALSE)
  y <- ifelse(md$sample_id %in% control_ids, "QC", "Sample")
  corrected <- qcrlscR::qc.rlsc.wrap(
    dat = dat,
    cls.qc = y,
    cls.bl = factor(md$batch),
    method = "divide",
    intra = TRUE,
    opti = FALSE,
    log10 = FALSE,
    outl = FALSE,
    shift = TRUE,
    span = span
  )
  corrected <- t(as.matrix(corrected))
  rownames(corrected) <- rownames(x_ord)
  colnames(corrected) <- md$sample_id
  corrected[, meta_df$sample_id, drop = FALSE]
}

run_qc_rfsc_with_controls <- function(x_mat, control_ids, meta_df, ntree = 500, coCV = 30, Frule = 0.8) {
  md <- align_meta_to_mat(x_mat, meta_df) |>
    dplyr::arrange(run_order)
  x_ord <- x_mat[, md$sample_id, drop = FALSE]

  control_ids <- intersect(control_ids, md$sample_id)
  if (length(control_ids) < 2L) {
    stop("QC-RFSC requires at least 2 QC controls.")
  }

  first_id <- md$sample_id[1]
  last_id <- md$sample_id[nrow(md)]
  if (!(first_id %in% control_ids) || !(last_id %in% control_ids)) {
    first_qc <- control_ids[which.min(md$run_order[match(control_ids, md$sample_id)])]
    last_qc <- control_ids[which.max(md$run_order[match(control_ids, md$sample_id)])]
    reordered_ids <- c(first_qc, setdiff(md$sample_id, c(first_qc, last_qc)), last_qc)
    md <- md[match(reordered_ids, md$sample_id), , drop = FALSE]
    x_ord <- x_ord[, reordered_ids, drop = FALSE]
  }

  qc_flag <- md$sample_id %in% control_ids
  stat_names <- ifelse(qc_flag, paste0(md$sample_id, "_QC"), md$sample_id)
  safe_feature_ids <- paste0("Feature_", seq_len(nrow(x_ord)))
  feature_map <- stats::setNames(rownames(x_ord), safe_feature_ids)

  samPeno <- data.frame(
    sample = stat_names,
    batch = md$batch,
    class = ifelse(qc_flag, NA_integer_, 1L),
    order = seq_len(nrow(md)),
    stringsAsFactors = FALSE
  )

  if (!is.na(samPeno$class[1]) || !is.na(samPeno$class[nrow(samPeno)])) {
    stop("QC-RFSC setup failed to place QC at first and last positions.")
  }

  samFile <- data.frame(name = safe_feature_ids, x_ord, check.names = FALSE)
  colnames(samFile)[-1] <- stat_names

  tmp <- tempfile("shiftcor_")
  dir.create(tmp)
  old <- getwd()
  setwd(tmp)
  on.exit(setwd(old), add = TRUE)

  write.csv(samPeno, "samPeno.csv", row.names = FALSE, quote = FALSE)
  write.csv(samFile, "samFile.csv", row.names = FALSE, quote = FALSE)

  run_silently(
    statTarget::shiftCor(
      samPeno = "samPeno.csv",
      samFile = "samFile.csv",
      Frule = Frule,
      MLmethod = "QCRFSC",
      ntree = ntree,
      coCV = coCV,
      plot = FALSE
    )
  )

  out_file <- file.path("statTarget", "shiftCor", "After_shiftCor", "shift_all_cor.csv")
  if (!file.exists(out_file)) stop("statTarget did not produce shift_all_cor.csv")

  corrected_df <- read.csv(out_file, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
  if ("class" %in% rownames(corrected_df)) {
    corrected_df <- corrected_df[rownames(corrected_df) != "class", , drop = FALSE]
  }

  corrected <- apply(corrected_df, 2, as.numeric)
  rownames(corrected) <- unname(feature_map[rownames(corrected_df)])
  colnames(corrected) <- sub("_QC$", "", colnames(corrected_df))
  corrected <- corrected / 1000
  missing_features <- setdiff(rownames(x_mat), rownames(corrected))
  missing_samples <- setdiff(colnames(x_mat), colnames(corrected))
  full <- x_mat
  common_features <- intersect(rownames(x_mat), rownames(corrected))
  common_samples <- intersect(colnames(x_mat), colnames(corrected))
  full[common_features, common_samples] <- corrected[
    common_features, common_samples, drop = FALSE
  ]
  attr(full, "qcrfsc_uncorrected_feature_ids") <- missing_features
  attr(full, "qcrfsc_uncorrected_sample_ids") <- missing_samples
  attr(full, "qcrfsc_uncorrected_feature_count") <- length(missing_features)
  attr(full, "qcrfsc_uncorrected_sample_count") <- length(missing_samples)
  full
}

run_tiger_with_holdout_qc <- function(
  x_mat,
  train_qc_ids,
  test_qc_ids,
  meta_df,
  use_injection = TRUE,
  mtry_percent = 0.4,
  nodesize_percent = 0.4,
  ntree = 5,
  parallel_cores = 1,
  duplicate_test_qc = TRUE,
  include_nonqc_test = TRUE
) {
  md <- align_meta_to_mat(x_mat, meta_df)
  safe_feature_ids <- paste0("Feature_", seq_len(nrow(x_mat)))
  feature_map <- stats::setNames(rownames(x_mat), safe_feature_ids)
  x_safe <- x_mat
  rownames(x_safe) <- safe_feature_ids

  sample_df <- tibble::tibble(
    sampleID = md$sample_id,
    sampleType = md$class,
    batchID = md$batch,
    injectionOrder = md$run_order
  ) |>
    dplyr::bind_cols(as.data.frame(t(x_safe), check.names = FALSE))

  train_qc_ids <- intersect(train_qc_ids, md$sample_id[md$class == "QC"])
  test_qc_ids <- intersect(test_qc_ids, md$sample_id[md$class == "QC"])
  if (length(train_qc_ids) < 2L) stop("TIGER holdout run requires at least 2 training QCs.")
  if (length(test_qc_ids) < 1L) stop("TIGER holdout run requires at least 1 test QC.")

  train_samples <- sample_df |>
    dplyr::filter(sampleID %in% train_qc_ids)

  if (isTRUE(duplicate_test_qc)) {
    nonqc_test <- if (isTRUE(include_nonqc_test)) {
      sample_df |>
        dplyr::filter(sampleType != "QC")
    } else {
      sample_df[0, , drop = FALSE]
    }

    qc_test_dup <- sample_df |>
      dplyr::filter(sampleID %in% test_qc_ids) |>
      dplyr::mutate(sampleID = paste0("QCtest__", sampleID), sampleType = "QC_test")

    test_samples <- dplyr::bind_rows(nonqc_test, qc_test_dup)
  } else {
    test_ids <- if (isTRUE(include_nonqc_test)) {
      unique(c(md$sample_id[md$class != "QC"], test_qc_ids))
    } else {
      unique(test_qc_ids)
    }
    test_samples <- sample_df |>
      dplyr::filter(sampleID %in% test_ids)
  }

  if (!isTRUE(use_injection)) {
    train_samples <- train_samples |>
      dplyr::select(-injectionOrder)
    test_samples <- test_samples |>
      dplyr::select(-injectionOrder)
  }

  corrected <- run_silently(
    TIGERr::run_TIGER(
      test_samples = test_samples,
      train_samples = train_samples,
      col_sampleID = "sampleID",
      col_sampleType = "sampleType",
      col_batchID = "batchID",
      col_order = if (isTRUE(use_injection)) "injectionOrder" else NULL,
      mtry_percent = mtry_percent,
      nodesize_percent = nodesize_percent,
      ntree = ntree,
      parallel.cores = parallel_cores
    )
  )

  meta_cols <- c("sampleID", "sampleType", "batchID")
  if (isTRUE(use_injection)) {
    meta_cols <- c(meta_cols, "injectionOrder")
  }

  feature_cols <- setdiff(colnames(corrected), meta_cols)
  m <- t(as.matrix(corrected[, feature_cols, drop = FALSE]))
  rownames(m) <- unname(feature_map[feature_cols])
  colnames(m) <- corrected$sampleID

  if (isTRUE(duplicate_test_qc)) {
    dup_cols <- grep("^QCtest__", colnames(m), value = TRUE)
    if (length(dup_cols) > 0L) {
      m_qc <- m[, dup_cols, drop = FALSE]
      colnames(m_qc) <- sub("^QCtest__", "", dup_cols)
      nondup_cols <- setdiff(colnames(m), dup_cols)
      m <- cbind(m[, nondup_cols, drop = FALSE], m_qc)
    }
  }

  output_ids <- md$sample_id[md$sample_id %in% colnames(m)]
  m[, output_ids, drop = FALSE]
}

run_tiger_all_corrected <- function(
  x_mat,
  control_ids,
  meta_df,
  use_injection = TRUE,
  mtry_percent = 0.4,
  nodesize_percent = 0.4,
  ntree = 5,
  parallel_cores = 1
) {
  run_tiger_with_holdout_qc(
    x_mat = x_mat,
    train_qc_ids = control_ids,
    test_qc_ids = control_ids,
    meta_df = meta_df,
    use_injection = use_injection,
    mtry_percent = mtry_percent,
    nodesize_percent = nodesize_percent,
    ntree = ntree,
    parallel_cores = parallel_cores,
    duplicate_test_qc = TRUE,
    include_nonqc_test = TRUE
  )
}

run_serrf_with_holdout_qc <- function(
  x_mat,
  train_qc_ids,
  test_qc_ids,
  meta_df,
  jitter_eps = 0,
  duplicate_test_qc = TRUE,
  include_nonqc_samples = TRUE
) {
  md <- align_meta_to_mat(x_mat, meta_df)

  if (!isTRUE(include_nonqc_samples)) {
    keep_ids <- md$sample_id[md$class == "QC"]
    md <- md[match(keep_ids, md$sample_id), , drop = FALSE]
    x_mat <- x_mat[, keep_ids, drop = FALSE]
  }

  train_qc_ids <- intersect(train_qc_ids, md$sample_id[md$class == "QC"])
  test_qc_ids <- intersect(test_qc_ids, md$sample_id[md$class == "QC"])
  if (length(train_qc_ids) < 2L) stop("SERRF holdout run requires at least 2 training QCs.")
  if (length(test_qc_ids) < 1L) stop("SERRF holdout run requires at least 1 test QC.")

  md_serrf <- md
  md_serrf$class <- as.character(md_serrf$class)
  md_serrf$class[md_serrf$class == "QC"] <- "QC_holdout"
  md_serrf$class[md_serrf$sample_id %in% train_qc_ids] <- "QC"

  x_serrf <- x_mat
  if (isTRUE(duplicate_test_qc)) {
    dup_ids <- paste0("QCtest__", test_qc_ids)
    dup_block <- x_mat[, test_qc_ids, drop = FALSE]
    colnames(dup_block) <- dup_ids
    x_serrf <- cbind(x_serrf, dup_block)

    dup_meta <- md[match(test_qc_ids, md$sample_id), , drop = FALSE]
    dup_meta$sample_id <- dup_ids
    dup_meta$class <- "QC_test"
    dup_meta$run_order <- dup_meta$run_order + seq_along(test_qc_ids) * 1e-4
    md_serrf <- dplyr::bind_rows(md_serrf, dup_meta)
  } else {
    md_serrf$class[md_serrf$sample_id %in% test_qc_ids] <- "QC_test"
  }

  md_serrf <- md_serrf[match(colnames(x_serrf), md_serrf$sample_id), , drop = FALSE]

  jitter_candidates <- unique(c(
    as.numeric(jitter_eps),
    if (isTRUE(jitter_eps <= 0)) c(1e-6, 1e-4) else c(as.numeric(jitter_eps) * 10)
  ))

  out <- NULL
  sample_map <- NULL
  feature_map <- NULL
  last_err <- NULL
  for (jit in jitter_candidates) {
    mal <- make_malbac_object(x_serrf, md_serrf, jitter_eps = jit)
    out_try <- tryCatch(
      run_silently(
        with_no_cluster(
          malbacR::bc_serrf(
            omicsData = mal$obj,
            sampletype_cname = "sample_type",
            test_val = "QC",
            group_cname = "sample_type"
          ),
          patch_do_parallel = TRUE
        )
      ),
      error = function(e) e
    )

    if (!inherits(out_try, "error")) {
      out <- out_try
      sample_map <- mal$sample_map
      feature_map <- mal$feature_map
      break
    }
    last_err <- out_try
  }

  if (is.null(out)) stop(last_err)

  m <- extract_malbac_matrix(out, sample_map, feature_map = feature_map)
  if (isTRUE(duplicate_test_qc)) {
    dup_cols <- paste0("QCtest__", test_qc_ids)
    dup_cols <- dup_cols[dup_cols %in% colnames(m)]
    if (length(dup_cols) > 0L) {
      m_qc <- m[, dup_cols, drop = FALSE]
      colnames(m_qc) <- sub("^QCtest__", "", dup_cols)
      nondup_cols <- setdiff(colnames(m), dup_cols)
      m <- cbind(m[, nondup_cols, drop = FALSE], m_qc)
    }
  }

  output_ids <- md$sample_id[md$sample_id %in% colnames(m)]
  m[, output_ids, drop = FALSE]
}

run_serrf_all_corrected <- function(
  x_mat,
  control_ids,
  meta_df,
  jitter_eps = 0
) {
  run_serrf_with_holdout_qc(
    x_mat = x_mat,
    train_qc_ids = control_ids,
    test_qc_ids = control_ids,
    meta_df = meta_df,
    jitter_eps = jitter_eps,
    duplicate_test_qc = TRUE,
    include_nonqc_samples = TRUE
  )
}

winn_ns_fun <- function(name) {
  get(name, envir = asNamespace("winn"))
}

run_winn_with_controls <- function(x_mat, control_ids = NULL, meta_df, auto_batch = FALSE, parameters = "auto") {
  md <- align_meta_to_mat(x_mat, meta_df)
  control_idx_local <- if (is.null(control_ids)) NULL else which(md$sample_id %in% control_ids)

  winn::winn(
    data = x_mat,
    batch = if (isTRUE(auto_batch)) NULL else md$batch,
    run_order = md$run_order,
    control_samples = control_idx_local,
    parameters = parameters
  )
}

run_winn_parameter_set <- function(
  x_mat,
  meta_df,
  control_ids,
  auto_batch = FALSE,
  test = "Ljung-Box",
  lag = NULL,
  acorr_fdr = 0.05,
  anova_fdr = 0.05,
  normalization = "shrink",
  spline_method = "conservative",
  scale_by_batch = FALSE,
  remove_batch_effects = "anova",
  pelt_penalty = NULL
) {
  md <- align_meta_to_mat(x_mat, meta_df)
  control_idx_local <- which(md$sample_id %in% control_ids)
  if (length(control_idx_local) < 2L) {
    stop("WiNN tuning/application requires at least 2 control samples.")
  }

  adjust_outliers <- winn_ns_fun("adjust_outliers_mad")
  auto_detect_batch <- winn_ns_fun(".auto_detect_batch")

  norm_data <- adjust_outliers(as.matrix(x_mat))
  batch_vec <- if (isTRUE(auto_batch)) {
    auto_detect_batch(norm_data, pelt_penalty = pelt_penalty)
  } else {
    md$batch
  }

  drift_corrected <- winn::autocorrelation_correct(
    norm_data,
    run_order = md$run_order,
    batch = batch_vec,
    lag = lag,
    test = test,
    detrend = "mean",
    fdr_threshold = acorr_fdr,
    spline_method = spline_method
  )

  batch_corrected <- if (identical(remove_batch_effects, "anova")) {
    winn::anova_batch_correction(drift_corrected, batch_vec, fdr_threshold = anova_fdr)
  } else {
    winn::combat_batch_correction(drift_corrected, batch_vec)
  }

  final_data <- if (identical(normalization, "none")) {
    batch_corrected
  } else {
    winn::normalize_by_dilution_factor(
      batch_corrected,
      processing = normalization,
      control_samples = control_idx_local
    )
  }

  if (isTRUE(scale_by_batch)) {
    final_data <- winn::scale_by_batch(final_data, batch_vec)
  }

  attr(final_data, "detected_batch") <- batch_vec
  attr(final_data, "pelt_penalty") <- attr(batch_vec, "pelt_penalty")
  final_data
}

run_winn_mode1_auto_qc <- function(x_mat, control_ids, meta_df, params = NULL) {
  if (is.null(params)) {
    return(run_winn_with_controls(x_mat, control_ids = control_ids, meta_df = meta_df, auto_batch = FALSE, parameters = "auto"))
  }

  run_winn_parameter_set(
    x_mat = x_mat,
    meta_df = meta_df,
    control_ids = control_ids,
    auto_batch = FALSE,
    test = as.character(params$test[1]),
    lag = if ("lag" %in% names(params) && is.finite(suppressWarnings(as.numeric(params$lag[1])))) as.integer(params$lag[1]) else NULL,
    acorr_fdr = as.numeric(params$acorr_fdr[1]),
    anova_fdr = as.numeric(params$anova_fdr[1]),
    normalization = as.character(params$normalization[1]),
    spline_method = as.character(params$spline_method[1]),
    scale_by_batch = as.logical(params$scale_by_batch[1]),
    pelt_penalty = NULL
  )
}

run_winn_mode2_auto_batch_qc <- function(x_mat, control_ids, meta_df, params = NULL) {
  if (is.null(params)) {
    return(run_winn_with_controls(x_mat, control_ids = control_ids, meta_df = meta_df, auto_batch = TRUE, parameters = "auto"))
  }

  pelt_penalty <- NULL
  if ("pelt_penalty" %in% names(params) && !is.na(params$pelt_penalty[1])) {
    pelt_penalty <- as.numeric(params$pelt_penalty[1])
  }

  run_winn_parameter_set(
    x_mat = x_mat,
    meta_df = meta_df,
    control_ids = control_ids,
    auto_batch = TRUE,
    test = as.character(params$test[1]),
    lag = if ("lag" %in% names(params) && is.finite(suppressWarnings(as.numeric(params$lag[1])))) as.integer(params$lag[1]) else NULL,
    acorr_fdr = as.numeric(params$acorr_fdr[1]),
    anova_fdr = as.numeric(params$anova_fdr[1]),
    normalization = as.character(params$normalization[1]),
    spline_method = as.character(params$spline_method[1]),
    scale_by_batch = as.logical(params$scale_by_batch[1]),
    pelt_penalty = pelt_penalty
  )
}

run_winn_mode3_default_no_qc <- function(x_mat, meta_df) {
  run_winn_with_controls(x_mat, control_ids = NULL, meta_df = meta_df, auto_batch = FALSE, parameters = "fixed")
}

weighted_pc_r2 <- function(mat, target, target_type = NULL, n_pcs = 5L) {
  if (!exists("weighted_pc_r2_explicit", mode = "function", inherits = TRUE)) {
    stop("Source scripts/weighted_pc_r2.R before using weighted_pc_r2().")
  }
  weighted_pc_r2_explicit(
    mat = mat,
    target = target,
    target_type = target_type,
    n_pcs = n_pcs
  )
}

cv_percent <- function(mat, sample_ids) {
  sample_ids <- sample_ids[sample_ids %in% colnames(mat)]
  if (length(sample_ids) < 2L) {
    return(stats::setNames(rep(NA_real_, nrow(mat)), rownames(mat)))
  }

  apply(mat[, sample_ids, drop = FALSE], 1, function(v) {
    ok <- is.finite(v)
    if (sum(ok) < 2L) return(NA_real_)
    mu <- mean(v[ok])
    if (!is.finite(mu) || mu <= 0) return(NA_real_)
    100 * stats::sd(v[ok]) / mu
  })
}

run_with_timeout <- function(code, timeout_sec = 600) {
  setTimeLimit(elapsed = timeout_sec, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
  force(code)
}

tune_rank <- function(x, higher_better = TRUE) {
  x <- as.numeric(x)
  if (higher_better) {
    x[!is.finite(x)] <- -Inf
    rank(x, ties.method = "average")
  } else {
    x[!is.finite(x)] <- Inf
    rank(-x, ties.method = "average")
  }
}

score_tuning_table <- function(df) {
  if (nrow(df) <= 1L) {
    df$score <- 1
    return(df)
  }

  df$score <-
    tune_rank(df$sample_identity_pc_r2, higher_better = TRUE) +
    tune_rank(df$plate_pc_r2, higher_better = FALSE) +
    tune_rank(df$control_cv_median, higher_better = FALSE) +
    tune_rank(df$replicate_icc_median, higher_better = TRUE) +
    tune_rank(df$feature_plate_icc_median, higher_better = TRUE)

  df
}

compute_pca_scores <- function(mat, meta_df, sample_ids = colnames(mat), n_pcs = 2L) {
  sample_ids <- sample_ids[sample_ids %in% colnames(mat)]
  if (length(sample_ids) < 3L) {
    return(tibble::tibble())
  }

  md <- meta_df[match(sample_ids, meta_df$sample_id), , drop = FALSE]
  x_log <- log1p(mat[, sample_ids, drop = FALSE])
  x_t <- t(x_log)

  keep <- apply(x_t, 2, function(v) {
    vv <- v[is.finite(v)]
    length(vv) >= 2L && stats::sd(vv) > 0
  })
  x_t <- x_t[, keep, drop = FALSE]
  if (nrow(x_t) < 3L || ncol(x_t) < 2L) {
    return(tibble::tibble())
  }

  pca <- stats::prcomp(x_t, center = TRUE, scale. = TRUE)
  pcs <- min(as.integer(n_pcs), ncol(pca$x))
  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

  tibble::as_tibble(pca$x[, seq_len(pcs), drop = FALSE]) |>
    dplyr::mutate(
      sample_id = sample_ids,
      sample = md$sample,
      role = md$role,
      role_plot = md$role_plot,
      batch = md$batch,
      plate = md$plate,
      run_order = md$run_order,
      pc1_var = var_exp[1],
      pc2_var = var_exp[min(2L, length(var_exp))]
    )
}

detect_winn_auto_batch <- function(data_mat, pelt_penalty = NULL) {
  detect_fun <- get(".auto_detect_batch", envir = asNamespace("winn"))
  detect_fun(data_mat, pelt_penalty = pelt_penalty)
}
