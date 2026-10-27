# Revision adapter: neutral IDs and training-QC-only abundance reference.
# This is a new analysis protocol; frozen published results are preserved.
run_qc_rfsc_review <- function(x_mat, control_ids, meta_df, ntree = 500, coCV = 30, Frule = 0.8, output_scale = c("training_qc_abundance", "relative_ratio")) {
  output_scale <- match.arg(output_scale)
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
  neutral_ids <- sprintf("Injection%07d", seq_len(nrow(md)))
  stat_names <- ifelse(qc_flag, paste0(neutral_ids, "_QC"), neutral_ids)
  sample_map <- stats::setNames(md$sample_id, stat_names)
  if (!identical(grepl("QC|qc|Qc", stat_names), qc_flag)) {
    stop("statTarget name-based training indices differ from the authorized controls.")
  }
  abundance_reference <- apply(x_ord[, qc_flag, drop = FALSE], 1L, stats::median)
  if (any(!is.finite(abundance_reference) | abundance_reference <= 0)) {
    stop("Training-QC abundance reference must be finite and positive for every feature.")
  }
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

  previous_options <- options(qcrfsc_authorized_names = stat_names[qc_flag])
  on.exit(options(previous_options), add = TRUE)
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

  corrected <- matrix(as.numeric(as.matrix(corrected_df)), nrow = nrow(corrected_df), ncol = ncol(corrected_df), dimnames = dimnames(corrected_df))
  rownames(corrected) <- unname(feature_map[rownames(corrected_df)])
  colnames(corrected) <- unname(sample_map[colnames(corrected_df)])
  if (anyNA(rownames(corrected)) || anyNA(colnames(corrected))) {
    stop("statTarget returned unknown feature or sample identifiers.")
  }
  corrected <- corrected / 1000
  if (identical(output_scale, "training_qc_abundance")) {
    corrected <- corrected * abundance_reference[rownames(corrected)]
  }
  missing_features <- setdiff(rownames(x_mat), rownames(corrected))
  missing_samples <- setdiff(colnames(x_mat), colnames(corrected))
  if (identical(output_scale, "relative_ratio") && (length(missing_features) || length(missing_samples))) {
    stop("Relative-ratio output cannot contain raw-abundance fallback entries.")
  }
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
  attr(full, "qcrfsc_output_scale") <- output_scale
  attr(full, "qcrfsc_sample_id_map") <- data.frame(original_id = md$sample_id, internal_id = stat_names, training_reference = qc_flag)
  attr(full, "qcrfsc_training_ids") <- md$sample_id[qc_flag]
  attr(full, "qcrfsc_abundance_reference") <- abundance_reference
  full
}
