# Chronological adapter: neutral IDs and training-QC-only abundance reference.
# Isolated copies of statTarget functions retain the native forest algorithm.
# Only REGfit's first/last-column QC guards are removed: RF predictions outside
# the training-QC range use the fitted terminal leaves without relocating QCs.
.qcrfsc_chronological_environment <- new.env(parent = asNamespace("statTarget"))
.qcrfsc_native_regfit <- getFromNamespace("REGfit", "statTarget")
.qcrfsc_statements <- as.list(body(.qcrfsc_native_regfit))[-1L]
.qcrfsc_remove <- vapply(.qcrfsc_statements, function(e) {
  s <- paste(deparse(e), collapse=" ")
  grepl("^st_QC <-|^ed_QC <-|^if \\(length\\(st_QC\\)|^if \\(length\\(ed_QC\\)", s)
}, logical(1))
stopifnot(sum(.qcrfsc_remove) == 4L)
.qcrfsc_assert <- quote({
  cn_audit <- colnames(x)
  expected_audit <- getOption("qcrfsc_authorized_names")
  found_audit <- cn_audit[grepl("QC|qc|Qc", cn_audit)]
  stopifnot(identical(cn_audit, getOption("qcrfsc_chronological_names")),
            identical(as.numeric(y), as.numeric(seq_along(cn_audit))),
            setequal(found_audit, expected_audit), length(found_audit)>=2L)
  if (exists(".rfsc_inner_audit", envir=.GlobalEnv, inherits=FALSE)) {
    audit_rows <- get(".rfsc_inner_audit", envir=.GlobalEnv)
    audit_rows[[length(audit_rows)+1L]] <- data.frame(call=length(audit_rows)+1L,
      available=length(cn_audit),authorized=length(expected_audit),used=length(found_audit),
      unauthorized=sum(!found_audit %in% expected_audit),chronological=TRUE)
    assign(".rfsc_inner_audit",audit_rows,envir=.GlobalEnv)
  }
})
body(.qcrfsc_native_regfit) <- as.call(c(list(as.name("{"),.qcrfsc_assert),
                                      .qcrfsc_statements[!.qcrfsc_remove]))
environment(.qcrfsc_native_regfit) <- .qcrfsc_chronological_environment
.qcrfsc_chronological_environment$REGfit <- .qcrfsc_native_regfit
.qcrfsc_shiftCor <- getFromNamespace("shiftCor", "statTarget")
environment(.qcrfsc_shiftCor) <- .qcrfsc_chronological_environment

run_qc_rfsc_review <- function(x_mat, control_ids, meta_df, ntree = 500, coCV = 30, Frule = 0.8, output_scale = c("training_qc_abundance", "relative_ratio")) {
  output_scale <- match.arg(output_scale)
  md <- align_meta_to_mat(x_mat, meta_df) |>
    dplyr::arrange(run_order)
  x_ord <- x_mat[, md$sample_id, drop = FALSE]

  control_ids <- intersect(control_ids, md$sample_id)
  if (length(control_ids) < 2L) {
    stop("QC-RFSC requires at least 2 QC controls.")
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

  samFile <- data.frame(name = safe_feature_ids, x_ord, check.names = FALSE)
  colnames(samFile)[-1] <- stat_names

  tmp <- tempfile("shiftcor_")
  dir.create(tmp)
  old <- getwd()
  setwd(tmp)
  on.exit(setwd(old), add = TRUE)

  write.csv(samPeno, "samPeno.csv", row.names = FALSE, quote = FALSE)
  write.csv(samFile, "samFile.csv", row.names = FALSE, quote = FALSE)

  previous_options <- options(qcrfsc_authorized_names = stat_names[qc_flag],
                              qcrfsc_chronological_names = stat_names)
  on.exit(options(previous_options), add = TRUE)
  run_silently(
    .qcrfsc_shiftCor(
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
  attr(full, "qcrfsc_order_policy") <- "chronological injection ranks; no QC relocation; native RF endpoint extrapolation"
  attr(full, "qcrfsc_outside_training_range_count") <- sum(seq_len(nrow(md)) < min(which(qc_flag)) | seq_len(nrow(md)) > max(which(qc_flag)))
  full
}
