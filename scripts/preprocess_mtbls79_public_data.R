# Manuscript-aligned preprocessing for the public MTBLS79 benchmark.
# This script implements the missing-value handling described in
# Supplementary Methods S7.

preprocess_mtbls79_public_data <- function(
  file_path = file.path("data", "public", "raw", "Dataset07__SFPM.xlsx"),
  output_dir = file.path("data", "public", "processed"),
  seed = 42L,
  reuse_existing = FALSE,
  verbose = interactive()
) {
  required_pkgs <- c("missForest", "openxlsx")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0L) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }

  if (!file.exists(file_path)) {
    stop(
      "Missing input file: ", file_path, "\n",
      "Download Dataset07__SFPM.xlsx from https://www.ebi.ac.uk/metabolights/MTBLS79 ",
      "and place it at the path above before running this script."
    )
  }

  workbook <- openxlsx::loadWorkbook(file_path)
  raw_data <- openxlsx::readWorkbook(workbook, sheet = "data")
  raw_meta <- openxlsx::readWorkbook(workbook, sheet = "meta")

  data_clean <- raw_data[-nrow(raw_data), , drop = FALSE]
  meta_clean <- raw_meta[-nrow(raw_meta), , drop = FALSE]

  sample_names <- data_clean[[1]]
  rownames(data_clean) <- sample_names
  data_clean <- data_clean[, -1, drop = FALSE]

  rownames(meta_clean) <- sample_names
  meta_clean <- meta_clean[, 1:4, drop = FALSE]

  data_clean[data_clean == 0] <- NA

  presence_threshold <- 0.8
  keep_features <- colSums(!is.na(data_clean)) >= (nrow(data_clean) * presence_threshold)
  filtered_data <- data_clean[, keep_features, drop = FALSE]

  colnames(filtered_data) <- paste0("metabolite", seq_len(ncol(filtered_data)))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  imputed_path <- file.path(output_dir, "MTBLS79_imputed_data.csv")
  metadata_path <- file.path(output_dir, "MTBLS79_metadata.csv")

  if (isTRUE(reuse_existing) && file.exists(imputed_path) && file.exists(metadata_path)) {
    imputed_data <- utils::read.csv(
      imputed_path,
      row.names = 1,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    meta_clean <- utils::read.csv(
      metadata_path,
      row.names = 1,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    set.seed(seed)
    imputation <- missForest::missForest(
      filtered_data,
      verbose = isTRUE(verbose)
    )
    imputed_data <- imputation$ximp

    utils::write.csv(imputed_data, imputed_path, row.names = TRUE)
    utils::write.csv(meta_clean, metadata_path, row.names = TRUE)
  }

  feature_summary <- list(
    n_samples = nrow(data_clean),
    n_raw_features = ncol(data_clean),
    n_retained_features = ncol(filtered_data),
    n_removed_features = ncol(data_clean) - ncol(filtered_data),
    presence_threshold = presence_threshold
  )

  if (isTRUE(verbose)) {
    message("Saved processed MTBLS79 files to: ", normalizePath(output_dir, winslash = "/"))
    message(
      "Features retained after the 80% rule: ",
      feature_summary$n_retained_features,
      " / ",
      feature_summary$n_raw_features
    )
  }

  invisible(list(
    raw_data = data_clean,
    raw_metadata = meta_clean,
    filtered_data = filtered_data,
    imputed_data = imputed_data,
    output_paths = c(
      imputed_data = imputed_path,
      metadata = metadata_path
    ),
    feature_summary = feature_summary
  ))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  file_path <- if (length(args) >= 1L) args[[1]] else file.path("data", "public", "raw", "Dataset07__SFPM.xlsx")
  output_dir <- if (length(args) >= 2L) args[[2]] else file.path("data", "public", "processed")
  preprocess_mtbls79_public_data(file_path = file_path, output_dir = output_dir, verbose = TRUE)
}
