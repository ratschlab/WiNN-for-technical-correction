#!/usr/bin/env Rscript
# Aggregate-only audit of native statTarget zero-to-missing preprocessing.
# This reads the frozen prepared inputs; it does not fit or modify any matrix.
args <- commandArgs(TRUE)
arg <- function(key, default) {
  z <- args[startsWith(args, paste0(key, "="))]
  if (length(z)) sub(paste0("^", key, "="), "", z[1]) else default
}
input_root <- arg("--input-root", ".")
output <- arg("--output", "qcrfsc_preprocessing_input_audit.csv")
stopifnot(as.character(utils::packageVersion("statTarget")) == "1.38.0")
filter_mv <- get("FilterMV", envir = asNamespace("statTarget"))
datasets <- c("mtbls79", "clinical_fiams", "batchcorr_set1", "sacurine", "waveica", sprintf("SIM%02d", 1:10))
rows <- lapply(datasets, function(dataset) {
  path <- file.path(input_root, paste0("input_", dataset, ".rds"))
  p <- readRDS(path)
  x <- p$x
  training <- intersect(p$training_ids, colnames(x))
  heldout <- intersect(unique(c(p$hidden_ids, p$external_ids)), colnames(x))
  groups <- list(all = colnames(x), training_qc = training,
                 heldout_reference = heldout,
                 other = setdiff(colnames(x), c(training, heldout)))
  native_input <- t(x)
  native_input[native_input <= 0] <- NA_real_
  native_class <- addNA(as.factor(ifelse(colnames(x) %in% training, NA_integer_, 1L)))
  filtered <- as.matrix(filter_mv(native_input, native_class, rule = 0.8))
  stopifnot(nrow(filtered) == ncol(x))
  do.call(rbind, lapply(names(groups), function(group) {
    z <- x[, groups[[group]], drop = FALSE]
    data.frame(dataset = dataset, injection_subset = group,
               n_features = nrow(z), n_injections = ncol(z), n_cells = length(z),
               n_zero = sum(z == 0, na.rm = TRUE), n_negative = sum(z < 0, na.rm = TRUE),
               n_na = sum(is.na(z)), n_nonfinite = sum(!is.finite(z)),
               n_features_with_zero = sum(rowSums(z == 0, na.rm = TRUE) > 0),
               n_native_missing_after_nonpositive_conversion = sum(is.na(z) | z <= 0, na.rm = TRUE),
               n_features_after_FilterMV = ncol(filtered),
               n_missing_after_FilterMV = sum(is.na(filtered[match(groups[[group]], colnames(x)), , drop = FALSE])),
               matrix_sha256 = digest::digest(x, algo = "sha256"),
               input_file_sha256 = digest::digest(file = path, algo = "sha256"))
  }))
})
out <- do.call(rbind, rows)
write.csv(out, output, row.names = FALSE)
print(out[out$injection_subset == "all", c("dataset", "n_cells", "n_zero", "n_negative", "n_na", "n_nonfinite")], row.names = FALSE)
