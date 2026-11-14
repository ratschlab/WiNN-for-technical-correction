#!/usr/bin/env Rscript
# Reproduce integer-vs-double CSV metadata, then reject each integrity failure.
root <- tempfile("hash_audit_")
base <- file.path(root, "results/frozen_endpoint_free_inputs")
datasets <- c("simulation", "mtbls79", "clinical_fiams", "batchcorr_set1", "sacurine", "waveica")
methods <- c("raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf", "winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc")
rows <- list()
for (dataset in datasets) for (method in methods) {
  rel <- file.path("results/primary", dataset, method, "corrected_matrix.rds")
  path <- file.path(base, sub("^results/", "", rel))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeBin(as.raw(1:32), path)
  rows[[length(rows) + 1L]] <- data.frame(category = "primary", path = rel, bytes = 32L,
                                         sha256 = digest::digest(file = path, algo = "sha256"))
}
pre <- do.call(rbind, rows)
dir.create(file.path(base, "provenance"))
manifest <- file.path(base, "provenance/matrix_checksums_pre.csv")
write.csv(pre, manifest, row.names = FALSE)
stopifnot(typeof(read.csv(manifest)$bytes) == "integer")
run <- function() {
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("analysis/audit_derivative_input_hashes.R", shQuote(root)), stdout = FALSE, stderr = FALSE)
  as.integer(status)
}
stopifnot(run() == 0L)
first <- file.path(base, sub("^results/", "", pre$path[1]))
writeBin(as.raw(32:1), first) # Same length, different bytes must fail the hash check.
stopifnot(run() != 0L)
writeBin(as.raw(1:32), first)
wrong_size <- pre; wrong_size$bytes[1] <- 33L
write.csv(wrong_size, manifest, row.names = FALSE)
stopifnot(run() != 0L)
write.csv(pre[-1, ], manifest, row.names = FALSE)
stopifnot(run() != 0L)
unlink(root, recursive = TRUE)
cat("PASS: integer byte metadata accepted; changed content, wrong exact size and incomplete set rejected.\n")
