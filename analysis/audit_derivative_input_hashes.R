#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1])
base <- file.path(root, "results/frozen_endpoint_free_inputs")
pre <- read.csv(file.path(base, "provenance/matrix_checksums_pre.csv"), stringsAsFactors = FALSE)
datasets <- c("simulation", "mtbls79", "clinical_fiams", "batchcorr_set1", "sacurine", "waveica")
methods <- c("raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf",
             "winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc")
expected <- file.path("results/primary", rep(datasets, each = length(methods)),
                      rep(methods, times = length(datasets)), "corrected_matrix.rds")
stopifnot(nrow(pre) == 54L, !anyDuplicated(pre$path), setequal(pre$path, expected),
          !anyNA(pre$sha256), all(grepl("^[0-9a-f]{64}$", pre$sha256)),
          all(is.finite(pre$bytes)), all(pre$bytes >= 0))
post <- pre
for (i in seq_len(nrow(pre))) {
  f <- file.path(base, sub("^results/", "", pre$path[i]))
  post$bytes[i] <- file.info(f)$size
  post$sha256[i] <- digest::digest(file = f, algo = "sha256")
}
# CSV integers and file.info() doubles encode the same byte counts. Compare
# exact values rather than storage types; never use a relative size tolerance.
stopifnot(identical(pre$path, post$path), identical(pre$sha256, post$sha256),
          all(is.finite(post$bytes)), all(pre$bytes == post$bytes))
write.csv(post, file.path(base, "provenance/matrix_checksums_post.csv"), row.names = FALSE)
cat("PASS: all", nrow(post), "derivative input matrices remain byte-identical.\n")
