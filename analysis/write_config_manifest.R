#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
config_root <- file.path(repo_root, "analysis", "config")
manifest_path <- file.path(config_root, "file_manifest.csv")
holdout_root <- file.path(repo_root, "config")

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("The digest package is required.", call. = FALSE)
}

config_files <- c(
  list.files(config_root, full.names = TRUE, recursive = FALSE),
  list.files(holdout_root, full.names = TRUE, recursive = TRUE)
)
config_files <- config_files[
  !file.info(config_files)$isdir &
    normalizePath(config_files, mustWork = TRUE) != normalizePath(manifest_path, mustWork = TRUE)
]
config_files <- sort(config_files)
observed <- data.frame(
  file = substring(config_files, nchar(repo_root) + 2L),
  bytes = as.numeric(file.info(config_files)$size),
  sha256 = unname(vapply(
    config_files, digest::digest, character(1), file = TRUE, algo = "sha256"
  )),
  stringsAsFactors = FALSE
)
rownames(observed) <- NULL

check_only <- "--check" %in% commandArgs(trailingOnly = TRUE)
if (check_only) {
  if (!file.exists(manifest_path)) stop("The configuration manifest is missing.", call. = FALSE)
  expected <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  expected$bytes <- as.numeric(expected$bytes)
  rownames(expected) <- NULL
  if (!identical(expected, observed)) {
    stop("One or more frozen configuration files do not match the manifest.", call. = FALSE)
  }
  message("Configuration manifest passed for ", nrow(observed), " files.")
} else {
  write.csv(observed, manifest_path, row.names = FALSE)
  message("Wrote configuration manifest for ", nrow(observed), " files.")
}
