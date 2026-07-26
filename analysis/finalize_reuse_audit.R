#!/usr/bin/env Rscript

release_root <- Sys.getenv("WINN_RELEASE_ROOT", unset = "")
if (!nzchar(release_root) || !dir.exists(release_root)) {
  stop("WINN_RELEASE_ROOT must identify the release directory.", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE) ||
    !requireNamespace("jsonlite", quietly = TRUE)) {
  stop("digest and jsonlite are required.", call. = FALSE)
}

audit_path <- file.path(
  release_root, "analysis", "config", "competitor_reuse_audit.tsv"
)
audit <- read.delim(audit_path, stringsAsFactors = FALSE, check.names = FALSE)
tasks <- read.csv(
  file.path(release_root, "analysis", "config", "primary_run_matrix.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
datasets <- read.csv(
  file.path(release_root, "analysis", "config", "dataset_manifest.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
expected_package_sha <-
  "71a0964cee2778b2e5789d20621147e074c7945e813cf76af2ceeb696104aae1"

checks <- lapply(seq_len(nrow(audit)), function(index) {
  row <- audit[index, , drop = FALSE]
  task <- tasks[
    tasks$dataset_key == row$dataset_key & tasks$method_id == row$method_id,
    , drop = FALSE
  ]
  if (nrow(task) != 1L) {
    stop("No unique primary task for ", row$dataset_key, " / ", row$method_id)
  }
  directory <- file.path(release_root, as.character(task$output_dir))
  required <- file.path(
    directory, c("complete.json", "run_manifest.json", "corrected_matrix.rds")
  )
  if (!all(file.exists(required))) {
    stop("Incomplete primary artifact: ", row$dataset_key, " / ", row$method_id)
  }
  manifest <- jsonlite::read_json(
    file.path(directory, "run_manifest.json"), simplifyVector = TRUE
  )
  matrix_path <- file.path(directory, "corrected_matrix.rds")
  value <- readRDS(matrix_path)
  recorded_decision <- as.character(manifest$reuse_decision)
  decision_ok <-
    identical(recorded_decision, as.character(row$decision)) ||
    identical(recorded_decision, as.character(row$public_reproduction_decision))
  package <- manifest$package
  expected_input_sha <- datasets$matrix_sha256[
    match(as.character(row$dataset_key), datasets$dataset_key)
  ]
  passed <-
    identical(as.character(manifest$status), "completed") &&
    decision_ok &&
    length(expected_input_sha) == 1L && !is.na(expected_input_sha) &&
    identical(as.character(manifest$input_object_sha256), expected_input_sha) &&
    identical(as.integer(dim(value)), c(as.integer(manifest$n_features), as.integer(manifest$n_samples))) &&
    all(is.finite(value)) && all(value >= 0) &&
    identical(
      digest::digest(value, algo = "sha256", serialize = TRUE),
      as.character(manifest$output_object_sha256)
    ) &&
    identical(
      digest::digest(file = matrix_path, algo = "sha256"),
      as.character(manifest$matrix_file_sha256)
    ) &&
    identical(as.character(package$version), "0.1.4") &&
    identical(as.character(package$source_archive_sha256), expected_package_sha)
  data.frame(
    dataset_key = as.character(row$dataset_key),
    method_id = as.character(row$method_id), passed = passed,
    stringsAsFactors = FALSE
  )
})
checks <- do.call(rbind, checks)
if (nrow(checks) != 36L || !all(checks$passed)) {
  failed <- checks[!checks$passed, , drop = FALSE]
  stop(
    "Competitor artifact validation failed for: ",
    paste(paste(failed$dataset_key, failed$method_id, sep = "/"), collapse = ", "),
    call. = FALSE
  )
}

audit$final_validation <- "passed"
temporary <- tempfile("competitor_reuse_audit_", tmpdir = dirname(audit_path))
on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
write.table(
  audit, temporary, sep = "\t", row.names = FALSE, col.names = TRUE,
  quote = FALSE, na = ""
)
if (!file.rename(temporary, audit_path)) {
  stop("Could not atomically update the reuse audit.", call. = FALSE)
}
message("Validated 36 primary Raw/competitor artifacts and finalized the reuse audit.")
