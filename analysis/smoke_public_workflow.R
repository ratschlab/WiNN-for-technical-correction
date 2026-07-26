#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
.libPaths(c(file.path(repo_root, "Rlib"), .libPaths()))

required <- c("digest", "jsonlite", "winn")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing package(s): ", paste(missing, collapse = ", "), call. = FALSE)

tarball <- file.path(repo_root, "package", "winn_0.1.4.tar.gz")
expected_sha <- "71a0964cee2778b2e5789d20621147e074c7945e813cf76af2ceeb696104aae1"
observed_sha <- digest::digest(file = tarball, algo = "sha256")
if (!identical(observed_sha, expected_sha)) stop("Frozen WiNN archive checksum mismatch.", call. = FALSE)
source(file.path(repo_root, "package", "assert_frozen_winn.R"), local = TRUE)
frozen_winn <- assert_frozen_winn(release_root = repo_root)
source(file.path(repo_root, "scripts", "weighted_pc_r2.R"), local = TRUE)

bundle <- file.path(repo_root, "data", "simulated", "canonical", "SIM01")
x <- readRDS(file.path(bundle, "raw_intensity.rds"))
truth <- readRDS(file.path(bundle, "clean_ground_truth.rds"))
meta <- read.csv(file.path(bundle, "sample_metadata.csv"), stringsAsFactors = FALSE)
meta <- meta[match(colnames(x), meta$sample_id), , drop = FALSE]
if (anyNA(meta$sample_id) || !identical(colnames(x), meta$sample_id)) {
  stop("Canonical simulation metadata do not align with the matrix.", call. = FALSE)
}

feature_index <- seq_len(min(30L, nrow(x)))
sample_index <- unlist(lapply(split(seq_len(nrow(meta)), meta$batch), function(index) {
  index[seq_len(min(24L, length(index)))]
}), use.names = FALSE)
x_small <- x[feature_index, sample_index, drop = FALSE]
meta_small <- meta[sample_index, , drop = FALSE]

fit_once <- function(matrix, metadata) {
  winn::winn(
    data = matrix, batch = metadata$batch, run_order = metadata$run_order,
    control_samples = NULL, parameters = "fixed", ljung_box_fitdf = 0L,
    return_details = TRUE
  )
}

baseline <- fit_once(x_small, meta_small)
set.seed(42)
permutation <- sample(seq_len(ncol(x_small)))
permuted <- fit_once(x_small[, permutation, drop = FALSE], meta_small[permutation, , drop = FALSE])
restored <- permuted$data[, order(permutation), drop = FALSE]
permutation_difference <- max(abs(baseline$data - restored))

simulation_hashes <- read.csv(
  file.path(repo_root, "analysis", "config", "simulation_bundle_hashes.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
sim01_hashes <- simulation_hashes[simulation_hashes$seed_id == "SIM01", , drop = FALSE]
simulation_identity_ok <- nrow(simulation_hashes) == 20L && nrow(sim01_hashes) == 2L &&
  digest::digest(x, algo = "sha256", serialize = TRUE) ==
    sim01_hashes$reference_object_sha256[sim01_hashes$object == "raw_intensity"] &&
  digest::digest(truth, algo = "sha256", serialize = TRUE) ==
    sim01_hashes$reference_object_sha256[sim01_hashes$object == "clean_ground_truth"]

selection <- read.csv(
  file.path(repo_root, "analysis", "config", "endpoint_free_selection_manifest.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
unavailable <- c(
  "hidden_qc_values_available", "biological_labels_available",
  "replicate_identities_available", "preservation_endpoints_available",
  "batch_r2_available", "final_rank_available"
)
selection_separation_ok <- all(unavailable %in% names(selection)) &&
  all(vapply(unavailable, function(field) !any(as.logical(selection[[field]])), logical(1)))
selection_provenance_fields <- c(
  "selection_seed", "source_package", "source_version",
  "information_available_during_selection", "departure_justification",
  "configuration_frozen_before_evaluation"
)
selection_provenance_ok <-
  nrow(selection) == 54L && all(selection_provenance_fields %in% names(selection)) &&
  all(is.finite(as.numeric(selection$selection_seed))) &&
  all(nzchar(selection$source_package)) && all(nzchar(selection$source_version)) &&
  all(selection$configuration_frozen_before_evaluation)

batch_guard_error <- try(
  weighted_pc_r2_explicit(x_small, as.numeric(factor(meta_small$batch))),
  silent = TRUE
)
batch_categorical_value <- weighted_pc_r2_explicit(
  x_small, meta_small$batch, target_type = "categorical"
)
categorical_batch_guard_ok <-
  inherits(batch_guard_error, "try-error") && is.finite(batch_categorical_value)

release_expressions <- parse(file.path(repo_root, "analysis", "release_helpers.R"))
validator_expression <- vapply(release_expressions, function(value) {
  is.call(value) && identical(value[[1L]], as.name("<-")) &&
    identical(value[[2L]], as.name("release_validate_output"))
}, logical(1))
if (sum(validator_expression) != 1L) {
  stop("Could not isolate the output-domain validator for its smoke test.", call. = FALSE)
}
domain_environment <- new.env(parent = baseenv())
eval(release_expressions[[which(validator_expression)]], envir = domain_environment)
domain_template <- matrix(
  c(1, 2), nrow = 1,
  dimnames = list("feature_1", c("sample_1", "sample_2"))
)
domain_candidate <- matrix(
  c(-2, 3), nrow = 1, dimnames = dimnames(domain_template)
)
domain_result <- domain_environment$release_validate_output(
  domain_candidate, list(x = domain_template), "smoke"
)
domain_record <- attr(domain_result, "intensity_floor", exact = TRUE)
nonfinite_error <- try(
  domain_environment$release_validate_output(
    matrix(c(NA_real_, 1), nrow = 1, dimnames = dimnames(domain_template)),
    list(x = domain_template), "smoke"
  ),
  silent = TRUE
)
output_domain_guard_ok <-
  min(domain_result) == 0 &&
  !is.null(domain_record) && domain_record$n_values == 1L &&
  domain_record$minimum_before_floor == -2 &&
  inherits(nonfinite_error, "try-error")

benchmark_expressions <- parse(file.path(repo_root, "scripts", "benchmark_helpers.R"))
tiger_fallback_expression <- vapply(benchmark_expressions, function(value) {
  is.call(value) && identical(value[[1L]], as.name("<-")) &&
    identical(value[[2L]], as.name("restore_tiger_nonfinite_features"))
}, logical(1))
if (sum(tiger_fallback_expression) != 1L) {
  stop("Could not isolate the TIGER non-finite fallback for its smoke test.", call. = FALSE)
}
tiger_environment <- new.env(parent = baseenv())
eval(benchmark_expressions[[which(tiger_fallback_expression)]], envir = tiger_environment)
tiger_input <- matrix(
  1:6, nrow = 2,
  dimnames = list(c("feature_1", "feature_2"), paste0("sample_", 1:3))
)
tiger_candidate <- tiger_input + 10
tiger_candidate[1, 2] <- NA_real_
tiger_result <- tiger_environment$restore_tiger_nonfinite_features(
  tiger_candidate, tiger_input
)
tiger_record <- attr(tiger_result, "tiger_nonfinite_fallback", exact = TRUE)
tiger_finite_result <- tiger_environment$restore_tiger_nonfinite_features(
  tiger_input + 10, tiger_input
)
tiger_nonfinite_fallback_ok <-
  identical(as.numeric(tiger_result[1, ]), as.numeric(tiger_input[1, ])) &&
  identical(as.numeric(tiger_result[2, ]), as.numeric(tiger_candidate[2, ])) &&
  isTRUE(tiger_record$applied) && tiger_record$n_nonfinite_values == 1L &&
  tiger_record$n_features == 1L && tiger_record$n_samples == 1L &&
  identical(as.vector(tiger_finite_result), as.vector(tiger_input + 10)) &&
  identical(dimnames(tiger_finite_result), dimnames(tiger_input))

checks <- data.frame(
  check = c(
    "frozen_tarball_sha256", "installed_version", "isolated_installation_path",
    "matrix_metadata_alignment",
    "finite_output", "dimensions_preserved", "permutation_invariance",
    "canonical_simulation_identity", "selection_endpoint_separation",
    "selection_software_seed_provenance", "categorical_batch_guard",
    "nonnegative_output_domain_guard", "tiger_nonfinite_feature_fallback"
  ),
  passed = c(
    observed_sha == expected_sha, as.character(utils::packageVersion("winn")) == "0.1.4",
    startsWith(
      frozen_winn$package_path,
      paste0(normalizePath(file.path(repo_root, "Rlib")), .Platform$file.sep)
    ),
    identical(colnames(x_small), meta_small$sample_id), all(is.finite(baseline$data)),
    identical(dim(baseline$data), dim(x_small)), permutation_difference < 1e-10,
    simulation_identity_ok, selection_separation_ok, selection_provenance_ok,
    categorical_batch_guard_ok, output_domain_guard_ok,
    tiger_nonfinite_fallback_ok
  ),
  detail = c(
    observed_sha, as.character(utils::packageVersion("winn")),
    frozen_winn$package_path,
    paste(nrow(x_small), "features x", ncol(x_small), "samples"),
    as.character(sum(!is.finite(baseline$data))), paste(dim(baseline$data), collapse = "x"),
    format(permutation_difference, scientific = TRUE),
    paste(nrow(simulation_hashes), "frozen matrix-identity records"),
    paste(length(unavailable), "evaluation fields unavailable during selection"),
    paste(nrow(selection), "selection records with software and seed provenance"),
    "batch metric requires explicit categorical target type",
    "negative values floored with provenance; non-finite values rejected",
    "one affected TIGER feature restored in full; finite features unchanged"
  ),
  stringsAsFactors = FALSE
)
output_dir <- file.path(repo_root, "results", "smoke")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(checks, file.path(output_dir, "smoke_checks.csv"), row.names = FALSE)
utils::capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
if (!all(checks$passed)) stop("Public smoke workflow failed.", call. = FALSE)
message("Public smoke workflow passed (", nrow(checks), " checks).")
