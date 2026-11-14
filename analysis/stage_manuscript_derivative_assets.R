#!/usr/bin/env Rscript

# Stage manuscript tables and figures generated from frozen derivative results.
# This script copies existing artifacts only; it never executes an analysis.

args <- commandArgs(trailingOnly = TRUE)
argument <- grep("^--manuscript-root=", args, value = TRUE)
if (length(argument) != 1L) {
  stop("Usage: Rscript analysis/stage_manuscript_derivative_assets.R --manuscript-root=/path/to/manuscript", call. = FALSE)
}
manuscript_root <- normalizePath(sub("^--manuscript-root=", "", argument), mustWork = TRUE)
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
if (!requireNamespace("digest", quietly = TRUE)) stop("Missing package: digest", call. = FALSE)

mapping <- data.frame(
  source = c(
    file.path("results/final/derivative/manuscript", c(
      "main_simulation_rows.tex", "main_mtbls79_rows.tex", "main_added_public_rows.tex",
      "main_clinical_rows.tex", "main_noqc_rows.tex", "full_benchmark_tables.tex",
      "supp_derivative_biology_rows.tex", "ablation_tables.tex", "supp_categorical_r2_trajectory.tex",
      "supp_reference_split_stability.tex", "supp_partial_confounding.tex"
    )),
    file.path("results/final/derivative/manuscript", c(
      "main_simulation_rows.csv", "main_mtbls79_rows.csv", "main_added_public_rows.csv",
      "main_clinical_rows.csv", "main_noqc_rows.csv", "full_benchmark_tables.csv"
    )),
    file.path("results/final/derivative/manuscript", c(
      "manuscript_key_biological_results.csv", "manuscript_association_calibration.csv",
      "manuscript_lobo_summary.csv", "supp_categorical_r2_trajectory.csv",
      "supp_reference_split_stability.csv", "supp_partial_confounding.csv"
    )),
    "results/final/tables/ablation_stage_impacts.csv",
    file.path("results/final/derivative/figures", c(
      "mtbls79_fixed_benchmark.pdf", "clinical_pair_aware_derivative.pdf",
      "waveica_technical_biological_tradeoff.pdf", "association_permutation_calibration.pdf",
      "leave_one_batch_out_effect_replication.pdf", "replicate_preservation.pdf",
      "categorical_batch_r2_trajectory.pdf"
    )),
    "results/final/figures/ablation_stage_impacts.pdf",
    "results/final/figures/selective_vs_forced.pdf",
    "results/final/figures/simulation_seed_stability.pdf",
    "results/final/figures/reference_split_stability.pdf",
    "results/final/figures/partial_confounding.pdf",
    file.path("results/final/derivative/figures/source_data", c(
      "mtbls79_fixed_benchmark.csv", "clinical_pair_aware_derivative.csv",
      "waveica_technical_biological_tradeoff.csv", "association_permutation_calibration.csv",
      "leave_one_batch_out_effect_replication.csv", "replicate_preservation.csv",
      "categorical_batch_r2_trajectory.csv"
    )),
    "results/final/figures/source_data/ablation_stage_impacts.csv",
    "results/final/figures/source_data/selective_vs_forced.csv",
    "results/final/figures/source_data/simulation_seed_stability.csv",
    "results/final/figures/source_data/reference_split_stability.csv",
    "results/final/figures/source_data/partial_confounding.csv"
  ),
  destination = c(
    file.path("generated/tables", c(
      "main_simulation_rows.tex", "main_mtbls79_rows.tex", "main_added_public_rows.tex",
      "main_clinical_rows.tex", "main_noqc_rows.tex", "full_benchmark_tables.tex",
      "supp_derivative_biology_rows.tex", "ablation_tables.tex", "supp_categorical_r2_trajectory.tex",
      "supp_reference_split_stability.tex", "supp_partial_confounding.tex"
    )),
    file.path("generated/source_data", c(
      "main_simulation_rows.csv", "main_mtbls79_rows.csv", "main_added_public_rows.csv",
      "main_clinical_rows.csv", "main_noqc_rows.csv", "full_benchmark_tables.csv"
    )),
    file.path("generated/source_data", c(
      "supp_derivative_biology_rows.csv", "supp_association_calibration.csv",
      "supp_lobo_summary.csv", "supp_categorical_r2_trajectory.csv",
      "supp_reference_split_stability.csv", "supp_partial_confounding.csv"
    )),
    "generated/source_data/ablation_stage_impacts_table.csv",
    file.path("generated/figures", c(
      "F2_public_data_derivative.pdf", "F3_clinical_pair_aware_derivative.pdf",
      "F4_waveica_tradeoff.pdf", "SF_association_permutation_calibration.pdf",
      "SF_lobo_effect_replication.pdf", "SF_replicate_preservation.pdf",
      "SF_categorical_batch_r2_trajectory.pdf"
    )),
    "generated/figures/SF_ablation_stage_impacts.pdf",
    "generated/figures/SF_selective_vs_forced.pdf",
    "generated/figures/SF_simulation_seed_stability_updated.pdf",
    "generated/figures/SF_reference_split_stability.pdf",
    "generated/figures/SF_partial_confounding.pdf",
    file.path("generated/source_data", c(
      "F2_public_data_derivative.csv", "F3_clinical_pair_aware_derivative.csv",
      "F4_waveica_tradeoff.csv", "SF_association_permutation_calibration.csv",
      "SF_lobo_effect_replication.csv", "SF_replicate_preservation.csv",
      "SF_categorical_batch_r2_trajectory.csv"
    )),
    "generated/source_data/SF_ablation_stage_impacts.csv",
    "generated/source_data/SF_selective_vs_forced.csv",
    "generated/source_data/SF_simulation_seed_stability_updated.csv",
    "generated/source_data/SF_reference_split_stability.csv",
    "generated/source_data/SF_partial_confounding.csv"
  ), stringsAsFactors = FALSE
)

source_paths <- file.path(repo_root, mapping$source)
missing <- mapping$source[!file.exists(source_paths)]
if (length(missing)) stop("Missing derivative artifact(s): ", paste(missing, collapse = ", "), call. = FALSE)
destination_paths <- file.path(manuscript_root, mapping$destination)
for (directory in unique(dirname(destination_paths))) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
copied <- file.copy(source_paths, destination_paths, overwrite = TRUE, copy.mode = TRUE)
if (!all(copied)) stop("One or more manuscript assets could not be staged.", call. = FALSE)

mapping$bytes <- as.numeric(file.info(destination_paths)$size)
mapping$sha256 <- vapply(destination_paths, digest::digest, character(1), file = TRUE, algo = "sha256")
mapping$source_sha256 <- vapply(source_paths, digest::digest, character(1), file = TRUE, algo = "sha256")
mapping$identical_to_source <- mapping$sha256 == mapping$source_sha256
mapping$correction_methods_executed <- FALSE
mapping$tuning_executed <- FALSE
mapping$selection_executed <- FALSE
if (!all(mapping$identical_to_source)) stop("Staged manuscript asset hash mismatch.", call. = FALSE)
write.csv(
  mapping, file.path(manuscript_root, "generated", "DERIVATIVE_ASSET_MANIFEST.csv"),
  row.names = FALSE, na = ""
)
cat("Staged ", nrow(mapping), " derivative manuscript assets with matching SHA-256 values.\n", sep = "")
