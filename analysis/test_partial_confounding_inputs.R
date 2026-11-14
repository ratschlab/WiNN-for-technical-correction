#!/usr/bin/env Rscript
# Reconstruct the publication grid and test the pre-fit portability contract.
args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value=TRUE))))
root <- dirname(here)
for (name in c("simulation_core.R", "partial_confounding_core.R", "canonical_cache.R", "partial_confounding_portable_hash.R", "partial_confounding_pilot_metrics.R"))
  source(file.path(root, "scripts", "robustness", name))
manifest <- read.csv(file.path(here, "config", "partial_confounding_task_manifest.csv"))
scenario_root <- file.path(root, "data", "simulated", "partial_confounding")
hashes <- read.csv(file.path(scenario_root, "config", "portable_matrix_hashes.csv"))
stopifnot(nrow(manifest) == 160L, length(unique(manifest$seed_id)) == 10L)
for (id in manifest$scenario_id) {
  bundle <- readRDS(file.path(scenario_root, "bundles", paste0(id, ".rds")))
  base <- generate_canonical_simulation(bundle$generation_seeds, include_artifact_matrices=TRUE)
  allocation <- bundle$evaluation_truth$phenotype_allocation
  truth <- bundle$evaluation_truth$feature_truth
  reconstructed <- construct_partial_confounding_matrices(base, allocation, truth, include_intensities=TRUE)
  row <- hashes[hashes$scenario_id == id, , drop=FALSE]
  stopifnot(partial_confounding_verify_frozen_inputs(reconstructed, row))
  study <- allocation[allocation$sample_type == "sample", ]
  stopifnot(nrow(study) == 450L, sum(study$phenotype_numeric) == 225L,
            sum(truth$responsive) == 100L,
            all(abs(truth$phenotype_effect_log[truth$responsive]) == .35))
  if (grepl("batch_1p00$", id))
    stopifnot(identical(as.integer(tapply(study$phenotype_numeric, study$plate, sum)), c(90L,90L,45L,0L,0L)))
}
# A changed diagnostic intensity hash must not reject identical log inputs.
diagnostic <- row
diagnostic$raw_intensity_quantized_9_sha256_diagnostic <- strrep("0", 64)
stopifnot(partial_confounding_verify_frozen_inputs(reconstructed, diagnostic))
# Exact-column and shape checks prevent silent R partial matching or recycling.
bad <- row; bad$raw_log_with_biology_quantized_9_sha256 <- strrep("0", 64)
stopifnot(!partial_confounding_verify_frozen_inputs(reconstructed, bad))
bad <- row; names(bad)[names(bad) == "raw_log_with_biology_quantized_9_sha256"] <- "raw_log_with_biology_quantized_9_sha256_diagnostic"
stopifnot(!partial_confounding_verify_frozen_inputs(reconstructed, bad),
          !partial_confounding_verify_frozen_inputs(reconstructed, rbind(row,row)))
# The phenotype coefficient agrees with ordinary lm, independently of the
# vectorized normal-equation implementation, for responsive and null features.
fit <- fit_partial_confounding_effects(reconstructed$clean_ground_truth, allocation, truth)
for (feature in c(which(truth$responsive)[1:3], which(!truth$responsive)[1:3])) {
  study$response <- log1p(reconstructed$clean_ground_truth[feature, study$sample_id])
  reference <- lm(response ~ phenotype_numeric + factor(plate) + within_plate_position_scaled, data=study)
  stopifnot(abs(coef(reference)[["phenotype_numeric"]] - fit$estimates$recovered_effect_log1p[feature]) < 1e-12)
}
cat("PASS: all 160 publication scenarios reproduce authoritative log hashes; balanced design, exact hash fields, diagnostic-only intensities, and independent effect estimates verified.\n")
