#!/usr/bin/env Rscript
# Build the revised RFSC/baseline input bundles from the public download pipeline.
# This creates local matrices; generated inputs must not be committed.
args <- commandArgs(TRUE)
arg <- function(key, default) {
  v <- args[startsWith(args, paste0(key, "="))]
  if (length(v)) sub(paste0("^", key, "="), "", v[1]) else default
}
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else commandArgs(FALSE)[match("-f", commandArgs(FALSE)) + 1L]
root <- dirname(dirname(normalizePath(script)))
source_root <- normalizePath(arg("--source-root", root), mustWork = TRUE)
out <- arg("--output-root", file.path(root, "results", "revision_inputs"))
requested <- arg("--dataset", "public")
datasets <- if (requested == "public") c(sprintf("SIM%02d", 1:10), "mtbls79", "batchcorr_set1", "sacurine", "waveica") else strsplit(requested, ",", fixed = TRUE)[[1]]
allowed <- c(sprintf("SIM%02d", 1:10), "mtbls79", "batchcorr_set1", "sacurine", "waveica", "clinical_fiams")
stopifnot(all(datasets %in% allowed))
source(file.path(root, "clinical_comparison_helpers.R"))
source(file.path(root, "scripts", "winn_ablation_helpers.R"))
source(file.path(root, "manuscript", "scripts", "benchmark_revision", "training_qc_tuning_engine.R"))
config <- file.path(root, "analysis", "config", "revision_20261027")
tq_output_root <- function(repo_root) config
tq_registry_path <- function(repo_root) file.path(root, "analysis", "config", "training_qc_candidate_registry.csv")
tq_fold_manifest_path <- function(repo_root) {
  if (dataset == "clinical_fiams") {
    path <- arg("--clinical-fold-manifest", "")
    if (!nzchar(path)) stop("Clinical inputs require the authorized private --clinical-fold-manifest path.")
    return(path)
  }
  file.path(config, "internal_qc_fold_manifest.csv")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- list()
for (dataset in datasets) {
  if (grepl("^SIM", dataset)) {
    bundle <- file.path(source_root, "data", "simulated", "canonical", dataset)
    x <- readRDS(file.path(bundle, "raw_intensity.rds"))
    truth <- readRDS(file.path(bundle, "clean_ground_truth.rds"))
    meta <- read.csv(file.path(bundle, "sample_metadata.csv"), check.names = FALSE)
    meta$is_study <- !as.logical(meta$is_qc)
    meta$biological_label <- ifelse(meta$is_study, "study", "QC")
    p <- list(x = x, truth = truth, meta = meta,
      training_ids = as.character(meta$sample_id[as.logical(meta$is_qc) & !as.logical(meta$is_hidden_reference)]),
      hidden_ids = as.character(meta$sample_id[as.logical(meta$is_hidden_reference)]), external_ids = character())
  } else {
    p <- tq_prepare_dataset(source_root, dataset)
    p$folds <- tq_fold_definitions(root, dataset, p$training_ids, p$meta)
    p$tuning_feature_ids <- tq_tuning_features(p$x, dataset)
    p$registry <- tq_read_registry(root, dataset, "QC-RFSC")
  }
  stopifnot(identical(colnames(p$x), as.character(p$meta$sample_id)), all(is.finite(p$x)), all(p$x >= 0),
    !length(intersect(p$training_ids, c(p$hidden_ids, p$external_ids))))
  saveRDS(p, file.path(out, paste0("input_", dataset, ".rds")))
  checks[[dataset]] <- data.frame(dataset = dataset, features = nrow(p$x), units = ncol(p$x),
    training = length(p$training_ids), withheld = length(c(p$hidden_ids, p$external_ids)),
    matrix_sha256 = digest::digest(p$x, algo = "sha256"), stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, checks), file.path(out, "input_manifest.csv"), row.names = FALSE)
cat("Prepared", length(datasets), "input bundles in", out, "\n")
