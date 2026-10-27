#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
seed_id <- sub("^--seed-id=", "", args[startsWith(args, "--seed-id=")][1])
stopifnot(seed_id %in% sprintf("SIM%02d",1:10))
f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script <- if(length(f)) sub("^--file=","",f[1]) else commandArgs(FALSE)[match('-f',commandArgs(FALSE))+1L]
here <- dirname(normalizePath(script))
arg_value <- function(key, default) {
  value <- args[startsWith(args, paste0(key, '='))]
  if(length(value)) sub(paste0('^', key, '='), '', value[1]) else default
}
input_root <- normalizePath(arg_value('--input-root', here), mustWork=TRUE)
matrix_root <- normalizePath(arg_value('--matrix-root', file.path(here, 'simulation_results')), mustWork=TRUE)
output_root <- arg_value('--output-root', file.path(here, 'simulation_evaluation'))
source_root <- here
for (name in c("original_benchmark_helpers.R","frozen_winn_ablation_helpers.R","frozen_weighted_pc_r2.R","frozen_validate_output.R","frozen_evaluation_helpers.R")) source(file.path(here,name))
suppressPackageStartupMessages({library(dplyr);library(tibble)})
p <- readRDS(file.path(input_root,paste0("input_",seed_id,".rds")))
x <- readRDS(file.path(matrix_root,seed_id,"corrected_training_qc_abundance.rds"))
y <- evaluate_release_matrix("simulation","qc_rfsc","QC-RFSC",x,p)
out <- file.path(output_root,seed_id)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
y$metrics$seed_id <- seed_id
write.csv(y$metrics,file.path(out,"metrics.csv"),row.names=FALSE)
writeLines("status=complete",file.path(out,"COMPLETE.txt"))
