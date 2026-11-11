#!/usr/bin/env Rscript
# Re-render the release overview figures from finished aggregate tables.
#
# Used after a WiNN-only re-run has been spliced into a previous release
# (analysis/splice_winn_rerun.R): the aggregate on the cluster only saw the re-run
# methods, so the all-method overview figures must be drawn again from the spliced
# tables. The plotting code is evaluated directly from analysis/aggregate_release.R
# (no copy), so the figures are identical in construction to a full aggregate. No
# correction, tuning or evaluation code is run.
#
# Usage: Rscript analysis/render_release_figures.R --final=<results/final>
# Requires <final>/tables/* and <final>/figures/source_data/ablation_step_impacts.csv.
args <- commandArgs(TRUE)
final_arg <- sub("^--final=", "", grep("^--final=", args, value = TRUE))
if (length(final_arg) != 1L) stop("Supply --final=<results/final>", call. = FALSE)
final_dir <- normalizePath(final_arg, mustWork = TRUE)
file_arg <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE)
aggregate_script <- file.path(dirname(normalizePath(file_arg)), "aggregate_release.R")

exprs <- parse(aggregate_script, keep.source = FALSE)
text <- vapply(exprs, function(e) paste(deparse(e, width.cutoff = 500L), collapse = " "), character(1))
first <- function(pattern) {
  hit <- grep(pattern, text)
  if (!length(hit)) stop("Pattern not found in aggregate_release.R: ", pattern, call. = FALSE)
  hit[1]
}
setup <- c(
  first("^suppressPackageStartupMessages"),
  first("^dataset_order <- "), first("^dataset_labels <- "),
  first("^method_order <- "), first("^method_palette <- ")
)
for (index in setup) eval(exprs[[index]], envir = globalenv())

table_dir <- file.path(final_dir, "tables")
figure_dir <- file.path(final_dir, "figures")
source_dir <- file.path(figure_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
read_table <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
primary <- read_table(file.path(table_dir, "primary_method_metrics_long.csv"))
stage_source <- read_table(file.path(source_dir, "ablation_step_impacts.csv"))
gate_ids <- c("G_SS", "G_AS", "G_SA", "G_AA")
gate_values <- read_table(file.path(table_dir, "selective_vs_forced_metrics_long.csv"))
seed_methods <- read_table(file.path(table_dir, "simulation_seed_method_metrics.csv"))
reference_metrics <- read_table(file.path(table_dir, "reference_split_method_metrics.csv"))
partial_summary <- read_table(file.path(table_dir, "partial_confounding_summary.csv"))

plot_start <- first("^theme_publication <- function")
plot_end <- first("^write.csv\\(figure_manifest, file.path\\(figure_dir")
for (index in seq.int(plot_start, plot_end)) eval(exprs[[index]], envir = globalenv())
message("Release overview figures re-rendered from tables in ", final_dir)
