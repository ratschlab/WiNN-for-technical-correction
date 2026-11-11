#!/usr/bin/env Rscript
# Combine the release validation of a WiNN-only re-run with the validation of the
# previous release whose competitor outputs it reuses.
#
# validate_release.R run with WINN_RERUN_METHODS audits the artifacts executed in the
# re-run (Raw and the WiNN modes) plus the spliced final tables. Checks whose detail
# counts executed artifacts (tasks, ledgers, failures, run manifests, matrices,
# reference separation, batch R2 and replicate-metric rows) therefore describe only
# the re-run. For those checks this script appends the previous release's detail for
# the reused competitor outputs, and requires both validations to have passed.
#
# Usage: Rscript analysis/combine_winn_rerun_validation.R --rerun=<release_validation.csv> \
#   --previous=<previous release_validation.csv> --out=<combined release_validation.csv>
args <- commandArgs(TRUE)
arg <- function(name) {
  hit <- sub(paste0("^--", name, "="), "", grep(paste0("^--", name, "="), args, value = TRUE))
  if (length(hit) != 1L) stop("Supply --", name, "=", call. = FALSE)
  hit
}
rerun <- read.csv(arg("rerun"), stringsAsFactors = FALSE)
previous <- read.csv(arg("previous"), stringsAsFactors = FALSE)
if (!all(rerun$passed) || !all(previous$passed)) stop("Both validations must pass.", call. = FALSE)
if (!setequal(rerun$check, previous$check)) stop("The two validations audit different checks.", call. = FALSE)
artifact_checks <- c(
  "all_prespecified_tasks_complete", "execution_ledgers_cover_prespecified_tasks_and_attempts",
  "no_unresolved_failure_artifacts", "all_analytical_runs_use_frozen_winn",
  "all_corrected_matrices_are_finite_and_nonnegative", "training_and_hidden_references_are_separated",
  "batch_r2_is_categorical", "replicate_metrics_have_distinct_definitions"
)
combined <- rerun
index <- combined$check %in% artifact_checks
old <- previous$detail[match(combined$check[index], previous$check)]
combined$detail[index] <- paste0(
  "WiNN 0.1.6 re-run: ", sub("^re-run methods \\([^)]*\\): ", "", combined$detail[index]),
  ". Reused competitor outputs (previous validated release): ", old
)
write.csv(combined, arg("out"), row.names = FALSE)
cat("Combined", nrow(combined), "checks;", sum(index), "details include the reused-output validation.\n")
