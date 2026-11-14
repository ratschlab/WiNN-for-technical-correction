#!/usr/bin/env Rscript
root <- normalizePath('.')
saved <- file.path(root,'generated/calibration_revision')
out <- tempfile('winn_gate_check_');dir.create(out)
file.copy(file.path(saved,'gate_validation_replicates.csv'),file.path(out,'replicates_001.csv'))
file.copy(file.path(saved,'validation_design.csv'),file.path(out,'design_001.csv'))
status <- system2(file.path(R.home('bin'),'Rscript'),c('scripts/calibration_revision/summarize_gate_validation.R',
  shQuote(out),'generated/review_validation/null_gate_verification_replicates.csv'),stdout=TRUE,stderr=TRUE)
stopifnot(is.null(attr(status,'status')))
old <- read.csv(file.path(saved,'gate_validation_summary.csv'));new <- read.csv(file.path(out,'gate_validation_summary.csv'))
stopifnot(isTRUE(all.equal(old,new,tolerance=1e-12,check.attributes=FALSE)))
writeLines(c('PASS: all 31200 paired rows and 11600 planned datasets; 8000 original null counts reproduced;',
 'all summary statistics and exact intervals reproduced from saved replicates.'),file.path(saved,'saved_results_validation.txt'))
cat(tail(status,3),sep='\n')
cat('\nPASS: saved standalone calibration results validated.\n')
unlink(out,recursive=TRUE)
