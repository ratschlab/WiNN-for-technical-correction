#!/usr/bin/env Rscript
# Reproduce the fixed-seed interval and deterministic scheduled-QC example.
a <- commandArgs(TRUE)
stopifnot(length(a) == 2L)
source(file.path(normalizePath(a[1]), "R", "winn.R"))
out <- normalizePath(a[2])
d <- read.csv(file.path(out, "fixed_seed_null_replicates.csv"))
stopifnot(nrow(d) == 1000L, !anyDuplicated(d$seed),
          setequal(d$seed, seq_len(1000L)), all(d$n == 50L),
          all(d$features == 100L), all(d$permutations == 999L),
          all(d$permutation_seed == 104729L),
          all(d$any_rejection == (d$rejected_features > 0L)))
ci <- binom.test(sum(d$any_rejection), nrow(d))$conf.int
write.csv(data.frame(n = 50L, features = 100L, replicates = nrow(d),
  permutation_seed = 104729L, any_false_selection = mean(d$any_rejection),
  lower = ci[1], upper = ci[2]),
  file.path(out, "fixed_seed_null_summary.csv"), row.names = FALSE)

z <- rep(1, 100L)
z[seq(10L, 100L, by = 10L)] <- 0
x <- matrix(expm1(z), nrow = 1L)
prepared <- adjust_outliers_mad(x)
z_prepared <- as.numeric(log1p(prepared))
legacy <- Box.test(z_prepared, lag = 10L, type = "Ljung-Box")
permutation <- .compute_permutation_ljung_box(prepared, 10L)
centered <- z_prepared - mean(z_prepared)
rho10 <- sum(centered[11:100] * centered[1:90]) / sum(centered^2)
write.csv(data.frame(n = 100L, qc_period = 10L, mad = mad(z),
  mad_entries_changed = sum(prepared != x), rho10 = rho10,
  q_statistic = unname(legacy$statistic), asymptotic_p = legacy$p.value,
  permutation_raw_p = permutation$raw_p_value,
  permutation_maxT_p = permutation$adjusted_p_value,
  selected = permutation$adjusted_p_value < 0.05,
  permutations = 999L, permutation_seed = 104729L),
  file.path(out, "periodic_qc_counterexample.csv"), row.names = FALSE)
cat("Reproduced the fixed-seed summary and deterministic periodic-QC example.\n")
