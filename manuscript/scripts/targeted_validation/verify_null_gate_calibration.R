#!/usr/bin/env Rscript
# Independent, larger null check of the exact production Ljung-Box test.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
root <- dirname(dirname(script_dir))
out <- file.path(root, "generated", "review_validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "source_snapshot", "winn", "R", "winn.R"))
rows <- list(); k <- 0L
for (n in c(20L, 50L, 100L, 200L)) for (seed in seq_len(1000)) {
  set.seed(900000 + 10000L * n + seed)
  z <- 3 + matrix(rnorm(100L * n, sd = .15), 100L)
  raw <- expm1(z)
  prepared <- adjust_outliers_mad(raw)
  lag <- .resolve_autocorrelation_lag(NULL, n, 0L)
  for (stage in c("unprepared", "mad_prepared")) {
    x <- if (stage == "unprepared") raw else prepared
    p <- .compute_ljung_box_pvalues(x, lag, 0L)
    # Independent base-R reference for every feature of the first 10 seeds.
    if (seed <= 10L) {
      reference <- apply(log1p(x), 1, function(v) Box.test(v, lag = lag, type = "Ljung-Box", fitdf = 0)$p.value)
      stopifnot(max(abs(p - reference)) < 1e-12)
    }
    rejected <- p.adjust(p, "BH") < .05
    k <- k + 1L
    rows[[k]] <- data.frame(n = n, lag = lag, seed = seed, stage = stage,
      rejected_features = sum(rejected), any_rejection = any(rejected),
      null_p_lt_005 = mean(p < .05), null_p_lt_0005 = mean(p < .0005))
  }
}
results <- do.call(rbind, rows)
write.csv(results, file.path(out, "null_gate_verification_replicates.csv"), row.names = FALSE)
summary <- aggregate(results[c("rejected_features", "any_rejection", "null_p_lt_005", "null_p_lt_0005")],
  results[c("n", "lag", "stage")], mean)
summary$replicates <- 1000L
interval <- t(apply(summary, 1, function(row) {
  binom.test(round(as.numeric(row[["any_rejection"]]) * 1000), 1000)$conf.int
}))
summary$any_rejection_lower <- interval[, 1]; summary$any_rejection_upper <- interval[, 2]
write.csv(summary, file.path(out, "null_gate_verification_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
