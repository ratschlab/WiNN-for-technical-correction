#!/usr/bin/env Rscript
# Targeted, artifact-free stress tests; these do not replace the frozen benchmark.
# Run from the metabolomics_winn workspace root with Rscript <this file>.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
root <- dirname(dirname(script_dir))
out <- file.path(root, "generated", "review_validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
src <- file.path(script_dir, "source_snapshot", "winn")
source(file.path(src, "R", "winn.R"))
source(file.path(src, "R", "winn-ablation.R"))
stopifnot(read.dcf(file.path(src, "DESCRIPTION"))[1, "Version"] == "0.1.6")
write.csv(data.frame(path = file.path("scripts", "targeted_validation", "source_snapshot", "winn", "R", c("winn.R", "winn-ablation.R")),
  md5 = tools::md5sum(c(file.path(src, "R", "winn.R"), file.path(src, "R", "winn-ablation.R")))),
  file.path(out, "biology_source_manifest.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "biology_sessionInfo.txt"))

# A deterministic counterexample copied from the review, on the log1p scale.
z <- rep(1, 100); z[seq(10, 100, 10)] <- 0
x <- matrix(expm1(z), nrow = 1)
prepared <- adjust_outliers_mad(x)
p <- Box.test(z, lag = 10, type = "Ljung-Box", fitdf = 0)$p.value
rho10 <- as.numeric(acf(z, lag.max = 10, plot = FALSE)$acf[11])
counter <- data.frame(n = 100, qc_period = 10, technical_artifacts = FALSE,
  mad = mad(x[1, ]), mad_entries_changed = sum(prepared != x), rho10 = rho10,
  q_lag10_contribution = 100 * 102 * rho10^2 / 90, ljung_box_p = p,
  selected = p < .05)
write.csv(counter, file.path(out, "periodic_qc_counterexample.csv"), row.names = FALSE)
stopifnot(counter$mad_entries_changed == 0, abs(counter$rho10 - .9) < 1e-10,
  counter$selected)

gate_stats <- function(p, truth) {
  reject <- is.finite(p) & p.adjust(p, "BH") < .05
  c(discoveries = sum(reject), sensitivity = if (any(truth)) mean(reject[truth]) else NA_real_,
    false_discovery_proportion = sum(reject & !truth) / max(1, sum(reject)),
    null_p_lt_005 = mean(p[!truth] < .05, na.rm = TRUE),
    null_p_q01 = unname(quantile(p[!truth], .01, na.rm = TRUE)),
    null_p_q05 = unname(quantile(p[!truth], .05, na.rm = TRUE)),
    null_p_median = median(p[!truth], na.rm = TRUE))
}

# Feature-wise random biological noise; scheduled QCs differ in mean but have
# no technical artifacts. Correlated case has a common random component.
grid <- expand.grid(n = c(20L, 50L, 100L, 200L), rho = c(0, .6),
  scenario = c("artifact_free_study_only", "artifact_free_periodic_qc", "mixed_known_drift"),
  seed = seq_len(100), stringsAsFactors = FALSE)
rows <- vector("list", nrow(grid) * 2L); cursor <- 0L
p_hist <- list(); h <- 0L
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]; set.seed(20261027 + g$seed)
  n <- g$n; nf <- 100L; is_qc <- seq_len(n) %% 10L == 0L
  noise <- sqrt(1 - g$rho) * matrix(rnorm(nf * n), nf) +
    sqrt(g$rho) * matrix(rep(rnorm(n), each = nf), nf)
  logx <- 3 + .15 * noise
  truth <- rep(FALSE, nf)
  if (g$scenario == "artifact_free_periodic_qc") logx[, is_qc] <- logx[, is_qc] - 1
  if (g$scenario == "mixed_known_drift") {
    truth[1:20] <- TRUE
    logx[truth, ] <- sweep(logx[truth, , drop = FALSE], 2, seq(-.6, .6, length.out = n), "+")
  }
  mat <- adjust_outliers_mad(expm1(logx))
  for (mode in c("all_injections", "study_only")) {
    keep <- if (mode == "study_only" && g$scenario == "artifact_free_periodic_qc") !is_qc else rep(TRUE, n)
    lag <- .resolve_autocorrelation_lag(NULL, sum(keep), 0L)
    pv <- .compute_ljung_box_pvalues(mat[, keep, drop = FALSE], lag, 0L)
    # Assert vectorized production test agrees with stats::Box.test.
    if (i == 1L) stopifnot(max(abs(pv - apply(log1p(mat[, keep, drop = FALSE]), 1,
      function(v) Box.test(v, lag = lag, type = "Ljung-Box", fitdf = 0)$p.value))) < 1e-12)
    cursor <- cursor + 1L
    rows[[cursor]] <- cbind(g, mode = mode, lag = lag, as.data.frame(as.list(gate_stats(pv, truth))))
    if (g$seed <= 10L) {
      h <- h + 1L
      p_hist[[h]] <- data.frame(g, mode = mode, feature = seq_len(nf),
        technical_drift = truth, p_value = pv, p_adj = p.adjust(pv, "BH"), row.names = NULL)
    }
  }
}
gate <- do.call(rbind, rows[seq_len(cursor)])
write.csv(gate, file.path(out, "drift_gate_calibration_replicates.csv"), row.names = FALSE)
write.csv(do.call(rbind, p_hist), file.path(out, "drift_gate_pvalues_first10seeds.csv"), row.names = FALSE)
gate_summary <- aggregate(gate[c("discoveries", "sensitivity", "false_discovery_proportion",
  "null_p_lt_005", "null_p_q01", "null_p_q05", "null_p_median")],
  gate[c("n", "rho", "scenario", "mode")], function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE))
write.csv(gate_summary, file.path(out, "drift_gate_calibration_summary.csv"), row.names = FALSE)

# Direct batch ANOVA gate, before/after fixed selective drift. Four balanced
# batches; each identity appears once in each batch in the repeated case.
# This explicitly evaluates repeated identity dependence with no batch effect.
batch_rows <- list(); k <- 0L
for (seed in seq_len(100)) for (repeated in c(FALSE, TRUE)) for (artifact in c(FALSE, TRUE)) {
  set.seed(300000 + seed); nf <- 100L; nb <- 4L; per_batch <- 25L
  ids <- if (repeated) rep(seq_len(per_batch), nb) else seq_len(nb * per_batch)
  b <- rep(seq_len(nb), each = per_batch)
  identity_effect <- matrix(rnorm(nf * length(unique(ids)), sd = .3), nf)
  logx <- 3 + identity_effect[, ids] + matrix(rnorm(nf * length(ids), sd = .1), nf)
  truth <- rep(FALSE, nf)
  if (artifact) {
    truth[1:20] <- TRUE
    logx[truth, ] <- sweep(logx[truth, , drop = FALSE], 2, c(-.4, -.1, .1, .4)[b], "+")
  }
  mat <- adjust_outliers_mad(expm1(logx))
  # Production drift decision and correction are executed, not approximated.
  drift <- autocorrelation_correct(mat, batch = b, run_order = seq_along(b),
    gate = "selective", return_diagnostics = TRUE)$data
  for (stage in c("before_drift", "after_drift")) {
    value <- if (stage == "before_drift") mat else drift
    pv <- .compute_batch_anova_pvalues(log1p(value), b)
    k <- k + 1L
    batch_rows[[k]] <- data.frame(seed = seed, repeated_identities = repeated,
      technical_batch_artifact = artifact, stage = stage, as.list(gate_stats(pv, truth)))
  }
}
batch <- do.call(rbind, batch_rows)
write.csv(batch, file.path(out, "batch_gate_calibration_replicates.csv"), row.names = FALSE)
write.csv(aggregate(batch[c("discoveries", "sensitivity", "false_discovery_proportion", "null_p_lt_005")],
  batch[c("repeated_identities", "technical_batch_artifact", "stage")],
  function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)),
  file.path(out, "batch_gate_calibration_summary.csv"), row.names = FALSE)

# Biological stress tests contain no technical artifacts. Half-up/half-down
# is included as the balanced reference. Responses are randomly allocated.
# Effects are compared with the exact clean contrast on the SAME log1p scale.
stages <- c("raw", "mad_only", "pqn_only", "mad_then_pqn")
bio_rows <- list(); k <- 0L
for (seed in seq_len(100)) for (scenario in c("balanced_sparse", "rare_subgroup", "broad_coordinated")) {
  set.seed(500000 + seed); nf <- 200L; ns <- 200L
  response <- rep(0, nf); response[1:40] <- rep(c(.7, -.7), each = 20)
  prevalence <- .5
  if (scenario == "rare_subgroup") { prevalence <- .05; response[1:40] <- 1.4 }
  if (scenario == "broad_coordinated") { prevalence <- .2; response[1:160] <- .7 }
  y <- integer(ns); y[sample.int(ns, round(ns * prevalence))] <- 1L
  latent <- rnorm(ns, sd = .08)
  baseline <- seq(2, 5, length.out = nf)
  logx <- matrix(rnorm(nf * ns, sd = .12), nf) + baseline +
    matrix(rep(latent, each = nf), nf) + outer(response, y)
  clean <- exp(logx)
  mad <- adjust_outliers_mad(clean)
  candidates <- list(raw = clean, mad_only = mad,
    pqn_only = normalize_by_dilution_factor(clean),
    mad_then_pqn = normalize_by_dilution_factor(mad))
  effect <- function(value) rowMeans(log1p(value[, y == 1L, drop = FALSE])) -
    rowMeans(log1p(value[, y == 0L, drop = FALSE]))
  ref <- effect(clean); responsive <- response != 0
  for (stage in stages) {
    value <- candidates[[stage]]; estimated <- effect(value)
    k <- k + 1L
    bio_rows[[k]] <- data.frame(seed = seed, scenario = scenario, stage = stage,
      response_prevalence = prevalence, responsive_features = sum(responsive),
      median_clean_effect = median(abs(ref[responsive])),
      median_recovery_ratio = median(estimated[responsive] / ref[responsive]),
      responsive_sign_concordance = mean(sign(estimated[responsive]) == sign(ref[responsive])),
      null_absolute_effect_change = median(abs(estimated[!responsive] - ref[!responsive])),
      responsive_absolute_effect_change = median(abs(estimated[responsive] - ref[responsive])),
      fraction_entries_changed = mean(abs(value - clean) > 1e-12),
      log1p_rmse_vs_clean = sqrt(mean((log1p(value) - log1p(clean))^2)))
  }
}
bio <- do.call(rbind, bio_rows)
write.csv(bio, file.path(out, "preparation_biology_replicates.csv"), row.names = FALSE)
write.csv(aggregate(bio[c("median_recovery_ratio", "responsive_sign_concordance", "null_absolute_effect_change",
  "responsive_absolute_effect_change", "fraction_entries_changed", "log1p_rmse_vs_clean")],
  bio[c("scenario", "stage")], mean), file.path(out, "preparation_biology_summary.csv"), row.names = FALSE)

cat("Targeted technical-gate and biological-preparation experiments complete.\n")
