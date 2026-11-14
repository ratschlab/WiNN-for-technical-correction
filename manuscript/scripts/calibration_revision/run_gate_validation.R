#!/usr/bin/env Rscript
# Locked paired comparison of legacy and permutation drift gates.
# Usage: Rscript run_gate_validation.R PACKAGE_SOURCE OUTPUT_DIR [SHARD=1] [SHARDS=1]
# Each shard is independent; merge only after every planned row is present.
a <- commandArgs(TRUE)
stopifnot(length(a) >= 2L)
pkg <- normalizePath(a[1]); out <- a[2]
shard <- if (length(a) >= 3L) as.integer(a[3]) else 1L
shards <- if (length(a) >= 4L) as.integer(a[4]) else 1L
stopifnot(shard >= 1L, shards >= shard)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(pkg, 'R', 'winn.R'))
stopifnot(exists('.compute_permutation_ljung_box'))
# Original null seeds and distributions are preserved exactly.
original <- expand.grid(n = c(20L, 50L, 100L, 200L), seed = seq_len(1000),
                        stringsAsFactors = FALSE)
original$experiment <- 'original_null'; original$rho <- 0
original$distribution <- 'gaussian'; original$scenario <- 'global_null'
original$features <- 100L
expanded <- expand.grid(n = c(20L, 50L, 100L, 200L), rho = c(0, .6),
  distribution = c('gaussian', 'student_t3', 'skewed'),
  scenario = c('global_null', 'mixed_linear_drift', 'periodic_qc'),
  seed = seq_len(100), stringsAsFactors = FALSE)
expanded$experiment <- 'robustness'; expanded$features <- 100L
large <- expand.grid(n = c(50L, 200L), rho = c(0, .6), seed = seq_len(100),
                     stringsAsFactors = FALSE)
large$experiment <- 'larger_family'; large$distribution <- 'gaussian'
large$scenario <- 'global_null'; large$features <- 1000L
cols <- c('experiment','n','rho','distribution','scenario','seed','features')
grid <- rbind(original[cols], expanded[cols], large[cols]); grid$task <- seq_len(nrow(grid))
write.csv(grid, file.path(out, sprintf('design_%03d.csv', shard)), row.names = FALSE)
selected <- which((grid$task - 1L) %% shards + 1L == shard)
metrics <- function(p, truth) {
  reject <- is.finite(p) & p < .05
  c(discoveries = sum(reject), false_discoveries = sum(reject & !truth),
    any_false_selection = any(reject & !truth),
    false_discovery_proportion = sum(reject & !truth)/max(1, sum(reject)),
    sensitivity = if (any(truth)) mean(reject[truth]) else NA_real_)
}
rows <- list(); cursor <- 0L
for (i in selected) {
  g <- grid[i, ]; n <- g$n; nf <- g$features
  data_seed <- if (g$experiment == 'original_null') 900000 + 10000L*n + g$seed else 7000000L + i
  set.seed(data_seed)
  noise <- function(k) switch(g$distribution,
    gaussian = rnorm(k), student_t3 = rt(k, df = 3)/sqrt(3),
    skewed = (rlnorm(k, sdlog = .7) - exp(.7^2/2))/sqrt((exp(.7^2)-1)*exp(.7^2)))
  if (g$experiment == 'original_null') {
    z <- 3 + matrix(rnorm(nf*n, sd=.15), nf)
  } else {
    z <- 3 + .15 * (sqrt(1-g$rho)*matrix(noise(nf*n), nf) +
      sqrt(g$rho)*matrix(rep(noise(n), each=nf), nf))
  }
  truth <- rep(FALSE, nf)
  if (g$scenario == 'mixed_linear_drift') {
    truth[seq_len(nf/5)] <- TRUE
    z[truth, ] <- sweep(z[truth, , drop=FALSE], 2, seq(-.6,.6,length.out=n), '+')
  }
  if (g$scenario == 'periodic_qc') z[, seq_len(n) %% 10L == 0L] <- z[, seq_len(n) %% 10L == 0L] - 1
  # Flooring ensures valid nonnegative intensities, identically across columns.
  # It does not alter the original Gaussian experiment in practice.
  x <- expm1(pmax(z, 0)); dim(x) <- c(nf,n)
  stages <- if (g$experiment == 'original_null') c('unprepared','mad_prepared') else 'mad_prepared'
  for (stage in stages) {
    prepared <- if (stage == 'mad_prepared') adjust_outliers_mad(x) else x
    lag <- .resolve_autocorrelation_lag(NULL,n,0L)
    legacy <- p.adjust(.compute_ljung_box_pvalues(prepared,lag,0L),'BH')
    # Independent, prespecified permutation seed for each generated dataset.
    perm_seed <- 104729L + i
    timing <- system.time(perm <- .compute_permutation_ljung_box(prepared,lag,
      n_permutations=999L, permutation_seed=perm_seed))
    for (algorithm in c('asymptotic_BH','permutation_maxT')) {
      adj <- if (algorithm == 'asymptotic_BH') legacy else perm$adjusted_p_value
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(g, stage=stage, algorithm=algorithm,
        data_seed=data_seed, permutation_seed=perm_seed, permutations=999L,
        lag=lag, as.list(metrics(adj,truth)),
        elapsed_permutation_seconds=unname(timing[['elapsed']]), row.names=NULL)
    }
  }
  if (length(rows) %% 40L == 0L) {
    saveRDS(rows,file.path(out,sprintf('checkpoint_%03d.rds',shard)))
    cat(sprintf('shard %d: task %d/%d\n',shard,i,nrow(grid))); flush.console()
  }
}
write.csv(do.call(rbind,rows),file.path(out,sprintf('replicates_%03d.csv',shard)),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,sprintf('session_%03d.txt',shard)))
