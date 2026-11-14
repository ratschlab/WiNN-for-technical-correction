#!/usr/bin/env Rscript
# Independent spot checks against base R, across every task and all three matrices.
args <- commandArgs(TRUE);stopifnot(length(args)==2)
records <- list()
for(task in 1:190) {
  d <- file.path(args[1],sprintf('task_%03d',task))
  stopifnot(file.exists(file.path(d,'complete.json')))
  m <- readRDS(file.path(d,'matrices.rds'))
  a <- read.csv(file.path(d,'allocation.csv'))
  e <- read.csv(gzfile(file.path(d,'feature_effects.csv.gz')))
  stopifnot(identical(dimnames(m$raw),dimnames(m$clean)),
    identical(dimnames(m$raw),dimnames(m$corrected)),
    identical(colnames(m$raw),a$sample_id))
  methods <- list(Raw=m$raw,WiNN=m$corrected,Truth=m$clean)
  worst <- 0
  for(method in names(methods)) {
    z <- methods[[method]];v <- e[e$method==method,]
    stopifnot(identical(rownames(z),v$feature_id),sum(v$responsive)==100,
      all(v$selected==(v$padj<.05)),max(abs(v$padj-p.adjust(v$p,'BH')))<1e-12)
    chosen <- c(head(which(v$responsive),2),tail(which(!v$responsive),2))
    for(j in chosen) {
      values <- log1p(z[j,]);y <- a$phenotype_numeric
      ref <- t.test(values[which(y==1)],values[which(y==0)])
      difference <- max(abs(v$p[j]-ref$p.value),abs(v$effect[j]-diff(rev(ref$estimate))))
      stopifnot(difference<1e-11);worst <- max(worst,difference)
    }
  }
  records[[task]] <- data.frame(task=task,dimnames_match=TRUE,bh_reproduced=TRUE,
    independent_welch_checks=12,max_absolute_difference=worst,passed=TRUE)
}
write.csv(do.call(rbind,records),args[2],row.names=FALSE)
cat('PASS: all 190 matrix/label alignments; BH across all features; 2280 independent Welch comparisons\n')
