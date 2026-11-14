#!/usr/bin/env Rscript
# Merge complete shards, assert the prespecified design, and plot paired results.
a <- commandArgs(TRUE); stopifnot(length(a)>=1L)
out <- normalizePath(a[1]); files <- list.files(out, '^replicates_[0-9]+[.]csv$', full.names=TRUE)
stopifnot(length(files)>0L)
d <- do.call(rbind,lapply(files,read.csv,stringsAsFactors=FALSE))
keys <- c('task','stage','algorithm')
stopifnot(!anyDuplicated(d[keys]))
design <- read.csv(file.path(out,'design_001.csv'),stringsAsFactors=FALSE)
expected <- sum(ifelse(design$experiment=='original_null',4L,2L))
stopifnot(nrow(d)==expected, setequal(d$task,design$task), all(d$permutations==999L))
d <- d[order(d$task,d$stage,d$algorithm),]
# Periodic QC shifts are ordered signals, although instrumental drift is absent.
# Their selections measure technical-effect specificity, not type-I error for
# the exchangeability null. Keep that distinction explicit in the source data.
d$null_model_applicable <- d$scenario != 'periodic_qc'
write.csv(d,file.path(out,'gate_validation_replicates.csv'),row.names=FALSE)
groups <- c('experiment','n','rho','distribution','scenario','features','stage','algorithm')
f <- interaction(d[groups],drop=TRUE,lex.order=TRUE)
summary <- do.call(rbind,lapply(split(d,f),function(x) {
  ci <- binom.test(sum(x$any_false_selection),nrow(x))$conf.int
  data.frame(x[1,groups],replicates=nrow(x),
    null_model_applicable=x$null_model_applicable[1],
    any_selection_rate=mean(x$discoveries>0),
    technical_false_selection_rate=mean(x$any_false_selection),
    mean_technical_fdp=mean(x$false_discovery_proportion),
    fwer=if(x$null_model_applicable[1])mean(x$any_false_selection)else NA_real_,
    fwer_lower=if(x$null_model_applicable[1])ci[1]else NA_real_,
    fwer_upper=if(x$null_model_applicable[1])ci[2]else NA_real_,
    mean_fdp=if(x$null_model_applicable[1])mean(x$false_discovery_proportion)else NA_real_,
    mean_discoveries=mean(x$discoveries),
    sensitivity=if(all(is.na(x$sensitivity))) NA_real_ else mean(x$sensitivity),
    sensitivity_se=if(all(is.na(x$sensitivity))) NA_real_ else sd(x$sensitivity)/sqrt(sum(is.finite(x$sensitivity))),row.names=NULL)
}))
write.csv(summary,file.path(out,'gate_validation_summary.csv'),row.names=FALSE)
# Reconstruct every original result from the old archived experiment.
old_path <- if(length(a)>=2L) a[2] else NA_character_
if(!is.na(old_path)) {
  old <- read.csv(old_path,stringsAsFactors=FALSE)
  legacy <- d[d$experiment=='original_null' & d$algorithm=='asymptotic_BH',]
  check <- merge(old,legacy,by=c('n','seed','stage'),suffixes=c('_old','_new'))
  stopifnot(nrow(check)==8000L, all(check$rejected_features==check$discoveries),
    all(check$any_rejection==check$any_false_selection))
}
cat('Verified complete prespecified design:',nrow(design),'datasets;',nrow(d),'paired stage/algorithm rows.\n')
print(summary[summary$experiment=='original_null' & summary$stage=='mad_prepared',],row.names=FALSE)
