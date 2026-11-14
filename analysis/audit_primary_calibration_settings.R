#!/usr/bin/env Rscript
args<-commandArgs(TRUE);stopifnot(length(args)==1L);root<-normalizePath(args[1])
paths<-list.files(file.path(root,"results/primary"),pattern="^method_details[.]rds$",recursive=TRUE,full.names=TRUE)
paths<-paths[grepl("/winn_[^/]+/",paths)]
stopifnot(length(paths)==18L)
rows<-lapply(paths,function(p){
 x<-readRDS(p)$selected_parameters
 threshold<-function(name)if(is.null(x[[name]]))x$fdr_threshold else x[[name]]
 data.frame(dataset=basename(dirname(dirname(p))),method_id=basename(dirname(p)),
  test=x$test,acorr_fdr=as.numeric(threshold("acorr_fdr")),anova_fdr=as.numeric(threshold("anova_fdr")),
  drift_calibration=x$drift_calibration,n_permutations=x$n_permutations,permutation_seed=x$permutation_seed)
})
x<-do.call(rbind,rows)
stopifnot(!anyDuplicated(x[c("dataset","method_id")]),all(x$test=="Ljung-Box"),
 all(x$drift_calibration=="permutation"),all(x$n_permutations==999L),all(x$permutation_seed==104729L),
 all(is.finite(x$acorr_fdr)),all(x$acorr_fdr%in%c(.01,.05,.1)),
 all(is.finite(x$anova_fdr)),all(x$anova_fdr%in%c(.01,.05,.1)))
out<-file.path(root,"results/final/validation");dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(x,file.path(out,"primary_selected_parameter_preflight.csv"),row.names=FALSE)
cat("PASS: all 18 primary modes retain calibrated settings and actual selected thresholds.\n")
print(x)
