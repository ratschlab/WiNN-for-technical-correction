#!/usr/bin/env Rscript
# Compare function bodies to verify unchanged preprocessing, smoothers and batch/PQN.
a<-commandArgs(TRUE);stopifnot(length(a)==3L)
old<-new.env();updated<-new.env()
for(f in c('winn.R','winn-ablation.R')) {
 source(file.path(a[1],'R',f),local=old)
 source(file.path(a[2],'R',f),local=updated)
}
nms<-intersect(ls(old,all.names=TRUE),ls(updated,all.names=TRUE))
nms<-nms[vapply(nms,function(n)is.function(old[[n]])&&is.function(updated[[n]]),logical(1))]
d<-data.frame(function_name=nms,body_unchanged=vapply(nms,function(n)
 identical(deparse(body(old[[n]])),deparse(body(updated[[n]]))),logical(1)))
allowed<-c('.compute_segment_autocorrelation_pvalues','.prepare_autocorrelation_correction',
 'autocorrelation_correct','winn','winn_ablation')
stopifnot(setequal(d$function_name[!d$body_unchanged],allowed))
write.csv(d,a[3],row.names=FALSE)
cat(sum(d$body_unchanged),'existing function bodies unchanged; only',length(allowed),'gate/call-path bodies changed.\n')
