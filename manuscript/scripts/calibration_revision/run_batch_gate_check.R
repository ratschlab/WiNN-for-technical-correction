#!/usr/bin/env Rscript
# Repeat only the batch stress test whose after-drift input changes.
a<-commandArgs(TRUE);stopifnot(length(a)>=2L)
source(file.path(normalizePath(a[1]),'R','winn.R'))
shard<-if(length(a)>=3L)as.integer(a[3]) else 1L
shards<-if(length(a)>=4L)as.integer(a[4]) else 1L
out<-a[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
rows<-list(); k<-0L
for(seed in seq_len(100L)) {
 if((seed-1L)%%shards+1L != shard)next
 for(repeated in c(FALSE,TRUE))for(artifact in c(FALSE,TRUE)) {
  set.seed(300000+seed);nf<-100L;nb<-4L;per_batch<-25L
  ids<-if(repeated)rep(seq_len(per_batch),nb) else seq_len(nb*per_batch)
  b<-rep(seq_len(nb),each=per_batch)
  identity_effect<-matrix(rnorm(nf*length(unique(ids)),sd=.3),nf)
  logx<-3+identity_effect[,ids]+matrix(rnorm(nf*length(ids),sd=.1),nf)
  truth<-rep(FALSE,nf)
  if(artifact) {
   truth[1:20]<-TRUE
   logx[truth,]<-sweep(logx[truth,,drop=FALSE],2,c(-.4,-.1,.1,.4)[b],'+')
  }
  mat<-adjust_outliers_mad(expm1(logx))
  drift<-autocorrelation_correct(mat,batch=b,run_order=seq_along(b),
    gate='selective',return_diagnostics=TRUE)$data
  for(stage in c('before_drift','after_drift')) {
   value<-if(stage=='before_drift')mat else drift
   pv<-.compute_batch_anova_pvalues(log1p(value),b)
   reject<-is.finite(pv)&p.adjust(pv,'BH')<.05;k<-k+1L
   rows[[k]]<-data.frame(seed=seed,repeated_identities=repeated,
     technical_batch_artifact=artifact,stage=stage,
     discoveries=sum(reject),sensitivity=if(any(truth))mean(reject[truth])else NA_real_,
     false_discovery_proportion=sum(reject&!truth)/max(1,sum(reject)),
     null_p_lt_005=mean(pv[!truth]<.05,na.rm=TRUE))
  }
 }
}
write.csv(do.call(rbind,rows),file.path(out,sprintf('batch_gate_%03d.csv',shard)),row.names=FALSE)

capture.output(sessionInfo(),file=file.path(out,sprintf('batch_session_%03d.txt',shard)))
