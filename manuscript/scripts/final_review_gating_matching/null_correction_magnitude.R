#!/usr/bin/env Rscript
# Reuse all 4000 existing Gaussian technical-null seeds, and measure ONLY the
# incremental production drift correction after MAD. No batch or PQN stage.
a<-commandArgs(FALSE); root<-dirname(dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',a,value=TRUE))))))
source(file.path(root,'scripts/targeted_validation/source_snapshot/winn/R/winn.R'))
out<-file.path(root,'generated/final_review_gating_matching');dir.create(out,recursive=TRUE,showWarnings=FALSE)
old<-read.csv(file.path(root,'generated/review_validation/null_gate_verification_replicates.csv'))
old<-old[old$stage=='mad_prepared',]
stopifnot(nrow(old)==4000L,!anyDuplicated(old[c('n','seed')]),
 setequal(unique(old$n),c(20L,50L,100L,200L)),
 all(vapply(c(20L,50L,100L,200L),function(n)identical(sort(as.integer(old$seed[old$n==n])),seq_len(1000L)),logical(1))))
rows<-vector('list',nrow(old))
for(i in seq_len(nrow(old))) {
 n<-old$n[i]; seed<-old$seed[i];set.seed(900000+10000L*n+seed)
 raw<-expm1(3+matrix(rnorm(100L*n,sd=.15),100L)); prepared<-adjust_outliers_mad(raw)
 # Isolate the incremental production drift stage on the prepared null matrix.
 vals<-lapply('selective',function(g) autocorrelation_correct(prepared,run_order=seq_len(n),batch=rep(1,n),
   lag=NULL,test='Ljung-Box',ljung_box_fitdf=0L,detrend='mean',fdr_threshold=.05,spline_method='conservative',gate=g,return_diagnostics=TRUE))
 rows[[i]]<-do.call(rbind,lapply(seq_along(vals),function(j) {
  diagnostic<-vals[[j]]$diagnostics
  stopifnot(nrow(diagnostic)==100L,all(is.finite(diagnostic$adjusted_p_value)),
    !anyNA(diagnostic$selective_significance))
  selected_features<-sum(diagnostic$selective_significance)
  stopifnot(selected_features==old$rejected_features[i])
  d<-log1p(vals[[j]]$data)-log1p(prepared)
  stopifnot(all(is.finite(d)))
  data.frame(n=n,seed=seed,policy='selective',selected_features=selected_features,
   changed_features=sum(rowSums(abs(d)>1e-12)>0),incremental_log1p_rmse=sqrt(mean(d^2)),
   mean_absolute_log1p_change=mean(abs(d)),max_absolute_log1p_change=max(abs(d)),
   noise_sd=.15)
 }))
 if(i%%100==0)cat(i,'/',nrow(old),'\n')
}
z<-do.call(rbind,rows);write.csv(z,file.path(out,'null_correction_magnitude_replicates.csv'),row.names=FALSE)
source(file.path(root,'scripts/final_review_gating_matching/summarize_null_magnitude.R'))
