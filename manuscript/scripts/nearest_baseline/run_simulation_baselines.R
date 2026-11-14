#!/usr/bin/env Rscript
# Same-information study-sample baselines on the exact canonical simulation inputs.
args <- commandArgs(TRUE)
opt <- function(k, default=NULL) { z <- args[startsWith(args,paste0(k,'='))]; if(length(z)) sub(paste0('^',k,'='),'',z[1]) else default }
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
here <- dirname(normalizePath(script)); root <- dirname(dirname(here))
input <- normalizePath(opt('--input-root'),mustWork=TRUE)
out <- opt('--output-root',file.path(root,'generated/review_validation'));dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(here,'doBC.R'))
source(file.path(root,'scripts/benchmark_revision/frozen_winn_ablation_helpers.R'))
hashes <- read.csv(file.path(root,'scripts/benchmark_revision/simulation_bundle_hashes.csv'))
reference <- read.csv(file.path(root,'generated/revision_panel_inputs/simulation_seed_method_metrics.csv'))
run_seed <- function(seed) {
 p <- readRDS(file.path(input,paste0('input_',seed,'.rds')))
 checks <- do.call(rbind,lapply(c('raw_intensity','clean_ground_truth'),function(object) {
  value <- if(object=='raw_intensity')p$x else p$truth
  expected <- hashes$reference_object_sha256[hashes$seed_id==seed & hashes$object==object]
  observed <- digest::digest(value,algo='sha256')
  data.frame(seed_id=seed,object=object,expected_sha256=expected,observed_sha256=observed,passed=identical(expected,observed))
 }))
 stopifnot(all(checks$passed))
 md<-p$meta;study<-md$sample_id[!md$is_qc];truth<-log1p(p$truth[,study,drop=FALSE]);z<-log1p(p$x)
 candidates <- list(raw=p$x)
 for(method in c('mean_batch_all','wehrens_ancova_all')) {
  zz <- t(apply(z,1,doBC,ref.idx=seq_len(ncol(z)),batch.idx=md$batch,
    seq.idx=if(method=='mean_batch_all')NULL else md$run_order,method='lm',minBsamp=if(method=='mean_batch_all')2 else 4))
  stopifnot(all(is.finite(zz)))
  x<-pmax(expm1(zz),0);dimnames(x)<-dimnames(p$x);attr(x,'floored_cells')<-sum(zz<0)
  candidates[[method]]<-x
 }
 rows<-do.call(rbind,lapply(names(candidates),function(method) {
  x<-candidates[[method]];value<-log1p(x[,study,drop=FALSE]);qc<-x[,p$hidden_ids,drop=FALSE]
  data.frame(seed_id=seed,method_id=method,heldout_qc_cv_mean=100*mean(apply(qc,1,function(v)sd(v)/abs(mean(v))),na.rm=TRUE),
   truth_log1p_rmse=sqrt(mean((value-truth)^2)),
   truth_metabolite_profile_icc_mean=mean(vapply(seq_len(nrow(value)),function(i)icc_a1_ablation(cbind(truth[i,],value[i,])),numeric(1)),na.rm=TRUE),
   truth_sample_profile_pearson_mean=mean(vapply(seq_len(ncol(value)),function(j)safe_cor_ablation(truth[,j],value[,j]),numeric(1)),na.rm=TRUE),
   floored_cells=if(is.null(attr(x,'floored_cells')))0 else attr(x,'floored_cells'))
 }))
 for(metric in c('heldout_qc_cv_mean','truth_log1p_rmse','truth_metabolite_profile_icc_mean','truth_sample_profile_pearson_mean')) {
  expected<-reference$value[reference$seed_id==seed & reference$method_id=='raw' & reference$metric==metric]
  stopifnot(length(expected)==1,abs(rows[rows$method_id=='raw',metric]-expected)<1e-10)
 }
 list(metrics=rows,checks=checks)
}
# --reuse-baseline-rows=<dir>: reuse the Raw and baseline per-seed rows and input
# checks already computed in <dir> (they do not depend on WiNN) and refresh only the
# WiNN rows and WiNN-versus-baseline comparisons from the reference table.
reuse<-opt('--reuse-baseline-rows')
if(!is.null(reuse)) {
 previous<-read.csv(file.path(reuse,'baseline_simulation_seed_metrics.csv'),stringsAsFactors=FALSE)
 rows<-previous[previous$method_id!='winn_fixed_no_qc',,drop=FALSE]
 checks<-read.csv(file.path(reuse,'baseline_simulation_input_checks.csv'),stringsAsFactors=FALSE)
 stopifnot(all(checks$passed),setequal(unique(rows$method_id),c('raw','mean_batch_all','wehrens_ancova_all')),nrow(rows)==30L)
 for(seed in sprintf('SIM%02d',1:10))for(metric in c('heldout_qc_cv_mean','truth_log1p_rmse','truth_metabolite_profile_icc_mean','truth_sample_profile_pearson_mean')) {
  expected<-reference$value[reference$seed_id==seed & reference$method_id=='raw' & reference$metric==metric]
  stopifnot(length(expected)==1,abs(rows[rows$seed_id==seed & rows$method_id=='raw',metric]-expected)<1e-10)
 }
} else {
cores<-as.integer(opt('--cores','4'))
result<-parallel::mclapply(sprintf('SIM%02d',1:10),run_seed,mc.cores=cores)
stopifnot(all(vapply(result,is.list,logical(1))))
rows<-do.call(rbind,lapply(result,'[[','metrics'));checks<-do.call(rbind,lapply(result,'[[','checks'))
}
# Use archived WiNN per-seed endpoints on the same verified inputs, not a new fit.
for(seed in sprintf('SIM%02d',1:10)) {
 r<-rows[rows$seed_id==seed & rows$method_id=='raw',];r$method_id<-'winn_fixed_no_qc';r$floored_cells<-NA
 for(metric in c('heldout_qc_cv_mean','truth_log1p_rmse','truth_metabolite_profile_icc_mean','truth_sample_profile_pearson_mean'))
  r[[metric]]<-reference$value[reference$seed_id==seed & reference$method_id=='winn_fixed_no_qc' & reference$metric==metric]
 rows<-rbind(rows,r)
}
write.csv(rows,file.path(out,'baseline_simulation_seed_metrics.csv'),row.names=FALSE)
write.csv(checks,file.path(out,'baseline_simulation_input_checks.csv'),row.names=FALSE)
summary<-do.call(rbind,lapply(split(rows,rows$method_id),function(d)do.call(rbind,lapply(c('heldout_qc_cv_mean','truth_log1p_rmse','truth_metabolite_profile_icc_mean','truth_sample_profile_pearson_mean'),function(metric) data.frame(method_id=d$method_id[1],metric=metric,n=nrow(d),mean=mean(d[[metric]]),se=sd(d[[metric]])/sqrt(nrow(d)))))))
write.csv(summary,file.path(out,'baseline_simulation_summary.csv'),row.names=FALSE)
diffs<-do.call(rbind,lapply(c('mean_batch_all','wehrens_ancova_all'),function(method)do.call(rbind,lapply(c('truth_log1p_rmse','truth_metabolite_profile_icc_mean','truth_sample_profile_pearson_mean'),function(metric) {
 a<-rows[rows$method_id=='winn_fixed_no_qc',];b<-rows[rows$method_id==method,];b<-b[match(a$seed_id,b$seed_id),]
 direction<-if(metric=='truth_log1p_rmse')-1 else 1;d<-direction*(a[[metric]]-b[[metric]])
 data.frame(comparator=method,metric=metric,mean_winn_advantage=mean(d),paired_se=sd(d)/sqrt(length(d)),n_seeds=length(d),winn_better_seeds=sum(d>0),ci_lower=mean(d)-qt(.975,df=9)*sd(d)/sqrt(10),ci_upper=mean(d)+qt(.975,df=9)*sd(d)/sqrt(10))
}))))
write.csv(diffs,file.path(out,'baseline_simulation_paired_differences.csv'),row.names=FALSE)
if(is.null(reuse))capture.output(sessionInfo(),file=file.path(out,'baseline_simulation_sessionInfo.txt'))
cat('Completed exact-input simulation baseline evaluation:',nrow(rows),'seed-method rows\n')
print(summary)
