#!/usr/bin/env Rscript
# Existing-table point endpoints, using the same frozen feature model.
full<-commandArgs(FALSE);here<-dirname(normalizePath(sub('^--file=','',grep('^--file=',full,value=TRUE)[1])))
source(file.path(here,'frozen_derivative_biology_helpers.R'))
p<-readRDS(file.path(here,'input_waveica.rds'))
x<-readRDS(file.path(here,'chronology_results_20261122/waveica/corrected_training_qc_abundance.rds'))
stopifnot(identical(dimnames(x),dimnames(p$x)))
md<-p$meta[p$meta$is_study,,drop=FALSE]
result<-do.call(rbind,lapply(c('raw','qc_rfsc'),function(method){
  value<-if(method=='raw')p$x else x
  fit<-two_group_model(log1p(value[,md$sample_id,drop=FALSE]),md$biological_group,'group_0','group_1')
  data.frame(dataset_key='waveica',method_id=method,biological_associated_features=sum(fit$significant),standardized_effect_magnitude=median(abs(fit$standardized_effect),na.rm=TRUE),within_group_variance=median(fit$within_group_variance,na.rm=TRUE),features=nrow(fit),study_samples=nrow(md),model='ordinary pooled-variance two-group model; BH threshold 0.05; no permutation null rerun')
}))
write.csv(result,file.path(here,'chronology_results_20261122/waveica/extra_biology_points.csv'),row.names=FALSE)
source(file.path(here,'frozen_winn_ablation_helpers.R'))
p<-readRDS(file.path(here,'input_mtbls79.rds'))
x<-readRDS(file.path(here,'chronology_results_20261122/mtbls79/corrected_training_qc_abundance.rds'))
result<-do.call(rbind,lapply(c('raw','qc_rfsc'),function(method){
  value<-if(method=='raw')p$x else x
  icc<-mtbls79_feature_cross_batch_icc(value,p$meta)
  data.frame(dataset_key='mtbls79',method_id=method,metabolite_icc_a1_median=median(icc,na.rm=TRUE),features=length(icc),endpoint='median across features of median pairwise-batch ICC(A,1), using shared extracts; point estimate')
}))
stopifnot(abs(result$metabolite_icc_a1_median[result$method_id=='raw']-0.548997409303439)<1e-12)
write.csv(result,file.path(here,'chronology_results_20261122/mtbls79/extra_biology_points.csv'),row.names=FALSE)
