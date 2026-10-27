#!/usr/bin/env Rscript
# Aggregate public phenotype placement; export no sample or participant identifiers.
args <- commandArgs(trailingOnly=TRUE)
getarg <- function(key, default) {x<-grep(paste0('^--',key,'='),args,value=TRUE);if(length(x))sub(paste0('^--',key,'='),'',x[1]) else default}
input_root <- getarg('input-root','.')
output_root <- getarg('output-root','.')
dir.create(output_root,recursive=TRUE,showWarnings=FALSE)
groups <- covariates <- list()
for(dataset in c('mtbls79','sacurine','waveica')) {
 p<-readRDS(file.path(input_root,paste0('input_',dataset,'.rds')));d<-p$meta
 stopifnot(identical(as.character(d$sample_id),colnames(p$x)))
 for(b in unique(d$batch)) {
  v<-d[d$batch==b,];v<-v[order(v$run_order),]
  v$quartile<-pmin(4L,ceiling(seq_len(nrow(v))*4/nrow(v)))
  v<-v[v$is_study,]
  v$group<-as.character(v[[switch(dataset,mtbls79='class',sacurine='gender',waveica='biological_group')]])
  g<-as.data.frame(table(quartile=factor(v$quartile,levels=1:4),group=v$group),stringsAsFactors=FALSE)
  groups[[length(groups)+1L]]<-data.frame(dataset=dataset,batch=b,acquisition_quartile=g$quartile,group=g$group,study_units=g$Freq)
  if(dataset=='sacurine')for(q in 1:4)for(covariate in c('age','bmi','osmolality')) {
   x<-v[[covariate]][v$quartile==q];x<-x[is.finite(x)]
   covariates[[length(covariates)+1L]]<-data.frame(dataset=dataset,batch=b,acquisition_quartile=q,covariate=covariate,n=length(x),median=median(x),q25=unname(quantile(x,.25)),q75=unname(quantile(x,.75)))
  }
 }
}
write.csv(do.call(rbind,groups),file.path(output_root,'design_group_acquisition_quartiles.csv'),row.names=FALSE)
write.csv(do.call(rbind,covariates),file.path(output_root,'design_covariate_acquisition_quartiles.csv'),row.names=FALSE)
cat('Public acquisition-quartile aggregates completed.\n')
