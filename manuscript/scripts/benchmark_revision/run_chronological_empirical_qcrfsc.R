#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
dataset <- sub('^--dataset=','',args[grepl('^--dataset=',args)])
if(length(dataset)!=1L) stop('Supply --dataset=<key>')
all_args<-commandArgs(FALSE);file_arg<-grep('^--file=',all_args,value=TRUE);script<-if(length(file_arg))sub('^--file=','',file_arg[1]) else all_args[match('-f',all_args)+1L]
root<-dirname(normalizePath(script))
arg_value<-function(key,default){v<-args[startsWith(args,paste0(key,'='))];if(length(v))sub(paste0('^',key,'='),'',v[1]) else default}
input_root<-normalizePath(arg_value('--input-root',root),mustWork=TRUE)
output_root<-arg_value('--output-root',file.path(root,'chronology_results_20261122'))
source(file.path(root,'original_benchmark_helpers.R'))
source(file.path(root,'training_qc_tuning_policy.R'))
source(file.path(root,'training_qc_tuning_engine.R'))
source(file.path(root,'qcrfsc_review_adapter.R'))
if (as.character(packageVersion('statTarget')) != '1.38.0' && !'--allow-version-smoke' %in% args) stop('Production rerun requires frozen statTarget 1.38.0')
suppressPackageStartupMessages(library(dplyr))
p <- readRDS(file.path(input_root,paste0('input_',dataset,'.rds')))
out <- file.path(output_root,dataset);dir.create(out,recursive=TRUE,showWarnings=FALSE)
log_line <- function(...) cat(format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),paste0(...,collapse=''),'\n')
log_line('Begin ',dataset,' statTarget ',packageVersion('statTarget'),' ',nrow(p$x),' features x ',ncol(p$x),' samples')
# Correction and candidate ranking receive only technical metadata and authorized controls.
technical <- tq_technical_metadata(p$meta)
registry <- p$registry
if('--smoke' %in% args) {registry<-registry[1,,drop=FALSE];p$tuning_feature_ids<-head(p$tuning_feature_ids,8);p$folds$rows<-head(p$folds$rows,1)}
x_tune <- p$x[p$tuning_feature_ids,,drop=FALSE]
write.csv(data.frame(feature_id=p$tuning_feature_ids),file.path(out,'tuning_feature_ids.csv'),row.names=FALSE)
write.csv(registry,file.path(out,'candidate_registry.csv'),row.names=FALSE)
.rfsc_inner_audit <- list()
# Candidate-fold jobs are independent and retain the exact frozen per-fit seed.
# Forking changes execution order only; each fit uses its own working directory.
tasks <- unlist(lapply(seq_len(nrow(registry)),function(j)
  lapply(p$folds$rows,function(fold)list(cfg=registry[j,,drop=FALSE],fold=fold))),recursive=FALSE)
cores <- as.integer(arg_value('--cores','1'))
fit_task <- function(task) {
  cfg<-task$cfg;fold<-task$fold
  .rfsc_inner_audit <<- list()
  seed<-tq_candidate_seed(dataset,'QC-RFSC',cfg$candidate_id,fold$fold_id,'candidate')
  set.seed(seed)
  fit_ids<-fold$fit_qc_ids[[1]]
  result<-tryCatch(run_qc_rfsc_review(x_tune,fit_ids,tq_method_metadata(technical,fit_ids),ntree=cfg$ntree,coCV=cfg$coCV,Frule=cfg$Frule),error=function(e)e)
  if(inherits(result,'error')) {
    score<-data.frame(n_validation_qc=length(fold$validation_qc_ids[[1]]),n_features_input=nrow(x_tune),n_features_scored=0,mean_training_qc_profile_correlation=NA_real_,mean_training_qc_cv_fraction=NA_real_,qc_only_score=NA_real_)
    error<-conditionMessage(result)
  } else {score<-training_qc_profile_score(result,fold$validation_qc_ids[[1]]);error<-''}
  list(score=cbind(data.frame(candidate_id=cfg$candidate_id,fold_id=fold$fold_id,seed=seed,error=error),score),audit=.rfsc_inner_audit)
}
results<-parallel::mclapply(tasks,fit_task,mc.cores=cores,mc.preschedule=FALSE,mc.set.seed=FALSE)
stopifnot(all(vapply(results,is.list,logical(1))))
fold_scores<-lapply(results,`[[`,'score')
.rfsc_inner_audit<-unlist(lapply(results,`[[`,'audit'),recursive=FALSE)
write.csv(do.call(rbind,fold_scores),file.path(out,'fold_scores.csv'),row.names=FALSE)
scores <- do.call(rbind,fold_scores)
summary <- do.call(rbind,lapply(split(scores,scores$candidate_id),function(z)data.frame(candidate_id=z$candidate_id[1],mean_qc_only_score=if(all(is.finite(z$qc_only_score)))mean(z$qc_only_score) else NA_real_,completed_folds=sum(is.finite(z$qc_only_score)),total_folds=nrow(z))))
summary <- summary[order(-summary$mean_qc_only_score,summary$candidate_id,na.last=TRUE),]
if(!is.finite(summary$mean_qc_only_score[1])) {write.csv(summary,file.path(out,'selection_summary.csv'),row.names=FALSE); stop('All candidates failed')}
selected <- registry[match(summary$candidate_id[1],registry$candidate_id),,drop=FALSE]
summary$selected <- summary$candidate_id==selected$candidate_id
write.csv(summary,file.path(out,'selection_summary.csv'),row.names=FALSE)
write.csv(selected,file.path(out,'selected_parameters.csv'),row.names=FALSE)
log_line('Selected candidate ',selected$candidate_id,'; full matrix fit')
fit_x <- if('--smoke'%in%args) x_tune else p$x
set.seed(tq_candidate_seed(dataset,'QC-RFSC',selected$candidate_id,'full_training','post_selection'))
corrected <- run_qc_rfsc_review(fit_x,p$training_ids,tq_method_metadata(technical,p$training_ids),ntree=selected$ntree,coCV=selected$coCV,Frule=selected$Frule)
saveRDS(corrected,file.path(out,'corrected_training_qc_abundance.rds'))
write.csv(attr(corrected,'qcrfsc_sample_id_map'),file.path(out,'sample_id_map.csv'),row.names=FALSE)
write.csv(do.call(rbind,.rfsc_inner_audit),file.path(out,'REGfit_training_index_assertions.csv'),row.names=FALSE)
stopifnot(all(vapply(.rfsc_inner_audit,function(z) z$chronological && z$unauthorized==0L,logical(1))))
# Reuse the immutable historical diagnostic matrix; it does not enter selection.
legacy_path<-file.path(root,'results',dataset,'legacy_reproduced_matrix.rds')
legacy<-readRDS(legacy_path)
if('--smoke'%in%args)legacy<-legacy[rownames(fit_x),colnames(fit_x),drop=FALSE]
saveRDS(legacy,file.path(out,'legacy_reproduced_matrix.rds'))
metric <- function(x,ids) {ids<-intersect(ids,colnames(x));if(length(ids)<2)return(NA_real_);mean(apply(x[,ids,drop=FALSE],1,function(z)sd(z)/abs(mean(z))),na.rm=TRUE)*100}
metrics <- do.call(rbind,lapply(c('raw','legacy','corrected'),function(method){x<-switch(method,raw=fit_x,legacy=legacy,corrected=corrected);data.frame(dataset=dataset,method=method,training_cv_percent=metric(x,p$training_ids),heldout_cv_percent=metric(x,c(p$hidden_ids,p$external_ids)),median_intensity=median(x),nonfinite=sum(!is.finite(x)),uncorrected_features=if(is.null(attr(x,'qcrfsc_uncorrected_feature_count')))0 else attr(x,'qcrfsc_uncorrected_feature_count'),uncorrected_samples=if(is.null(attr(x,'qcrfsc_uncorrected_sample_count')))0 else attr(x,'qcrfsc_uncorrected_sample_count'))}))
write.csv(metrics,file.path(out,'revision_diagnostic_metrics.csv'),row.names=FALSE)
writeLines(c('status=complete',paste0('dataset=',dataset),paste0('statTarget=',packageVersion('statTarget')),paste0('mode=',if('--smoke'%in%args)'smoke' else 'full'),paste0('input_sha256=',digest::digest(p$x,algo='sha256')),paste0('output_sha256=',digest::digest(corrected,algo='sha256'))),file.path(out,'COMPLETE.txt'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
log_line('Completed ',dataset)
