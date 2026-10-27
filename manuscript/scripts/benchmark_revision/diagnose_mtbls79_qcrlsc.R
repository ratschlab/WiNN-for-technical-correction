args<-commandArgs(TRUE);opt<-function(k,d){v<-args[startsWith(args,paste0(k,'='))];if(length(v))sub(paste0('^',k,'='),'',v[1]) else d};all_args<-commandArgs(FALSE);file_arg<-grep('^--file=',all_args,value=TRUE);script<-if(length(file_arg))sub('^--file=','',file_arg[1]) else all_args[match('-f',all_args)+1L];here<-dirname(normalizePath(script));input_root<-normalizePath(opt('--input-root',here));out<-opt('--output-root',file.path(here,'qcrlsc_diagnostic'));dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(dplyr))
source(file.path(here,'original_benchmark_helpers.R'));source(file.path(here,'public_benchmark_method_helpers.R'))
p<-readRDS(file.path(input_root,'input_mtbls79.rds'));md<-p$meta
warnings<-character()
x<-withCallingHandlers(run_qc_rlsc_id_safe(p$x,p$training_ids,md,span=.7,degree=1,shift_batches=TRUE,min_qcs_per_batch=3),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
saveRDS(x,file.path(out,'mtbls79_qcrlsc_diagnostic_matrix.rds'))
cv<-function(m)apply(m[,p$hidden_ids,drop=FALSE],1,function(z)sd(z)/mean(z)*100)
a<-data.frame(feature_id=rownames(x),raw_heldout_cv_percent=cv(p$x),corrected_heldout_cv_percent=cv(x))
a<-a[order(-a$corrected_heldout_cv_percent),];write.csv(a,file.path(out,'mtbls79_qcrlsc_feature_diagnostics.csv'),row.names=FALSE)
writeLines(unique(warnings),file.path(out,'mtbls79_qcrlsc_warnings.txt'))
trace<-list();coverage<-list()
for(batch in unique(md$batch)){
 idx<-which(md$batch==batch);q<-which(md$sample_id[idx]%in%p$training_ids);ord<-seq_along(idx)
 coverage[[as.character(batch)]]<-data.frame(batch=batch,n_injections=length(idx),training_qcs=length(q),first_qc_position=min(q),last_qc_position=max(q),left_extrapolation_injections=sum(ord<min(q)),right_extrapolation_injections=sum(ord>max(q)),span=.7,degree=1)
 for(feature in head(a$feature_id,5)){
  yy<-log1p(p$x[feature,idx]);fit<-suppressWarnings(loess(yy[q]~q,span=.7,degree=1,control=loess.control(surface='direct')));pred<-suppressWarnings(predict(fit,data.frame(q=ord)))
  trace[[paste(batch,feature)]]<-data.frame(batch=batch,feature_id=feature,within_batch_order=ord,sample_id=md$sample_id[idx],training_qc=md$sample_id[idx]%in%p$training_ids,heldout_qc=md$sample_id[idx]%in%p$hidden_ids,raw_log1p=yy,fitted_log1p=pred,within_batch_corrected_log1p=yy-pred+mean(yy[q]))
 }
}
write.csv(do.call(rbind,trace),file.path(out,'mtbls79_qcrlsc_failing_profile_traces.csv'),row.names=FALSE)
write.csv(do.call(rbind,coverage),file.path(out,'mtbls79_qcrlsc_coverage.csv'),row.names=FALSE)
write.csv(data.frame(method=c('raw','qcrlsc_selected'),heldout_mean_cv_percent=c(mean(cv(p$x)),mean(cv(x))),warning_count=c(0,length(warnings)),qcrlsc_version=as.character(packageVersion('qcrlscR'))),file.path(out,'mtbls79_qcrlsc_summary.csv'),row.names=FALSE)
