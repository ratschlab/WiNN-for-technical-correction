# Official Wehrens ANCOVA and mean-only batch baselines, evaluated descriptively.
# All-injection fits receive the same intensity/run-order/batch information as fixed WiNN.
args<-commandArgs(TRUE)
value<-function(key,default=NULL){v<-args[startsWith(args,paste0(key,'='))];if(!length(v))return(default);sub(paste0('^',key,'='),'',v[1])}
input_root<-value('--input-root');if(is.null(input_root))stop('Supply --input-root=/path/containing/input_<dataset>.rds')
input_root<-normalizePath(input_root,mustWork=TRUE)
out<-value('--output-root',file.path(getwd(),'nearest_baseline_results'));dir.create(out,recursive=TRUE,showWarnings=FALSE)
out<-normalizePath(out)
script_arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)[1]
script_dir<-dirname(normalizePath(sub('^--file=','',script_arg)))
source(file.path(script_dir,'doBC.R'))
source(file.path(script_dir,'weighted_pc_r2.R'))
datasets<-strsplit(value('--datasets','mtbls79,clinical_fiams,batchcorr_set1,sacurine,waveica'),',',fixed=TRUE)[[1]]
metrics<-list()
for(dataset in datasets) {
 p<-readRDS(file.path(input_root,paste0('input_',dataset,'.rds')));md<-p$meta
 study<-which(md$is_study);all_idx<-seq_len(ncol(p$x));z<-log1p(p$x)
 methods<-list(raw=p$x)
 for(spec in c('mean_batch_all','wehrens_ancova_all','wehrens_ancova_study')) {
  idx<-if(spec=='wehrens_ancova_study')study else all_idx
  order<-if(spec=='mean_batch_all')NULL else md$run_order
  zz<-t(apply(z,1,doBC,ref.idx=idx,batch.idx=md$batch,seq.idx=order,method='lm',minBsamp=if(spec=='mean_batch_all')2 else 4))
  dimnames(zz)<-dimnames(p$x)
  if(any(!is.finite(zz)))stop('Official baseline returned nonfinite values: ',dataset,'/',spec)
  x<-pmax(expm1(zz),0);dimnames(x)<-dimnames(p$x)
  attr(x,'floored_cells')<-sum(zz<0)
  methods[[spec]]<-x
  dir.create(file.path(out,'nearest_baseline_matrices',dataset),recursive=TRUE,showWarnings=FALSE)
  saveRDS(x,file.path(out,'nearest_baseline_matrices',dataset,paste0(spec,'.rds')))
 }
 for(method in names(methods)) {
  x<-methods[[method]];held<-intersect(c(p$hidden_ids,p$external_ids),colnames(x));q<-x[,held,drop=FALSE]
  cv<-mean(apply(q,1,function(z)sd(z)/abs(mean(z))),na.rm=TRUE)*100
  r2<-weighted_pc_r2_explicit(x[,study,drop=FALSE],md$batch[study],target_type='categorical')
  group<-if(dataset=='waveica')md$biological_group else if(dataset=='mtbls79')md$class else NULL
  group_r2<-if(!is.null(group))weighted_pc_r2_explicit(x[,study,drop=FALSE],group[study],target_type='categorical') else NA_real_
  floors<-attr(x,'floored_cells');if(is.null(floors))floors<-0
  metrics[[paste(dataset,method)]]<-data.frame(dataset=dataset,method=method,heldout_cv_percent=cv,study_batch_weighted_pc_r2=r2,study_group_weighted_pc_r2=group_r2,floored_cells=floors,information=if(method=='wehrens_ancova_study')'intensity, order, batch, reference/study identity' else 'intensity, order, batch',output_scale='nonnegative intensity; log1p regression; expm1 and floor at zero')
 }
 write.csv(do.call(rbind,metrics),file.path(out,'nearest_baseline_metrics.csv'),row.names=FALSE)
 cat('Completed official Wehrens and mean-only baseline:',dataset,'\n')
}
