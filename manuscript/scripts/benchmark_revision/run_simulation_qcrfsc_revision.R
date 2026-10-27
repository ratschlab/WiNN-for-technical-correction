#!/usr/bin/env Rscript
args<-commandArgs(TRUE);opt<-function(k,d=NULL){v<-args[startsWith(args,paste0(k,'='))];if(length(v))sub(paste0('^',k,'='),'',v[1]) else d}
seed_id<-opt('--seed-id');if(is.null(seed_id)||!seed_id%in%sprintf('SIM%02d',1:10))stop('Supply --seed-id=SIM01 through SIM10')
all_args<-commandArgs(FALSE);f<-grep('^--file=',all_args,value=TRUE);script<-if(length(f))sub('^--file=','',f[1])else all_args[match('-f',all_args)+1L];here<-dirname(normalizePath(script))
input_root<-normalizePath(opt('--input-root',here));out<-file.path(opt('--output-root',file.path(here,'simulation_results')),seed_id);dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(here,'original_benchmark_helpers.R'));source(file.path(here,'qcrfsc_review_adapter.R'));source(file.path(here,'frozen_winn_ablation_helpers.R'))
suppressPackageStartupMessages(library(dplyr))
if(as.character(packageVersion('statTarget'))!='1.38.0')stop('Frozen statTarget1.38.0 required')
p<-readRDS(file.path(input_root,paste0('input_',seed_id,'.rds')))
hashes<-read.csv(file.path(here,'simulation_bundle_hashes.csv'))
checks<-do.call(rbind,lapply(c('raw_intensity','clean_ground_truth'),function(obj){x<-if(obj=='raw_intensity')p$x else p$truth;expected<-hashes[hashes$seed_id==seed_id&hashes$object==obj,];observed<-digest::digest(x,algo='sha256');data.frame(seed_id=seed_id,object=obj,observed_sha256=observed,expected_sha256=expected$reference_object_sha256,passed=identical(observed,expected$reference_object_sha256))}))
write.csv(checks,file.path(out,'input_identity_checks.csv'),row.names=FALSE);if(!all(checks$passed))stop('Frozen simulation identity mismatch')
md<-p$meta;technical<-md[,c('sample_id','batch','run_order','within_batch_order')]
rng<-2026072610L+as.integer(sub('SIM','',seed_id))*100L+4L
set.seed(rng);legacy<-run_qc_rfsc_with_controls(p$x,p$training_ids,technical,ntree=500,coCV=30,Frule=.8)
.audit<-list()
invisible(capture.output(trace('REGfit',where=asNamespace('statTarget'),print=FALSE,tracer=quote({
 cn0<-colnames(x);used0<-cn0[grepl('QC|qc|Qc',cn0)];expected0<-intersect(cn0,getOption('qcrfsc_authorized_names'));stopifnot(setequal(used0,expected0));a<-get('.audit',envir=.GlobalEnv);a[[length(a)+1L]]<-data.frame(call=length(a)+1L,available=length(cn0),used=length(used0),authorized=length(expected0),unauthorized=sum(!used0%in%expected0));assign('.audit',a,envir=.GlobalEnv)
}))))
set.seed(rng);corrected<-run_qc_rfsc_review(p$x,p$training_ids,technical,ntree=500,coCV=30,Frule=.8)
untrace('REGfit',where=asNamespace('statTarget'))
write.csv(do.call(rbind,.audit),file.path(out,'REGfit_training_index_assertions.csv'),row.names=FALSE)
saveRDS(corrected,file.path(out,'corrected_training_qc_abundance.rds'));saveRDS(legacy,file.path(out,'legacy_reproduced_matrix.rds'))
study<-md$sample_id[!md$is_qc];truth<-log1p(p$truth[,study,drop=FALSE])
results<-do.call(rbind,lapply(c('raw','legacy','corrected'),function(method){x<-switch(method,raw=p$x,legacy=legacy,corrected=corrected);candidate<-log1p(x[,study,drop=FALSE]);feature_icc<-vapply(seq_len(nrow(candidate)),function(j)icc_a1_ablation(cbind(truth[j,],candidate[j,])),numeric(1));sample_cor<-vapply(seq_len(ncol(candidate)),function(j)safe_cor_ablation(truth[,j],candidate[,j]),numeric(1));q<-x[,p$hidden_ids,drop=FALSE];heldcv<-mean(apply(q,1,function(z)sd(z)/mean(z)),na.rm=TRUE)*100;fallback<-attr(x,'qcrfsc_uncorrected_feature_count');if(is.null(fallback))fallback<-0;data.frame(seed_id=seed_id,method=method,truth_log1p_rmse=sqrt(mean((candidate-truth)^2)),truth_metabolite_profile_icc_mean=mean(feature_icc,na.rm=TRUE),truth_sample_profile_pearson_mean=mean(sample_cor,na.rm=TRUE),heldout_qc_cv_mean=heldcv,uncorrected_features=fallback,rng=rng,ntree=500,coCV=30,Frule=.8)}))
write.csv(results,file.path(out,'simulation_truth_metrics.csv'),row.names=FALSE)
writeLines(c('status=complete',paste0('seed_id=',seed_id),'statTarget=1.38.0','exact_frozen_input_hashes_verified=TRUE'),file.path(out,'COMPLETE.txt'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
