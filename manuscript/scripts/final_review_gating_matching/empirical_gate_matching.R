#!/usr/bin/env Rscript
# Production 0.1.6 snapshot; same prepared inputs and gate settings as manuscript.
a<-commandArgs(FALSE);root<-dirname(dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',a,value=TRUE))))))
args<-commandArgs(TRUE);input<-if(length(args))args[1]else'/tmp/winn-final-review-inputs'
source(file.path(root,'scripts/targeted_validation/source_snapshot/winn/R/winn.R'))
source(file.path(root,'scripts/targeted_validation/source_snapshot/winn/R/winn-ablation.R'))
source(file.path(root,'scripts/nearest_baseline/derivative_biology_helpers.R'))
out<-file.path(root,'generated/final_review_gating_matching');dir.create(out,recursive=TRUE,showWarnings=FALSE)
results<-differences<-protocol<-manifests<-list()
for(dataset in c('batchcorr_set1','clinical_fiams')) {
 p<-readRDS(file.path(input,paste0('input_',dataset,'.rds')));md<-p$meta;study<-which(md$is_study)
 stopifnot(identical(colnames(p$x),as.character(md$sample_id)))
 matrices<-list(raw=p$x)
 for(policy in c('selective','forced')) {
  cache<-file.path(input,paste0(dataset,'_',policy,'.rds'))
  if(file.exists(cache))v<-readRDS(cache) else {
   cat('Fitting',dataset,policy,'\n')
   v<-winn_ablation(p$x,batch=md$batch,run_order=md$run_order,parameters='fixed',
    use_outlier_shrinkage=TRUE,drift_gate=if(policy=='selective')'selective'else'all',batch_gate='selective',
    pqn_mode='shrink',fdr_threshold=.05,test='Ljung-Box',lag=NULL,ljung_box_fitdf=0L,
    spline_method='conservative',return_intermediates=FALSE,return_diagnostics=FALSE)$data
   saveRDS(v,cache)
  }
  stopifnot(identical(dimnames(v),dimnames(p$x)),all(is.finite(v)),all(v>=0))
  matrices[[policy]]<-v
 }
 field<-if(dataset=='batchcorr_set1')'accession_id'else'sample'
 smd<-md[study,,drop=FALSE];id<-as.character(smd[[field]]);units<-unique(id)
 queries<-lapply(matrices,function(m)cross_batch_matching(log1p(m[,study,drop=FALSE]),smd,field)$query)
 stopifnot(identical(queries$raw$query_index,queries$selective$query_index),identical(queries$raw$query_index,queries$forced$query_index))
 weights<-as.numeric(table(factor(queries$raw$identity,levels=units)))
 draws<-cluster_bootstrap_draws(units,NULL,2000,2026072701L+if(dataset=='batchcorr_set1')199L else 1594L)
 counts<-t(vapply(draws,function(d)tabulate(match(d,units),nbins=length(units)),numeric(length(units))))
 boot<-points<-list()
 for(policy in names(matrices)) {
  q<-queries[[policy]]; numer<-as.numeric(tapply(q$top1,factor(q$identity,levels=units),sum)); numer[is.na(numer)]<-0
  boot[[policy]]<-as.numeric(counts%*%numer)/as.numeric(counts%*%weights);points[[policy]]<-mean(q$top1)
  results[[length(results)+1]]<-data.frame(dataset=dataset,policy=policy,accuracy=points[[policy]],ci_lower=quantile(boot[[policy]],.025),ci_upper=quantile(boot[[policy]],.975),queries=nrow(q),identities=sum(weights>0))
  sim<-cor(log1p(matrices[[policy]][,study,drop=FALSE]));bat<-as.character(smd$batch)
  audit<-do.call(rbind,lapply(q$query_index,function(i){can<-which(bat!=bat[i]&is.finite(sim[i,])); nmatch<-sum(id[can]==id[i]);data.frame(candidates=length(can),same_identity_candidates=nmatch,chance=nmatch/length(can),top_tie=sum(sim[i,can]==max(sim[i,can]))>1)}))
  protocol[[length(protocol)+1]]<-data.frame(dataset=dataset,policy=policy,study_profiles=length(study),eligible_queries=nrow(q),eligible_identities=sum(weights>0),all_identities=length(units),candidate_min=min(audit$candidates),candidate_max=max(audit$candidates),matching_candidates_min=min(audit$same_identity_candidates),matching_candidates_max=max(audit$same_identity_candidates),chance_accuracy=mean(audit$chance),exact_top_ties=sum(audit$top_tie))
  manifests[[length(manifests)+1]]<-data.frame(dataset=dataset,policy=policy,matrix_sha256=digest::digest(matrices[[policy]],algo='sha256'),features=nrow(matrices[[policy]]),profiles=ncol(matrices[[policy]]))
 }
 archived<-read.csv(file.path(root,'generated/review_validation/frozen_identity_retrieval_estimates.csv'))
 for(policy in c('raw','selective')) {
  target<-archived$point_estimate[archived$dataset==dataset & archived$method_id==(if(policy=='raw')'raw'else'winn_fixed_no_qc') & archived$metric=='cross_batch_top1_accuracy']
  stopifnot(length(target)==1,abs(points[[policy]]-target)<1e-12)
  archived_row<-archived[archived$dataset==dataset & archived$method_id==(if(policy=='raw')'raw'else'winn_fixed_no_qc') & archived$metric=='cross_batch_top1_accuracy',]
  stopifnot(max(abs(quantile(boot[[policy]],c(.025,.975))-unlist(archived_row[c('ci_lower','ci_upper')])))<1e-12)
 }
 for(comparator in c('raw','forced')) {
  d<-boot$selective-boot[[comparator]]
  differences[[length(differences)+1]]<-data.frame(dataset=dataset,contrast=paste0('selective_minus_',comparator),difference=points$selective-points[[comparator]],ci_lower=quantile(d,.025),ci_upper=quantile(d,.975),bootstrap_replicates=2000)
 }
 cat(dataset,'complete\n')
}
for(spec in list(c('results','empirical_gate_matching_estimates.csv'),c('differences','empirical_gate_matching_differences.csv'),c('protocol','matching_protocol_audit.csv'),c('manifests','empirical_gate_matching_matrix_manifest.csv')))write.csv(do.call(rbind,get(spec[1])),file.path(out,spec[2]),row.names=FALSE)
print(do.call(rbind,results));print(do.call(rbind,differences))
