#!/usr/bin/env Rscript
# Evaluate existing same-information baseline fits, and align retrieval intervals
# with their query-weighted point estimand. No correction is refitted here.
args <- commandArgs(TRUE)
opt <- function(k,d=NULL){z<-args[startsWith(args,paste0(k,'='))];if(length(z))sub(paste0('^',k,'='),'',z[1])else d}
here<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
root<-dirname(dirname(here));input<-normalizePath(opt('--input-root'),mustWork=TRUE)
matrices_root<-normalizePath(opt('--matrix-root',file.path(input,'nearest_baseline_matrices')),mustWork=TRUE)
derivative<-normalizePath(opt('--derivative-root'),mustWork=TRUE)
out<-opt('--output-root',file.path(root,'generated/review_validation'));dir.create(out,recursive=TRUE,showWarnings=FALSE)
datasets<-strsplit(opt('--datasets','batchcorr_set1,clinical_fiams,mtbls79,waveica,sacurine'),',',fixed=TRUE)[[1]]
if(!length(datasets)||anyDuplicated(datasets)||!all(datasets%in%c('batchcorr_set1','clinical_fiams','mtbls79','waveica','sacurine')))stop('Invalid --datasets selection')
source(file.path(here,'derivative_biology_helpers.R'))
read_input_csv<-function(path,...)read.csv(if(file.exists(path))path else paste0(path,'.gz'),...)
labels<-c(raw='Raw',combat='ComBat',qc_rlsc='QC-RLSC',tiger='TIGER',serrf='SERRF',
 winn_auto_qc='WiNN (QC-tuned)',winn_auto_batch_qc='WiNN (inferred segments)',winn_fixed_no_qc='WiNN',
 mean_batch_all='Mean-only batch',wehrens_ancova_all='Study-sample ANCOVA')
old_path<-file.path(root,'generated/review_validation/frozen_identity_retrieval_estimates_before_query_weighting.csv')
old_est<-read_input_csv(if(file.exists(old_path))old_path else file.path(root,'generated/review_validation/frozen_identity_retrieval_estimates.csv'))
for(filename in c('frozen_identity_retrieval_estimates.csv','frozen_identity_retrieval_paired_differences.csv')) {
 archive<-file.path(out,sub('.csv$','_before_query_weighting.csv',filename));if(!file.exists(archive))file.copy(file.path(root,'generated/review_validation',filename),archive)
}
checks<-list();check<-function(label,test,value=''){checks[[length(checks)+1L]]<<-data.frame(check=label,passed=test,detail=value);if(!test)stop(label,': ',value)}
input_manifest<-list();estimates<-list();differences<-list();bootstrap_rows<-list();protocol<-list();effect_rows<-list();lobo_rows<-list()
for(dataset in datasets) {
 p<-readRDS(file.path(input,paste0('input_',dataset,'.rds')));md<-p$meta;check(paste(dataset,'metadata alignment'),identical(colnames(p$x),as.character(md$sample_id)))
 values<-list(raw=p$x)
 for(method in c('mean_batch_all','wehrens_ancova_all')) {
  path<-file.path(matrices_root,dataset,paste0(method,'.rds'));x<-readRDS(path)
  check(paste(dataset,method,'dimensions and identifiers'),identical(dimnames(x),dimnames(p$x)) && all(is.finite(x)))
  values[[method]]<-x
  input_manifest[[length(input_manifest)+1L]]<-data.frame(dataset=dataset,object=method,sha256=digest::digest(x,algo='sha256'),features=nrow(x),samples=ncol(x),origin='official doBC all-injection baseline')
 }
 input_manifest[[length(input_manifest)+1L]]<-data.frame(dataset=dataset,object='raw_input',sha256=digest::digest(p$x,algo='sha256'),features=nrow(p$x),samples=ncol(p$x),origin=if(dataset=='clinical_fiams')'clinical pair reconstruction'else'prepared benchmark input')
 study<-which(md$is_study);study_md<-md[study,,drop=FALSE]
 if(dataset %in% c('batchcorr_set1','clinical_fiams')) {
  field<-if(dataset=='batchcorr_set1')'accession_id'else'sample'
  identities<-as.character(study_md[[field]]);units<-unique(identities)
  matches<-lapply(values,function(x)cross_batch_matching(log1p(x[,study,drop=FALSE]),study_md,field))
  query_counts<-table(matches$raw$query$identity)
  unit_weights<-as.numeric(query_counts[match(units,names(query_counts))]);unit_weights[is.na(unit_weights)]<-0
  check(paste(dataset,'raw queries distinct'),!anyDuplicated(matches$raw$query$query_index))
  identity_key<-if(dataset=='clinical_fiams')hash_private_identity(units)else units
  frozen_path<-file.path(derivative,dataset,'matching_metrics_by_identity.csv')
  if(file.exists(frozen_path)||file.exists(paste0(frozen_path,'.gz'))) frozen<-read_input_csv(frozen_path,stringsAsFactors=FALSE) else {
   # Clinical identity tables were intentionally omitted from the archived public
   # release. Its complete bootstrap series is a linear transformation of the
   # identity summaries. Recover those summaries using the exact archived RNG
   # seed and unit order, then validate every series and the direct raw result.
   check('Clinical-only matching reconstruction',dataset=='clinical_fiams')
   legacy<-read_input_csv(file.path(derivative,dataset,'bootstrap_replicates.csv'))
   legacy_draws<-cluster_bootstrap_draws(units,NULL,2000,2026072701L+1594L)
   legacy_counts<-t(vapply(legacy_draws,function(draw)tabulate(match(draw,units),nbins=length(units)),numeric(length(units))))
   eligible<-unit_weights>0
   design<-legacy_counts[,eligible,drop=FALSE]/rowSums(legacy_counts[,eligible,drop=FALSE])
   decomp<-qr(design)
   check('Clinical archived bootstrap system has full column rank',decomp$rank==sum(eligible),decomp$rank)
   frozen<-do.call(rbind,lapply(setdiff(unique(legacy$method_id),'qc_rfsc'),function(method) {
    recovered<-sapply(c('cross_batch_top1_accuracy','cross_batch_mean_reciprocal_rank'),function(metric) {
     selected<-legacy[legacy$method_id==method & legacy$metric==metric,]
     response<-selected$value[match(1:2000,selected$bootstrap)]
     z<-qr.coef(decomp,response)
     residual<-max(abs(design %*% z-response))
     check(paste('Clinical',method,metric,'bootstrap reconstruction'),residual<1e-11,format(residual,digits=15))
     check(paste('Clinical',method,metric,'recovered metric bounds'),all(z >= -1e-10 & z <= 1+1e-10))
     if(metric=='cross_batch_top1_accuracy')z<-round(z*unit_weights[eligible])/unit_weights[eligible]
     pmin(1,pmax(0,z))
    })
    data.frame(identity=identity_key[eligible],top1=recovered[,1],reciprocal_rank=recovered[,2],method_id=method,method=labels[[method]])
   }))
  }
  frozen<-frozen[frozen$method_id!='qc_rfsc',]
  frozen_methods<-unique(frozen$method_id)
  identity_values<-lapply(frozen_methods,function(method){v<-frozen[frozen$method_id==method,];i<-match(identity_key,v$identity);z<-v[i,c('top1','reciprocal_rank')];z[unit_weights==0,]<-0;check(paste(dataset,method,'all eligible identity summaries'),all(is.finite(as.matrix(z[unit_weights>0,]))));as.matrix(z)})
  names(identity_values)<-frozen_methods
  for(method in names(matches)) {
   count<-table(matches[[method]]$query$identity)
   check(paste(dataset,method,'common query set'),identical(matches[[method]]$query$query_index,matches$raw$query$query_index))
   v<-matches[[method]]$identity;i<-match(units,v$identity);z<-as.matrix(v[i,c('top1','reciprocal_rank')]);z[unit_weights==0,]<-0
   if(method=='raw')check(paste(dataset,'recomputed raw retrieval matches archived summaries'),max(abs(z-identity_values$raw))<1e-10)else identity_values[[method]]<-z
  }
  for(method in frozen_methods)for(j in 1:2) {
   metric<-c('cross_batch_top1_accuracy','cross_batch_mean_reciprocal_rank')[j]
   original<-old_est$point_estimate[old_est$dataset==dataset & old_est$method_id==method & old_est$metric==metric]
   actual<-sum(unit_weights*identity_values[[method]][,j])/sum(unit_weights)
   check(paste(dataset,method,metric,'reproduces query-weighted point'),length(original)==1 && abs(actual-original)<1e-10,actual)
  }
  draws<-cluster_bootstrap_draws(units,NULL,2000,2026072701L+if(dataset=='batchcorr_set1')199L else 1594L)
  counts<-t(vapply(draws,function(draw)tabulate(match(draw,units),nbins=length(units)),numeric(length(units))))
  weighted_counts<-sweep(counts,2,unit_weights,'*');denominator<-rowSums(weighted_counts)
  for(j in 1:2) {
   metric<-c('cross_batch_top1_accuracy','cross_batch_mean_reciprocal_rank')[j]
   bs<-sapply(identity_values,function(z)as.vector(weighted_counts %*% z[,j])/denominator)
   point<-vapply(identity_values,function(z)sum(unit_weights*z[,j])/sum(unit_weights),numeric(1))
   for(method in names(identity_values)) {
    v<-bs[,method];ci<-quantile(v,c(.025,.975));estimates[[length(estimates)+1L]]<-data.frame(method_id=method,method=labels[[method]],metric=metric,direction='higher',point_estimate=point[[method]],bootstrap_mean=mean(v),ci_lower=ci[1],ci_upper=ci[2],valid_bootstrap=length(v),dataset=dataset,uncertainty_scope='paired 2000-draw identity-cluster bootstrap retaining query-count weights; fitted correction and candidate set fixed')
    bootstrap_rows[[length(bootstrap_rows)+1L]]<-data.frame(dataset=dataset,bootstrap=seq_along(v),method_id=method,metric=metric,value=v)
    if(method!='winn_fixed_no_qc') {
     delta<-bs[,'winn_fixed_no_qc']-v;ci_d<-quantile(delta,c(.025,.975));differences[[length(differences)+1L]]<-data.frame(metric=metric,direction='higher',winn_method_id='winn_fixed_no_qc',winn_method='WiNN',comparator_method_id=method,comparator_method=labels[[method]],point_direction_normalized_difference=point[['winn_fixed_no_qc']]-point[[method]],bootstrap_mean_difference=mean(delta),ci_lower=ci_d[1],ci_upper=ci_d[2],probability_winn_advantage=mean(delta>1e-12),valid_bootstrap=length(delta),difference_definition='WiNN minus comparator; positive favors WiNN',dataset=dataset,uncertainty_scope='paired identity resampling with query-count weights; fitted correction and candidate set fixed')
    }
   }
  }
  # All historical methods used the same frozen common feature/sample panel.
  # Their identity summaries reproduce the original query-weighted points.
  for(method in setdiff(names(identity_values),names(matches))) {
   protocol[[length(protocol)+1L]]<-data.frame(dataset=dataset,method_id=method,study_injections=nrow(study_md),eligible_queries=sum(unit_weights),eligible_identities=sum(unit_weights>0),all_identities=length(units),candidate_min=NA,candidate_median=NA,candidate_max=NA,exact_top_ties=NA,similarity='Pearson across features on log1p intensity',candidates='historical frozen common study panel; same other-batch eligibility',tie_rule='stable decreasing similarity sort; first in original input-column order wins',bootstrap_replicates=2000)
  }
  # Candidate-set geometry and exact top-score ties for newly evaluated matrices.
  for(method in names(matches)) {
   v<-log1p(values[[method]][,study,drop=FALSE]);sim<-cor(v);diag(sim)<-NA
   qi<-matches[[method]]$query$query_index
   candidate_n<-vapply(qi,function(i)sum(study_md$batch!=study_md$batch[i] & is.finite(sim[i,])),integer(1))
   ties<-vapply(qi,function(i){eligible<-which(study_md$batch!=study_md$batch[i] & is.finite(sim[i,]));sum(sim[i,eligible]==max(sim[i,eligible]))>1},logical(1))
   protocol[[length(protocol)+1L]]<-data.frame(dataset=dataset,method_id=method,study_injections=nrow(study_md),eligible_queries=length(qi),eligible_identities=length(query_counts),all_identities=length(units),candidate_min=min(candidate_n),candidate_median=median(candidate_n),candidate_max=max(candidate_n),exact_top_ties=sum(ties),similarity='Pearson across features on log1p intensity',candidates='all finite-similarity study profiles in other batches; query requires same identity among candidates',tie_rule='stable decreasing similarity sort; first in original input-column order wins',bootstrap_replicates=2000)
  }
 } else {
  z<-lapply(values,function(x)log1p(x[,study,drop=FALSE]));md<-study_md
  if(dataset=='mtbls79'){
   units<-unique(md$sample_rep);group<-vapply(units,function(u)unique(md$class[md$sample_rep==u])[1],character(1));collapsed<-lapply(z,function(x){v<-vapply(units,function(u)rowMeans(x[,md$sample_rep==u,drop=FALSE]),numeric(nrow(x)));rownames(v)<-rownames(x);v})
   fits<-lapply(collapsed,function(x)list(cow_vs_sheep=two_group_model(x,group,'S','C')))
  } else if(dataset=='waveica')fits<-lapply(z,function(x)list(biological_group=two_group_model(x,md$biological_group,'group_0','group_1')))
  else {
   z<-lapply(values,function(x)log1p(sweep(x[,study,drop=FALSE],2,md$osmolality,'/')))
   md$age_z<-as.numeric(scale(md$age));md$bmi_z<-as.numeric(scale(md$bmi));md$gender<-factor(md$gender);md$sampling_factor<-factor(md$sampling)
   design<-model.matrix(~age_z+bmi_z+gender+sampling_factor,md);coefs<-c(age='age_z',bmi='bmi_z',gender=grep('^gender',colnames(design),value=TRUE)[1])
   fits<-lapply(z,function(x)lapply(coefs,function(coef){t<-ols_feature_model(x,design,match(coef,colnames(design)))$table;t$effect<-t$coefficient;t}))
  }
  frozen_effects<-read_input_csv(file.path(derivative,dataset,if(dataset=='sacurine')'sacurine_covariate_associations.csv'else'association_features.csv'))
  if(dataset=='sacurine')frozen_effects$effect<-frozen_effects$coefficient else frozen_effects$endpoint<-names(fits$raw)[1]
  fits$winn_fixed_no_qc<-list()
  for(endpoint in names(fits$raw)) {
   old_raw<-frozen_effects[frozen_effects$method_id=='raw' & frozen_effects$endpoint==endpoint,]
   raw<-fits$raw[[endpoint]]
   check(paste(dataset,endpoint,'raw feature effects match frozen export'),max(abs(raw$effect-old_raw$effect[match(raw$feature_id,old_raw$feature_id)]))<1e-10)
   frozen_fit<-frozen_effects[frozen_effects$method_id=='winn_fixed_no_qc' & frozen_effects$endpoint==endpoint,]
   fits$winn_fixed_no_qc[[endpoint]]<-frozen_fit[match(raw$feature_id,frozen_fit$feature_id),]
  }
  for(method in names(fits))for(endpoint in names(fits[[method]])) {
   fit<-fits[[method]][[endpoint]];raw<-fits$raw[[endpoint]];keep<-is.finite(fit$effect)&is.finite(raw$effect)
   check(paste(dataset,method,endpoint,'matched feature effects'),identical(fit$feature_id,raw$feature_id))
   effect_rows[[length(effect_rows)+1L]]<-data.frame(dataset=dataset,method_id=method,endpoint=endpoint,features=sum(keep),slope_vs_raw=sum(raw$effect[keep]*fit$effect[keep])/sum(raw$effect[keep]^2),sign_agreement=mean(sign(fit$effect[keep])==sign(raw$effect[keep])),pearson_vs_raw=cor(fit$effect[keep],raw$effect[keep]),significant_features=sum(fit$significant),median_absolute_effect=median(abs(fit$effect[keep])))
  }
  # Same leave-one-batch-out coefficient comparison as the original evaluation;
  # matrices are jointly corrected before the batch-specific regressions.
  for(batch in unique(md$batch))for(method in names(z)) {
   train<-md$batch!=batch;test<-!train;x<-z[[method]]
   if(dataset=='sacurine') {
    dtrain<-droplevels(md[train,]);dtest<-droplevels(md[test,]);a<-model.matrix(~age_z+bmi_z+gender+sampling_factor,dtrain);b<-model.matrix(~age_z+bmi_z+gender+sampling_factor,dtest)
    ep<-lapply(names(coefs),function(endpoint){coef<-coefs[[endpoint]];if(!coef%in%colnames(a)||!coef%in%colnames(b))return(NULL);list(endpoint=endpoint,a=ols_feature_model(x[,train,drop=FALSE],a,match(coef,colnames(a)))$table$coefficient,b=ols_feature_model(x[,test,drop=FALSE],b,match(coef,colnames(b)))$table$coefficient)})
   } else {
    if(dataset=='mtbls79') {
     collapse<-function(selected){u<-unique(md$sample_rep[selected]);m<-vapply(u,function(id)rowMeans(x[,selected & md$sample_rep==id,drop=FALSE]),numeric(nrow(x)));rownames(m)<-rownames(x);list(x=m,group=vapply(u,function(id)unique(md$class[md$sample_rep==id])[1],character(1)))}
     ca<-collapse(train);cb<-collapse(test);a<-two_group_model(ca$x,ca$group,'S','C');b<-two_group_model(cb$x,cb$group,'S','C')
    } else {a<-two_group_model(x[,train,drop=FALSE],md$biological_group[train],'group_0','group_1');b<-two_group_model(x[,test,drop=FALSE],md$biological_group[test],'group_0','group_1')}
    ep<-if(is.null(a)||is.null(b))list()else list(list(endpoint=if(dataset=='mtbls79')'cow_vs_sheep'else'biological_group',a=a$effect,b=b$effect))
   }
   for(e in ep)if(!is.null(e)) {v<-effect_concordance(e$a,e$b);lobo_rows[[length(lobo_rows)+1L]]<-data.frame(dataset=dataset,method_id=method,endpoint=e$endpoint,heldout_batch=batch,pearson=v['pearson'],spearman=v['spearman'],sign_concordance=v['sign_concordance'],n_features=v['n_features'])}
  }
  historical_lobo<-read_input_csv(file.path(derivative,dataset,'leave_one_batch_out_effects.csv'))
  historical_lobo<-historical_lobo[historical_lobo$method_id=='winn_fixed_no_qc' & historical_lobo$feature_panel=='all_common_features',]
  if(!'endpoint'%in%names(historical_lobo))historical_lobo$endpoint<-names(fits$raw)[1]
  historical_lobo$dataset<-dataset
  lobo_rows[[length(lobo_rows)+1L]]<-historical_lobo[,c('dataset','method_id','endpoint','heldout_batch','pearson','spearman','sign_concordance','n_features')]
 }
}
if(length(estimates))write.csv(do.call(rbind,estimates),file.path(out,'frozen_identity_retrieval_estimates.csv'),row.names=FALSE)
if(length(differences))write.csv(do.call(rbind,differences),file.path(out,'frozen_identity_retrieval_paired_differences.csv'),row.names=FALSE)
if(length(bootstrap_rows))write.csv(do.call(rbind,bootstrap_rows),gzfile(file.path(out,'retrieval_query_weighted_bootstrap_replicates.csv.gz')),row.names=FALSE)
if(length(protocol))write.csv(do.call(rbind,protocol),file.path(out,'baseline_identity_protocol_counts.csv'),row.names=FALSE)
if(length(effect_rows))write.csv(do.call(rbind,effect_rows),file.path(out,'baseline_empirical_feature_effects.csv'),row.names=FALSE)
if(length(lobo_rows))write.csv(do.call(rbind,lobo_rows),file.path(out,'baseline_empirical_lobo_effects.csv'),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(out,'baseline_biology_validation_checks.csv'),row.names=FALSE)
write.csv(do.call(rbind,input_manifest),file.path(out,'baseline_empirical_input_manifest.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'baseline_biology_sessionInfo.txt'))
cat('Completed retrieval and baseline effect evaluations with',length(checks),'checks.\n')
