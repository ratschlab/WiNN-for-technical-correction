#!/usr/bin/env Rscript
args<-commandArgs(FALSE);here<-dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1])))
root<-dirname(dirname(here));src<-file.path(root,'generated/review_validation');dest<-file.path(root,'generated/supplementary_tables')
labels<-c(raw='Raw',combat='ComBat',qc_rlsc='QC-RLSC',tiger='TIGER',serrf='SERRF',winn_auto_qc='WiNN (QC-tuned)',winn_auto_batch_qc='WiNN (inferred segments)',winn_fixed_no_qc='WiNN',mean_batch_all='Mean-only batch',wehrens_ancova_all='Study-sample ANCOVA')
ci<-function(point,lo,hi)sprintf('%.3f [%.3f, %.3f]',ifelse(abs(point)<.0005,0,point),ifelse(abs(lo)<.0005,0,lo),ifelse(abs(hi)<.0005,0,hi))
start<-function(columns,caption,label,header,landscape=FALSE)c(if(landscape)'\\begin{landscape}', '\\footnotesize',paste0('\\begin{longtable}{',columns,'}'),paste0('\\caption{',caption,'}\\label{',label,'}\\\\'),'\\toprule',paste0(header,' \\\\'),'\\midrule','\\endfirsthead','\\toprule',paste0(header,' \\\\'),'\\midrule','\\endhead','\\bottomrule','\\endlastfoot')
finish<-function(lines,path,landscape=FALSE)writeLines(c(lines,'\\end{longtable}',if(landscape)'\\end{landscape}'),file.path(dest,path))
a<-read.csv(file.path(src,'frozen_identity_retrieval_estimates.csv'));d<-read.csv(file.path(src,'frozen_identity_retrieval_paired_differences.csv'))
s<-start('p{2.7cm}p{3.4cm}p{4.3cm}p{4.3cm}p{4.8cm}', 'Cross-batch identity discrimination. Estimates and paired biological-identity bootstrap percentile intervals retain eligible-query counts within every resampled identity (2,000 draws). The correlation candidate set and fitted corrections remain fixed. BatchCorr uses 693 eligible queries from 346 identities; clinical FIA-MS uses 536 correction-unit queries from 267 identities. Top-1 accuracy and mean reciprocal rank (MRR) use Pearson similarity of log1p profiles in other batches. The final column is WiNN minus comparator.','tab:supp_identity_retrieval','Dataset & Method & Top-1 accuracy [95\\% interval] & MRR [95\\% interval] & Paired top-1 difference [95\\% interval]',TRUE)
for(dataset in c('batchcorr_set1','clinical_fiams'))for(method in names(labels)){
 q<-a[a$dataset==dataset & a$method_id==method,];top<-q[q$metric=='cross_batch_top1_accuracy',];mrr<-q[q$metric=='cross_batch_mean_reciprocal_rank',];delta<-d[d$dataset==dataset & d$comparator_method_id==method & d$metric=='cross_batch_top1_accuracy',]
 if(!nrow(top))next
 s<-c(s,paste0(paste(c(if(dataset=='batchcorr_set1')'BatchCorr Set 1'else'Clinical FIA-MS',labels[[method]],ci(top$point_estimate,top$ci_lower,top$ci_upper),ci(mrr$point_estimate,mrr$ci_lower,mrr$ci_upper),if(method=='winn_fixed_no_qc')'Reference'else ci(delta$point_direction_normalized_difference,delta$ci_lower,delta$ci_upper)),collapse=' & '),' \\\\'))
}
finish(s,'table_s23_identity_retrieval.tex',TRUE)
a<-read.csv(file.path(src,'baseline_simulation_summary.csv'));d<-read.csv(file.path(src,'baseline_simulation_paired_differences.csv'))
s<-start('p{3.3cm}p{2.4cm}p{4.1cm}p{4.1cm}p{4.1cm}', 'Same-information baseline validation on the ten canonical simulation seeds. Mean-only batch correction and study-sample ANCOVA use all injections without reference or biological labels. Upper rows give mean $\\pm$ standard error across seeds. Lower rows give paired WiNN advantages and 95\\% $t_9$ intervals; positive values favor WiNN (comparator minus WiNN for RMSE, WiNN minus comparator for ICC and correlation). WiNN improves all three truth endpoints in every paired seed.','tab:supp_baseline_truth','Method & Hidden-reference CV (\\%) & Truth log1p RMSE & Feature ICC & Sample correlation',TRUE)
for(method in c('raw','mean_batch_all','wehrens_ancova_all','winn_fixed_no_qc')){
 metrics<-c('heldout_qc_cv_mean','truth_log1p_rmse','truth_metabolite_profile_icc_mean','truth_sample_profile_pearson_mean')
 vals<-vapply(metrics,function(metric){z<-a[a$method_id==method & a$metric==metric,];sprintf(if(metric=='heldout_qc_cv_mean')'$%.3f \\pm %.3f$'else'$%.4f \\pm %.4f$',z$mean,z$se)},character(1))
 s<-c(s,paste0(paste(c(labels[[method]],vals),collapse=' & '),' \\\\'))
}
s<-c(s,'\\midrule','\\multicolumn{5}{l}{Paired WiNN advantage over comparator [95\\% interval]} \\\\')
for(method in c('mean_batch_all','wehrens_ancova_all')){
 vals<-vapply(c('truth_log1p_rmse','truth_metabolite_profile_icc_mean','truth_sample_profile_pearson_mean'),function(metric){z<-d[d$comparator==method & d$metric==metric,];ci(z$mean_winn_advantage,z$ci_lower,z$ci_upper)},character(1))
 s<-c(s,paste0(paste(c(labels[[method]],'---',vals),collapse=' & '),' \\\\'))
}
finish(s,'table_s27_baseline_simulation.tex',TRUE)
a<-read.csv(file.path(src,'baseline_empirical_feature_effects.csv'));l<-read.csv(file.path(src,'baseline_empirical_lobo_effects.csv'))
s<-start('p{3cm}p{3.4cm}rrrrr', 'Same-information baseline comparisons of empirical feature effects. Effects are matched by feature against raw input; slope is the through-origin slope, sign agreement is the fraction with the same sign, and Pearson correlation compares effect vectors. MTBLS79 models use 20 independent extracts after averaging repeated injections. Sacurine models use log1p intensity after osmolality division and simultaneous age, BMI, sex and sampling covariates. LOBO reports median train--held-out batch effect correlation, with correction fitted jointly across batches. FDR counts use Benjamini--Hochberg at 0.05 within each method and endpoint. WiNN coefficients and LOBO results come from the frozen derivative export. Raw effects are a descriptive reference, not known biological truth.','tab:supp_baseline_effects','Dataset / endpoint & Method & Slope & Sign agreement & Pearson & LOBO Pearson & FDR count',TRUE)
for(i in seq_len(nrow(a))){z<-a[i,];lo<-l[l$dataset==z$dataset & l$method_id==z$method_id & l$endpoint==z$endpoint,];name<-if(z$dataset=='mtbls79')'MTBLS79 cow--sheep'else if(z$dataset=='waveica')'WaveICA groups'else paste('Sacurine',z$endpoint)
 s<-c(s,paste0(paste(c(name,labels[[z$method_id]],sprintf('%.3f',c(z$slope_vs_raw,z$sign_agreement,z$pearson_vs_raw,median(lo$pearson,na.rm=TRUE))),z$significant_features),collapse=' & '),' \\\\'))}
finish(s,'table_s28_baseline_empirical_effects.tex',TRUE)
