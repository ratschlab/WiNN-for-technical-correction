# Can be sourced by the runner or rerun standalone on the saved 4000 seeds.
if(!exists('out')) out<-'generated/final_review_gating_matching'
z<-read.csv(file.path(out,'null_correction_magnitude_replicates.csv'))
stopifnot(nrow(z)==4000,all(z$changed_features==z$selected_features))
z$affected_profile_rmse<-ifelse(z$changed_features>0,z$incremental_log1p_rmse*sqrt(100/z$changed_features),NA_real_)
write.csv(z,file.path(out,'null_correction_magnitude_replicates.csv'),row.names=FALSE)
s<-do.call(rbind,lapply(split(z,list(z$n,z$policy)),function(x)data.frame(n=x$n[1],policy=x$policy[1],replicates=nrow(x),
 probability_any_change=mean(x$changed_features>0),mean_changed_features=mean(x$changed_features),
 mean_incremental_rmse=mean(x$incremental_log1p_rmse),p95_incremental_rmse=quantile(x$incremental_log1p_rmse,.95),
 conditional_mean_rmse=mean(x$incremental_log1p_rmse[x$changed_features>0]),
 mean_rmse_percent_noise_sd=100*mean(x$incremental_log1p_rmse)/.15,
 mean_max_absolute_change=mean(x$max_absolute_log1p_change),
 affected_profile_mean_rmse=mean(x$affected_profile_rmse,na.rm=TRUE),
 affected_profile_p95_rmse=quantile(x$affected_profile_rmse,.95,na.rm=TRUE),
 affected_profile_max_rmse=max(x$affected_profile_rmse,na.rm=TRUE))))
write.csv(s,file.path(out,'null_correction_magnitude_summary.csv'),row.names=FALSE);print(s)
