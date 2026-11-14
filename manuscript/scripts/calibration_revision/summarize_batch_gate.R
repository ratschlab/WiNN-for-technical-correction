#!/usr/bin/env Rscript
# Merge new after-drift checks and assert unchanged pre-drift ANOVA results.
a<-commandArgs(TRUE);stopifnot(length(a)==3L)
files<-list.files(a[1],'^batch_gate_[0-9]+[.]csv$',full.names=TRUE)
d<-do.call(rbind,lapply(files,read.csv));stopifnot(nrow(d)==800L)
keys<-c('seed','repeated_identities','technical_batch_artifact','stage')
stopifnot(!anyDuplicated(d[keys]))
old<-read.csv(a[2]);before<-merge(old[old$stage=='before_drift',],
 d[d$stage=='before_drift',],by=keys,suffixes=c('_old','_new'))
stopifnot(nrow(before)==400L)
for(field in c('discoveries','sensitivity','false_discovery_proportion','null_p_lt_005')) {
 x<-before[[paste0(field,'_old')]];y<-before[[paste0(field,'_new')]]
 stopifnot(identical(is.na(x),is.na(y)),all(abs(x-y)<1e-12,na.rm=TRUE))
}
dir.create(a[3],recursive=TRUE,showWarnings=FALSE)
write.csv(d,file.path(a[3],'batch_gate_calibration_replicates.csv'),row.names=FALSE)
s<-aggregate(d[c('discoveries','sensitivity','false_discovery_proportion','null_p_lt_005')],
 d[c('repeated_identities','technical_batch_artifact','stage')],
 function(x)if(all(is.na(x)))NA_real_ else mean(x,na.rm=TRUE))
write.csv(s,file.path(a[3],'batch_gate_calibration_summary.csv'),row.names=FALSE)
cat('800 rows complete; all400 pre-drift rows reproduce the frozen ANOVA check.\n')
