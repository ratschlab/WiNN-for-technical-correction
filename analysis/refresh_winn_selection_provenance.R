#!/usr/bin/env Rscript
# Correct descriptive selection provenance for the new installed package.
# This does not read or alter corrected matrices, scores, or selected candidates.
a<-commandArgs(TRUE);stopifnot(length(a)==1L);root<-normalizePath(a[1])
config<-read.csv(file.path(root,"analysis/config/endpoint_free_selection_manifest.csv"),stringsAsFactors=FALSE,check.names=FALSE)
rows<-list()
for(family in c("primary","simulation_seeds","reference_splits")) {
 paths<-list.files(file.path(root,"results",family),pattern="^selected_configuration[.]csv$",recursive=TRUE,full.names=TRUE)
 for(path in paths) {
  x<-read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
  if(!all(c("source_version","method_id","dataset_key")%in%names(x)) || !nrow(x) || !all(grepl("^winn_",x$method_id)))next
  before<-digest::digest(path,file=TRUE,algo="sha256")
  for(i in seq_len(nrow(x))) {
   y<-config[config$method_id==x$method_id[i]&config$dataset_key==x$dataset_key[i],]
   stopifnot(nrow(y)==1L)
   x$source_version[i]<-y$source_version
   x$candidate_count[i]<-y$candidate_count
   x$selected_parameters[i]<-y$selected_parameters
  }
  write.csv(x,path,row.names=FALSE,quote=TRUE)
  rows[[length(rows)+1L]]<-data.frame(path=substring(path,nchar(root)+2L),before_sha256=before,after_sha256=digest::digest(path,file=TRUE,algo="sha256"),scope="WiNN descriptive package version, calibration settings and nominal grid count; no fitted outputs changed")
 }
}
out<-file.path(root,"results/final/validation");dir.create(out,recursive=TRUE,showWarnings=FALSE)
record<-file.path(out,"selection_provenance_refresh.csv")
initial<-file.path(out,"selection_provenance_refresh_initial.csv")
# A retry must not erase the first before/after provenance transition.
if(file.exists(record)&&!file.exists(initial))file.copy(record,initial)
write.csv(do.call(rbind,rows),record,row.names=FALSE)
if(!file.exists(initial))file.copy(record,initial)
