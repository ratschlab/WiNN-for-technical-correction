#!/usr/bin/env Rscript
a<-commandArgs(TRUE);stopifnot(length(a)==1L);root<-normalizePath(a[1]);base<-file.path(root,"results/final")
p<-list.files(base,recursive=TRUE,full.names=TRUE)
p<-p[!file.info(p)$isdir & !grepl("/derivative/",p) & basename(p)!="file_manifest.csv"]
x<-data.frame(path=substring(p,nchar(base)+2L),bytes=file.info(p)$size,sha256=vapply(p,digest::digest,character(1),file=TRUE,algo="sha256"))
write.csv(x,file.path(root,"results/frozen_endpoint_free_inputs/provenance/pre_derivative_file_manifest.csv"),row.names=FALSE)
