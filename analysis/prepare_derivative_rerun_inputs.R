#!/usr/bin/env Rscript
# Use new WiNN matrices and the exact previous Raw/comparator files.
a <- commandArgs(TRUE); stopifnot(length(a)==1L)
root <- normalizePath(a[1]); base <- file.path(root,"results/frozen_endpoint_free_inputs")
registry <- read.csv(file.path(root,"analysis/config/method_manifest.csv"),stringsAsFactors=FALSE)
methods <- c("raw","combat","qc_rlsc","qc_rfsc","tiger","serrf","winn_auto_qc","winn_auto_batch_qc","winn_fixed_no_qc")
datasets <- c("simulation","mtbls79","clinical_fiams","batchcorr_set1","sacurine","waveica")
old <- read.csv(file.path(base,"provenance/matrix_checksums_pre.csv"),stringsAsFactors=FALSE)
rows <- list();k<-0L
for(d in datasets)for(m in methods) {
 rel <- file.path("results/primary",d,m,"corrected_matrix.rds")
 target <- file.path(base,"primary",d,m,"corrected_matrix.rds")
 if(grepl("^winn_",m)) {
  src<-file.path(root,rel);stopifnot(file.exists(src));dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  stopifnot(file.copy(src,target,overwrite=TRUE))
 } else {
  expected<-old$sha256[match(rel,old$path)]
  stopifnot(length(expected)==1L,!is.na(expected),identical(digest::digest(file=target,algo="sha256"),expected))
 }
 k<-k+1L;rows[[k]]<-data.frame(category="primary",path=rel,bytes=file.info(target)$size,sha256=digest::digest(file=target,algo="sha256"))
}
manifest<-do.call(rbind,rows)
stopifnot(nrow(manifest)==54L)
write.csv(manifest,file.path(base,"provenance/matrix_checksums_pre.csv"),row.names=FALSE)
unlink(file.path(base,"provenance/matrix_checksums_post.csv"))
cat("Staged 18 new WiNN matrices; 36 Raw/comparator file hashes preserved.\n")
