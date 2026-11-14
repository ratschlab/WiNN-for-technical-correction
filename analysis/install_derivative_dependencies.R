#!/usr/bin/env Rscript
# Install the archived PCA solver used by the frozen derivative workflow.
args<-commandArgs(TRUE)
lib<-if(length(args))args[1] else Sys.getenv("R_LIBS_USER", "Rlib")
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
lib<-normalizePath(lib)
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),".."))
manifest<-read.csv(file.path(root,"analysis/config/calibration_revision/derivative_dependencies.csv"),stringsAsFactors=FALSE)
stopifnot(requireNamespace("digest",quietly=TRUE))
for(i in seq_len(nrow(manifest))) {
 row<-manifest[i,]
 if(requireNamespace(row$package,quietly=TRUE) && as.character(utils::packageVersion(row$package))==row$version)next
 archive<-tempfile(fileext=".tar.gz")
 utils::download.file(row$source_archive_url,archive,mode="wb",quiet=TRUE)
 stopifnot(identical(digest::digest(archive,file=TRUE,algo="sha256"),row$sha256))
 utils::install.packages(archive,repos=NULL,type="source",lib=lib)
 stopifnot(as.character(utils::packageVersion(row$package,lib.loc=lib))==row$version)
 unlink(archive)
}
cat("Pinned derivative dependencies are available.\n")
