#!/usr/bin/env Rscript
# Test the public entry point against the single canonical chronological adapter.
a <- commandArgs(FALSE); here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", a, value=TRUE))))
root <- dirname(here)
loaded_before <- "statTarget" %in% loadedNamespaces()
source(file.path(root,"scripts/benchmark_helpers.R"))
stopifnot(identical("statTarget" %in% loadedNamespaces(), loaded_before))
stopifnot(as.character(packageVersion("statTarget"))=="1.38.0")
native_before <- digest::digest(list(body(getFromNamespace("REGfit","statTarget")),body(statTarget::shiftCor)),algo="sha256")
canonical <- new.env(parent=globalenv())
sys.source(file.path(root,"manuscript/scripts/benchmark_revision/qcrfsc_review_adapter.R"),canonical)
set.seed(281); n <- 40L; f <- 12L
x <- outer(exp(seq(2,9,length.out=f)),exp(rnorm(n,0,.03))) * matrix(exp(rnorm(f*n,0,.01)),f,n)
rownames(x) <- paste0("metabolite",seq_len(f)); colnames(x) <- paste0("sample",seq_len(n))
md <- data.frame(sample_id=colnames(x),batch=rep(1:2,each=20),run_order=seq_len(n))
controls <- colnames(x)[c(5,10,15,20,21,25,30,35)]
.rfsc_inner_audit <- list()
run <- function(fn,x,md,ids,scale="training_qc_abundance") {
  set.seed(9124)
  invisible(capture.output(y <- fn(x,ids,md,ntree=30,coCV=100,Frule=.8,output_scale=scale)))
  y
}
a <- run(run_qc_rfsc_with_controls,x,md,controls)
b <- run(canonical$run_qc_rfsc_review,x,md,controls)
stopifnot(isTRUE(all.equal(a,b,tolerance=0)), identical(dimnames(a),dimnames(x)),all(is.finite(a)),all(a>=0),
  identical(attr(a,"qcrfsc_sample_id_map")$original_id,as.character(md$sample_id)),
  attr(a,"qcrfsc_outside_training_range_count")==9L)
names2 <- colnames(x); names2[c(3,7,23,28)] <- c("hidden_QC03","qc_patient7","Qc_reference23","ordinary_QC28")
xx <- x; colnames(xx) <- names2; mm <- md; mm$sample_id <- names2
renamed <- run(run_qc_rfsc_with_controls,xx,mm,names2[match(controls,colnames(x))])
stopifnot(isTRUE(all.equal(unname(a),unname(renamed),check.attributes=FALSE,tolerance=0)))
r <- run(run_qc_rfsc_with_controls,x,md,controls,"relative_ratio")
stopifnot(isTRUE(all.equal(unname(a),unname(r*apply(x[,controls,drop=FALSE],1,median)),check.attributes=FALSE,tolerance=1e-12)))
# Fully bracketed input must reproduce the unmodified native algorithm exactly.
bracketed <- colnames(x)[c(1,5,10,15,20,21,25,30,35,40)]
full <- run(run_qc_rfsc_with_controls,x,md,bracketed)
canonical$.qcrfsc_shiftCor <- statTarget::shiftCor
native <- run(canonical$run_qc_rfsc_review,x,md,bracketed)
stopifnot(isTRUE(all.equal(unname(full),unname(native),check.attributes=FALSE,tolerance=0)),
  length(.rfsc_inner_audit)>0L,all(vapply(.rfsc_inner_audit,function(z)z$chronological,logical(1))),
  native_before==digest::digest(list(body(getFromNamespace("REGfit","statTarget")),body(statTarget::shiftCor)),algo="sha256"))
cat("PASS: lazy public wrapper; exact canonical equivalence; chronological positions with unbracketed controls; hidden-ID invariance; training-QC abundance scale; native equivalence for bracketed controls; installed statTarget unchanged.\n")
