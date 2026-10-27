# Integration checks with the real statTarget REGfit training-index implementation.
args <- commandArgs(FALSE);file_arg<-grep('^--file=',args,value=TRUE);script<-if(length(file_arg))sub('^--file=','',file_arg[1]) else args[match('-f',args)+1L];here<-dirname(normalizePath(script))
source(file.path(here,'original_benchmark_helpers.R'))
source(file.path(here,'qcrfsc_review_adapter.R'))
suppressPackageStartupMessages(library(dplyr))
set.seed(281)
n<-40L;f<-12L
x<-outer(exp(seq(2,9,length.out=f)),exp(rnorm(n,0,.03))) * matrix(exp(rnorm(f*n,0,.01)),f,n)
rownames(x)<-paste0('metabolite',seq_len(f));colnames(x)<-paste0('sample',seq_len(n))
md<-data.frame(sample_id=colnames(x),batch=rep(1:2,each=20),run_order=seq_len(n))
controls<-colnames(x)[c(1,5,10,15,20,21,25,30,35,40)]
.adversarial_checked <- 0L
invisible(capture.output(trace('REGfit',where=asNamespace('statTarget'),print=FALSE,tracer=quote({
  used<-colnames(x)[grepl('QC|qc|Qc',colnames(x))]
  allowed<-intersect(colnames(x),getOption('qcrfsc_authorized_names'))
  stopifnot(setequal(used,allowed))
  assign('.adversarial_checked',get('.adversarial_checked',envir=.GlobalEnv)+1L,envir=.GlobalEnv)
}))))
run<-function(x,md,controls,scale='training_qc_abundance'){
 set.seed(9124)
 invisible(capture.output(y<-run_qc_rfsc_review(x,controls,md,ntree=30,coCV=100,Frule=.8,output_scale=scale)))
 y
}
a<-run(x,md,controls)
stopifnot(is.matrix(a),identical(dimnames(a),dimnames(x)),all(is.finite(a)),all(a>=0))
# Change only held-out/study identifiers to strings that statTarget normally treats as QC.
adversarial_names<-colnames(x);adversarial_names[c(3,7,23,28)]<-c('hidden_QC03','qc_patient7','Qc_reference23','ordinary_QC28')
xx<-x;colnames(xx)<-adversarial_names;mm<-md;mm$sample_id<-adversarial_names
b<-run(xx,mm,adversarial_names[match(controls,colnames(x))])
stopifnot(isTRUE(all.equal(unname(a),unname(b),check.attributes=FALSE,tolerance=0)))
map<-attr(b,'qcrfsc_sample_id_map')
stopifnot(!anyDuplicated(map$internal_id),setequal(map$original_id,colnames(xx)),identical(grepl('QC|qc|Qc',map$internal_id),map$training_reference))
# Scale restoration uses only raw training-QC feature medians.
r<-run(x,md,controls,'relative_ratio')
expected<-r*apply(x[,controls,drop=FALSE],1,median)
stopifnot(isTRUE(all.equal(unname(a),unname(expected),check.attributes=FALSE,tolerance=1e-12)))
stopifnot(.adversarial_checked>0L,identical(sort(attr(a,'qcrfsc_training_ids')),sort(controls)))
untrace('REGfit',where=asNamespace('statTarget'))
cat('PASS: real REGfit uses exactly authorized controls; adversarial hidden IDs cannot leak; mapping and dimensions preserved; non-oracle abundance restoration verified.\n')
cat('statTarget version:',as.character(packageVersion('statTarget')),'; REGfit calls asserted:',.adversarial_checked,'\n')
