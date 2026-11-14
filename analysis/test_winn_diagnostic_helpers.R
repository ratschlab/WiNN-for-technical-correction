#!/usr/bin/env Rscript
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),".."))
source(file.path(root,"scripts/public_benchmark_method_helpers.R"))
stopifnot(parse_winn_auto_messages("Autocorrelation test: Ljung-Box (FDR: 0.01)")$acorr_fdr==0.01)
stopifnot(parse_winn_auto_messages("Autocorrelation test: Ljung-Box (gate threshold: 0.1)")$acorr_fdr==0.1)
stopifnot(is.na(parse_winn_auto_messages("No selected threshold")$acorr_fdr))
stopifnot(utils::packageVersion("winn")=="0.1.7.9999")
set.seed(20261114)
x<-matrix(rexp(20L*40L),20L,dimnames=list(paste0("F",1:20),paste0("S",1:40)))
batch<-rep(c("A","B"),each=20L);order<-rep(20:1,2)
x[1,]<-exp(rep(seq(-1,2,length.out=20),2))
for(calibration in c("permutation","asymptotic")) {
 expected<-winn::autocorrelation_correct(x,batch=batch,run_order=order,return_diagnostics=TRUE,
   fdr_threshold=0.1,drift_calibration=calibration,n_permutations=49L,permutation_seed=104729L)$diagnostics
 observed<-winn_drift_test_table(x,batch,order,"Ljung-Box",0.1,
   drift_calibration=calibration,n_permutations=49L,permutation_seed=104729L)
 expected<-expected[match(paste(observed$feature_id,observed$batch),paste(expected$feature_id,expected$segment)),]
 stopifnot(isTRUE(all.equal(observed$p_adj,expected$adjusted_p_value,tolerance=0)),
   identical(observed$selected_for_drift,expected$selective_significance))
}
cat("PASS: legacy/current threshold parsing and gate diagnostics agree with correction.\n")
