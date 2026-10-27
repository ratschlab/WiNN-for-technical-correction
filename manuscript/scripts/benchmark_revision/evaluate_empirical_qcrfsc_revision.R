#!/usr/bin/env Rscript
args<-commandArgs(TRUE);opt<-function(k,d=NULL){v<-args[startsWith(args,paste0(k,'='))];if(length(v))sub(paste0('^',k,'='),'',v[1]) else d}
dataset<-opt('--dataset');if(is.null(dataset))stop('Supply --dataset')
all_args<-commandArgs(FALSE);f<-grep('^--file=',all_args,value=TRUE);script<-if(length(f))sub('^--file=','',f[1])else all_args[match('-f',all_args)+1L];here<-dirname(normalizePath(script))
source_root<-here
source(file.path(here,'original_benchmark_helpers.R'));source(file.path(here,'frozen_winn_ablation_helpers.R'));source(file.path(here,'frozen_weighted_pc_r2.R'));source(file.path(here,'frozen_validate_output.R'));source(file.path(here,'frozen_evaluation_helpers.R'))
suppressPackageStartupMessages({library(dplyr);library(tibble)})
p<-readRDS(file.path(opt('--input-root',here),paste0('input_',dataset,'.rds')))
input_results<-opt('--matrix-root',file.path(here,'results'));out<-file.path(opt('--output-root',file.path(here,'empirical_evaluation')),dataset);dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(method in c('legacy','corrected')) {
 filename<-if(method=='legacy')'legacy_reproduced_matrix.rds' else 'corrected_training_qc_abundance.rds'
 x<-readRDS(file.path(input_results,dataset,filename))
 value<-evaluate_release_matrix(dataset,paste0('qc_rfsc_',method),paste('QC-RFSC',method),x,p)
 # Private feature/identity records stay on Euler. Only aggregate metrics enter the manuscript package.
 method_dir<-file.path(out,method);dir.create(method_dir,recursive=TRUE,showWarnings=FALSE)
 for(field in names(value))write.csv(value[[field]],file.path(method_dir,paste0(field,'.csv')),row.names=FALSE)
}
writeLines(c('status=complete',paste0('dataset=',dataset),'scope=frozen endpoint definitions; refitted QC-RFSC only; uncertainty and other-method matrix re-evaluation remain'),file.path(out,'COMPLETE.txt'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
