#!/usr/bin/env Rscript
# Evaluate the saved corrected RFSC fit with the manuscript's sample-mean estimand.
args <- commandArgs(TRUE)
opt <- function(key) { v <- args[startsWith(args, paste0(key, "="))]; if (length(v)!=1L) stop("Supply ",key); sub(paste0("^",key,"="),"",v) }
full <- commandArgs(FALSE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=",full,value=TRUE)))))
source(file.path(root,"manuscript/scripts/nearest_baseline/derivative_biology_helpers.R"))
input_path <- opt("--input"); matrix_path <- opt("--matrix"); output <- opt("--output")
p <- readRDS(input_path); x <- readRDS(matrix_path)
stopifnot(identical(dimnames(x),dimnames(p$x)), all(is.finite(x)), all(x>=0))
md <- p$meta[match(colnames(x),p$meta$sample_id),,drop=FALSE]
stopifnot(identical(as.character(md$sample_id),colnames(x)))
md <- md[md$is_study,,drop=FALSE]; units <- unique(md$sample_rep)
group <- vapply(units,function(u) {g<-unique(md$class[md$sample_rep==u]);stopifnot(length(g)==1L);g},character(1))
stopifnot(nrow(md)==134L,length(units)==20L,all(table(group)==10L))
rows <- lapply(c("raw","qc_rfsc"),function(method) {
  value <- if(method=="raw")p$x else x
  z <- log1p(value[,md$sample_id,drop=FALSE])
  means <- vapply(units,function(u)rowMeans(z[,md$sample_rep==u,drop=FALSE]),numeric(nrow(z)))
  data.frame(dataset_key="mtbls79",method_id=method,
    batch_pc_variance_percent=100*pca_r2_from_log(z,as.character(md$batch),"categorical"),
    biological_pc_variance_percent=100*pca_r2_from_log(means,group,"categorical"),
    batch_n=ncol(z),biological_n=ncol(means),features=nrow(z),
    biological_unit="mean log1p profile per sample extract",uncertainty="not recomputed")
})
result <- do.call(rbind,rows)
reference <- read.csv(file.path(root,"manuscript/generated/revision_panel_inputs/eccb2026_public_effect_tradeoffs_v2.csv"))
raw <- reference[reference$dataset_key=="mtbls79" & reference$method_id=="raw",]
stopifnot(nrow(raw)==1L,abs(result$batch_pc_variance_percent[1]-raw$batch_pc_variance_percent)<1e-8,
          abs(result$biological_pc_variance_percent[1]-raw$biological_pc_variance_percent)<1e-8)
previous <- read.csv(file.path(root,"manuscript/generated/review_validation/qcrfsc_completed_empirical_all_endpoints.csv"))
old_batch <- previous$value[previous$dataset_key=="mtbls79" & previous$metric=="batch_weighted_pc_r2_categorical"]
stopifnot(length(old_batch)==1L,abs(result$batch_pc_variance_percent[2]-100*old_batch)<1e-8)
dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
write.csv(result,output,row.names=FALSE)
write.csv(data.frame(input=c("prepared_input","corrected_matrix"),sha256=vapply(c(input_path,matrix_path),function(f)digest::digest(file=f,algo="sha256"),character(1))),sub("[.]csv$","_input_hashes.csv",output),row.names=FALSE)
print(result)
cat("PASS: raw endpoints reproduce manuscript points; corrected batch endpoint reproduces archived RFSC evaluation.\n")
