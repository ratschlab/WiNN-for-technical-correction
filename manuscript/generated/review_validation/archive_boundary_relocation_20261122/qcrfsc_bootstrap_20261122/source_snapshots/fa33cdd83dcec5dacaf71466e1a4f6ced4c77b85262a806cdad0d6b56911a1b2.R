#!/usr/bin/env Rscript
# Conditional uncertainty for the completed statTarget 1.38.0 corrected fits.
# No correction, tuning, feature selection or parameter search is performed.
args <- commandArgs(TRUE)
opt <- function(k, default=NULL) { z <- args[startsWith(args,paste0(k,"="))]; if(length(z)) sub(paste0("^",k,"="),"",z[1]) else default }
here <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(here,"frozen_derivative_biology_helpers.R"))
dataset <- opt("--dataset"); input_root <- opt("--input-root",here)
out <- file.path(opt("--output-root",file.path(here,"bootstrap_20261122")),dataset)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
nboot <- as.integer(opt("--nboot","2000")); cores <- as.integer(opt("--cores","1"))
input_path <- file.path(input_root,paste0("input_",dataset,".rds"))
matrix_path <- file.path(opt("--matrix-root",file.path(here,"results")),dataset,"corrected_training_qc_abundance.rds")
p <- readRDS(input_path); x <- readRDS(matrix_path)
stopifnot(identical(dimnames(x),dimnames(p$x)), all(is.finite(x)), all(x>=0),
          identical(as.character(p$meta$sample_id),colnames(x)))
md <- p$meta[p$meta$is_study,,drop=FALSE]
z <- log1p(x[,md$sample_id,drop=FALSE])
seed <- 2026072701L + switch(dataset, mtbls79=79L, sacurine=113L, batchcorr_set1=199L, clinical_fiams=1594L, waveica=0L)
pc <- function(v,g,type="categorical") pca_r2_from_log(v,g,type)
if(dataset %in% c("batchcorr_set1","clinical_fiams")) {
  field <- if(dataset=="batchcorr_set1") "accession_id" else "sample"
  identities <- as.character(md[[field]]); units <- unique(identities)
  index <- split(seq_len(nrow(md)),identities)
  reps <- replicate_identity_metrics(z,md,field)
  matching <- cross_batch_matching(z,md,field)
  # Current manuscript estimand is query weighted, not equal identity weighted.
  counts <- table(matching$query$identity)
  weights <- as.numeric(counts[match(units,names(counts))]); weights[is.na(weights)] <- 0
  mt <- matching$identity[match(units,matching$identity$identity),]
  mt[weights==0,c("top1","reciprocal_rank")] <- 0
  draws <- cluster_bootstrap_draws(units,NULL,nboot,seed)
  metric <- function(draw) {
    j <- match(draw,units); cols <- unlist(index[draw],use.names=FALSE)
    r <- reps[match(draw,reps$identity),]
    c(batch_weighted_pc_r2_categorical=pc(z[,cols,drop=FALSE],factor(md$batch[cols])),
      sample_replicate_pearson_median=median(r$pearson,na.rm=TRUE),
      genuine_replicate_icc_median=median(r$icc_a1,na.rm=TRUE),
      cross_batch_top1_accuracy=sum(weights[j]*mt$top1[j])/sum(weights[j]),
      cross_batch_mean_reciprocal_rank=sum(weights[j]*mt$reciprocal_rank[j])/sum(weights[j]))
  }
  point <- metric(units)
  # Save aggregate candidate geometry, never clinical identity-level data.
  sim <- cor(z); diag(sim)<-NA_real_
  qi <- matching$query$query_index
  cn <- vapply(qi,function(i) sum(md$batch != md$batch[i] & is.finite(sim[i,])),integer(1))
  chance <- vapply(qi,function(i) { j<-which(md$batch!=md$batch[i] & is.finite(sim[i,])); mean(identities[j]==identities[i]) },numeric(1))
  ties <- vapply(qi,function(i) {j<-which(md$batch!=md$batch[i] & is.finite(sim[i,]));sum(sim[i,j]==max(sim[i,j]))>1L},logical(1))
  write.csv(data.frame(dataset_key=dataset,method_id="qc_rfsc",eligible_queries=length(qi),eligible_identities=sum(weights>0),all_identities=length(units),candidate_min=min(cn),candidate_median=median(cn),candidate_max=max(cn),chance_top1=mean(chance),exact_top_ties=sum(ties)),file.path(out,"matching_protocol.csv"),row.names=FALSE)
} else if(dataset=="mtbls79") {
  units <- unique(md$sample_rep)
  group <- vapply(units,function(u) { g<-unique(md$class[md$sample_rep==u]);stopifnot(length(g)==1L);g },character(1))
  index <- split(seq_len(nrow(md)),md$sample_rep)
  draws <- cluster_bootstrap_draws(units,group,nboot,seed)
  metric <- function(draw) {
    cols <- unlist(index[draw],use.names=FALSE)
    means <- vapply(draw,function(u)rowMeans(z[,index[[u]],drop=FALSE]),numeric(nrow(z)))
    c(batch_weighted_pc_r2_categorical=pc(z[,cols,drop=FALSE],factor(md$batch[cols])),
      biological_weighted_pc_r2=pc(means,factor(group[match(draw,units)])))
  }
  point <- metric(units)
} else if(dataset=="waveica") {
  units <- as.character(md$sample_id); group <- factor(md$biological_group,levels=c("group_0","group_1"))
  draws <- cluster_bootstrap_draws(units,interaction(md$batch,group,drop=TRUE),nboot,seed)
  metric <- function(draw) {
    j<-match(draw,units); basis<-pca_basis_from_log(z[,j,drop=FALSE])
    c(batch_weighted_pc_r2_categorical=pca_r2_from_basis(basis,factor(md$batch[j]),"categorical"),
      biological_weighted_pc_r2=pca_r2_from_basis(basis,factor(group[j]),"categorical"))
  }
  point <- metric(units)
} else if(dataset=="sacurine") {
  units <- as.character(md$sample_id)
  biological <- log1p(sweep(x[,md$sample_id,drop=FALSE],2,md$osmolality,"/"))
  draws <- cluster_bootstrap_draws(units,md$batch,nboot,seed)
  metric <- function(draw) {
    j<-match(draw,units); raw_basis<-pca_basis_from_log(z[,j,drop=FALSE]); bio_basis<-pca_basis_from_log(biological[,j,drop=FALSE])
    c(batch_weighted_pc_r2_categorical=pca_r2_from_basis(raw_basis,factor(md$batch[j]),"categorical"),
      age_weighted_pc_r2=pca_r2_from_basis(bio_basis,md$age[j],"continuous"),
      bmi_weighted_pc_r2=pca_r2_from_basis(bio_basis,md$bmi[j],"continuous"),
      gender_weighted_pc_r2=pca_r2_from_basis(bio_basis,factor(md$gender[j]),"categorical"))
  }
  point <- metric(units)
} else stop("Unsupported dataset")
write.csv(data.frame(dataset_key=dataset,method_id="qc_rfsc",metric=names(point),point_estimate=point),file.path(out,"point_estimates.csv"),row.names=FALSE)
# Save checksums before work; draws remain private because they contain IDs.
write.csv(data.frame(dataset_key=dataset,object=c("prepared_input","corrected_matrix","bootstrap_draws","helper_source","evaluation_source"),sha256=c(vapply(c(input_path,matrix_path),function(f)digest::digest(file=f,algo="sha256"),character(1)),digest::digest(draws,algo="sha256"),vapply(file.path(here,c("frozen_derivative_biology_helpers.R","bootstrap_empirical_qcrfsc.R")),function(f)digest::digest(file=f,algo="sha256"),character(1))),nboot=nboot,seed=seed),file.path(out,"input_manifest.csv"),row.names=FALSE)
values <- parallel::mclapply(seq_along(draws),function(i) metric(draws[[i]]),mc.cores=cores,mc.preschedule=TRUE,mc.set.seed=FALSE)
stopifnot(all(vapply(values,is.numeric,logical(1))))
bs <- do.call(rbind,values)
stopifnot(nrow(bs)==nboot, all(is.finite(bs)))
replicates <- do.call(rbind,lapply(seq_along(point),function(j)data.frame(dataset_key=dataset,method_id="qc_rfsc",bootstrap=seq_len(nboot),metric=names(point)[j],value=bs[,j])))
summary <- do.call(rbind,lapply(seq_along(point),function(j)data.frame(dataset_key=dataset,method_id="qc_rfsc",method="QC-RFSC",metric=names(point)[j],direction=if(names(point)[j]=="batch_weighted_pc_r2_categorical")"lower" else "higher",point_estimate=point[j],bootstrap_mean=mean(bs[,j]),bootstrap_se=sd(bs[,j]),ci_lower=unname(quantile(bs[,j],.025)),ci_upper=unname(quantile(bs[,j],.975)),valid_bootstrap=nboot)))
write.csv(replicates,gzfile(file.path(out,"bootstrap_replicates.csv.gz")),row.names=FALSE)
write.csv(summary,file.path(out,"bootstrap_estimates.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
writeLines(c("status=complete","correction_refitted=FALSE","tuning_repeated=FALSE",paste0("bootstrap_replicates=",nboot),paste0("bootstrap_seed=",seed)),file.path(out,"COMPLETE.txt"))
print(summary)
