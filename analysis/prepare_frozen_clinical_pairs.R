#!/usr/bin/env Rscript
# Stage authorized frozen pair-level inputs; never include these files in Git.
a <- commandArgs(TRUE)
stopifnot(length(a)==4L)
raw <- readRDS(a[1]); md <- read.csv(a[2],stringsAsFactors=FALSE)
assignments <- read.csv(a[3],stringsAsFactors=FALSE)
out <- a[4]
# The saved Raw method has execution provenance as an attribute, not data.
for (name in setdiff(names(attributes(raw)),c("dim","dimnames"))) attr(raw,name)<-NULL
expected <- "38d77f7fe5f065af7f6b3188c1e5d6f00998327ebba349075e83a226865211cd"
stopifnot(identical(digest::digest(raw,algo="sha256"),expected),
          identical(as.character(md$sample_id),colnames(raw)))
md$within_batch_order <- ave(md$run_order,md$batch,FUN=seq_along)
md$is_qc <- md$role!="clinical";md$is_study <- !md$is_qc
md$biological_label <- md$sample
p <- list(x=raw,meta=md,training_ids=as.character(md$sample_id[md$role=="control"]),
          external_ids=as.character(md$sample_id[md$role=="NIST1950"]),
          hidden_ids=character(),correction_unit="adjacent_dip_pair")
stopifnot(nrow(raw)==1594L,ncol(raw)==660L,length(p$training_ids)==86L,length(p$external_ids)==38L)
assignments <- assignments[assignments$dataset=="clinical_fiams" & assignments$seed_id %in% sprintf("SPLIT%02d",1:10),,drop=FALSE]
stopifnot(nrow(assignments)==1240L,length(unique(assignments$seed_id))==10L,
          all(assignments$sample_id%in%colnames(raw)))
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=","",script_arg))))
expected_counts <- read.csv(file.path(root,"analysis/config/reference_split_clinical_summary.csv"),stringsAsFactors=FALSE)
counts <- aggregate(rep(1L,nrow(assignments)),assignments[c("dataset","seed_id","assignment","batch")],sum)
names(counts)[5] <- "n_units"
key <- function(d) do.call(paste,c(d[c("dataset","seed_id","assignment","batch")],sep="|"))
stopifnot(setequal(key(counts),key(expected_counts)),
          identical(as.integer(counts$n_units[match(key(expected_counts),key(counts))]),as.integer(expected_counts$n_units)))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
saveRDS(p,file.path(out,"prepared_pair_level.rds"),version=3)
write.csv(md,file.path(out,"pair_metadata.csv"),row.names=FALSE)
write.csv(assignments,file.path(out,"reference_split_assignments.csv"),row.names=FALSE)
cat("PASS: exact frozen clinical matrix, 660 pair units and ten reference splits staged.\n")
