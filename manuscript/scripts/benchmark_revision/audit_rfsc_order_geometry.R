#!/usr/bin/env Rscript
# Aggregate acquisition-coordinate audit of the documented boundary-QC adapter.
args<-commandArgs(TRUE)
opt<-function(k){v<-args[startsWith(args,paste0(k,"="))];stopifnot(length(v)==1);sub(paste0("^",k,"="),"",v)}
input<-opt("--input-root"); output<-opt("--output")
rows<-list()
for(ds in c("mtbls79","clinical_fiams","batchcorr_set1","sacurine","waveica")) {
  p<-readRDS(file.path(input,paste0("input_",ds,".rds")))
  md<-p$meta[order(p$meta$run_order),,drop=FALSE]; ids<-as.character(md$sample_id)
  scopes<-list(full_training=p$training_ids)
  for(fold in p$folds$rows) scopes[[paste0("tuning_",fold$fold_id)]]<-fold$fit_qc_ids[[1]]
  for(scope in names(scopes)) {
    controls<-intersect(scopes[[scope]],ids); old_qc<-which(ids%in%controls)
    first<-min(old_qc);last<-max(old_qc)
    moved<-!(1L%in%old_qc && length(ids)%in%old_qc)
    reordered<-if(moved)c(first,setdiff(seq_along(ids),c(first,last)),last)else seq_along(ids)
    new_position<-match(seq_along(ids),reordered)
    rows[[length(rows)+1L]]<-data.frame(dataset_key=ds,scope=scope,injections=length(ids),training_qcs=length(controls),first_qc_position=first,last_qc_position=last,relocation_used=moved,relocated_reference_count=if(moved)sum(c(first!=1L,last!=length(ids)))else 0L,positions_changed=sum(new_position!=seq_along(ids)),max_abs_position_change=max(abs(new_position-seq_along(ids))),training_positions_changed=sum(new_position[old_qc]!=old_qc),study_positions_changed=sum(new_position[md$is_study]!=which(md$is_study)))
  }
}
write.csv(do.call(rbind,rows),output,row.names=FALSE)
