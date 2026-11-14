#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args)==4L)
source_root <- normalizePath(args[1]); out <- args[2]; task <- as.integer(args[3]); core_dir <- args[4]
source(file.path(source_root,'scripts/robustness/simulation_core.R'))
source(file.path(source_root,'scripts/robustness/partial_confounding_core.R'))
source(file.path(source_root,'scripts/robustness/partial_confounding_pilot_metrics.R'))
source(file.path(core_dir,'core.R'))
source(file.path(source_root,'package/assert_frozen_winn.R'))
frozen <- assert_frozen_winn(source_root)
stopifnot(as.character(packageVersion('winn'))=='0.1.6')
grid <- figure2_grid(); stopifnot(task>=1,task<=nrow(grid)*10)
seed_num <- (task-1L)%/%nrow(grid)+1L; g <- grid[(task-1L)%%nrow(grid)+1L,]
seed_id <- sprintf('CONF%02d',seed_num)
bundle <- readRDS(file.path(source_root,'data/simulated/partial_confounding/bundles',paste0(seed_id,'_none_0p00.rds')))
base <- generate_canonical_simulation(bundle$generation_seeds,include_artifact_matrices=TRUE)
truth <- bundle$evaluation_truth$feature_truth
priority <- partial_confounding_allocation_priority(base$metadata,bundle$component_seeds[['allocation_priority']])
a <- figure2_allocate(base$metadata,priority,bundle$master_seed+911L,g$axis,g$strength)
x <- figure2_construct(base,a,truth)
diag <- partial_confounding_diagnostics(a)
diag$summary$overlap_samples <- sum(vapply(split(a$phenotype_numeric,a$plate),function(v) {
  v <- v[!is.na(v)]; if(length(unique(v))==2) length(v) else 0L
},integer(1)))
task_out <- file.path(out,sprintf('task_%03d',task));dir.create(task_out,recursive=TRUE,showWarnings=FALSE)
set.seed(bundle$component_seeds[['method_WINN_fixed_no_QC']])
capture <- capture.output(timing <- system.time(fit <- winn::winn(x$raw,
  batch=base$metadata$plate,run_order=base$metadata$run_order,control_samples=NULL,
  parameters='fixed',ljung_box_fitdf=0L,return_details=TRUE)),type='output')
writeLines(capture,file.path(task_out,'correction.log'))
stopifnot(identical(dim(fit$data),dim(x$raw)),all(is.finite(fit$data)),all(fit$data>=0))
matrices <- list(Raw=x$raw,WiNN=fit$data,Truth=x$clean)
clean_marginal <- figure2_marginal(x$clean,a,truth)
clean_adjusted <- fit_partial_confounding_effects(x$clean,a,truth)$estimates
convert_adjusted <- function(e) data.frame(feature_id=e$feature_id,responsive=e$responsive,
  effect=e$recovered_effect_log1p,selected=e$significant_bh_0p05)
clean_adjusted2 <- convert_adjusted(clean_adjusted)
rows <- list(); effects <- list()
for(method in names(matrices)) {
  marginal <- figure2_marginal(matrices[[method]],a,truth)
  adjusted <- fit_partial_confounding_effects(matrices[[method]],a,truth)$estimates
  for(endpoint in c('marginal','adjusted')) {
    e <- if(endpoint=='marginal') marginal else convert_adjusted(adjusted)
    eclean <- if(endpoint=='marginal') clean_marginal else clean_adjusted2
    estimable <- endpoint=='marginal' || all(adjusted$estimable)
    r <- figure2_summarize_effect(e,eclean,estimable)
    r$endpoint <- endpoint;r$method <- method;r$estimable <- estimable
    rows[[paste(method,endpoint)]] <- r
  }
  marginal$method <- method;effects[[method]] <- marginal
}
metrics <- do.call(rbind,rows);metrics$seed_id <- seed_id;metrics$axis <- g$axis;metrics$strength <- g$strength
write.csv(metrics,file.path(task_out,'metrics.csv'),row.names=FALSE)
con <- gzfile(file.path(task_out,'feature_effects.csv.gz'),'wt');write.csv(do.call(rbind,effects),con,row.names=FALSE);close(con)
write.csv(a,file.path(task_out,'allocation.csv'),row.names=FALSE)
write.csv(diag$summary,file.path(task_out,'design_diagnostics.csv'),row.names=FALSE)
write.csv(diag$batch_counts,file.path(task_out,'batch_counts.csv'),row.names=FALSE)
saveRDS(list(raw=x$raw,clean=x$clean,corrected=fit$data),file.path(task_out,'matrices.rds'))
jsonlite::write_json(list(task=task,seed_id=seed_id,axis=g$axis,strength=g$strength,
  package=frozen,technical_sha256=digest::digest(x$technical,algo='sha256'),
  baseline_sha256=digest::digest(base$artifact_matrices$clean_log,algo='sha256'),
  feature_truth_sha256=digest::digest(truth,algo='sha256'),
  generation_seeds=as.list(bundle$generation_seeds),plate_seed=bundle$master_seed+911L,
  runtime_seconds=unname(timing[['elapsed']]),phenotype_given_to_correction=FALSE,
  protocol_sha256=digest::digest(file=file.path(core_dir,'PROTOCOL.md'),algo='sha256')),
  file.path(task_out,'complete.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(task_out,'sessionInfo.txt'))
cat('Completed',task,seed_id,g$axis,g$strength,'\n')
