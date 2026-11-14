#!/usr/bin/env Rscript
source('scripts/robustness/simulation_core.R')
source('scripts/robustness/partial_confounding_core.R')
source('scripts/robustness/partial_confounding_pilot_metrics.R')
source('analysis/figure2_confounding/core.R')
g <- figure2_grid(); stopifnot(nrow(g)==19, sum(g$strength==0)==1)
checks <- 0L
for(seed in 1:10) {
  b <- readRDS(sprintf('data/simulated/partial_confounding/bundles/CONF%02d_none_0p00.rds',seed))
  base <- generate_canonical_simulation(b$generation_seeds,TRUE)
  p <- partial_confounding_allocation_priority(base$metadata,b$component_seeds[['allocation_priority']])
  allocations <- vector('list',nrow(g))
  for(i in seq_len(nrow(g))) {
    a <- figure2_allocate(base$metadata,p,b$master_seed+911L,g$axis[i],g$strength[i])
    stopifnot(sum(a$phenotype_numeric==1,na.rm=TRUE)==180,
      sum(a$phenotype_numeric==0,na.rm=TRUE)==270,
      all(is.na(a$phenotype_numeric[a$sample_type=='control'])))
    d <- partial_confounding_diagnostics(a)$summary
    if(g$strength[i]==1 && g$axis[i] %in% c('batch','joint'))
      stopifnot(!d$design_full_rank, abs(d$cramers_v_phenotype_batch-1)<1e-12)
    else stopifnot(d$design_full_rank)
    if(g$axis[i]=='run_order') {
      stopifnot(all(tapply(a$phenotype_numeric,a$plate,sum,na.rm=TRUE)==36))
      if(g$strength[i]==1) for(v in split(a[a$sample_type=='sample',],a$plate[a$sample_type=='sample']))
        stopifnot(max(v$study_order_rank[v$phenotype_numeric==1])<min(v$study_order_rank[v$phenotype_numeric==0]))
    }
    allocations[[i]] <- a;checks <- checks+1L
  }
  batch <- allocations[[which(g$axis=='batch' & g$strength==1)]]
  joint <- allocations[[which(g$axis=='joint' & g$strength==1)]]
  stopifnot(identical(batch$phenotype_numeric,joint$phenotype_numeric))
  x <- figure2_construct(base,batch,b$evaluation_truth$feature_truth)
  stopifnot(max(abs(rowSums(x$biological)))<1e-10)
  e <- figure2_marginal(x$clean,batch,b$evaluation_truth$feature_truth)
  y <- batch$phenotype_numeric
  for(j in 1:10) {
    z <- log1p(x$clean[j,]); ref <- t.test(z[which(y==1)],z[which(y==0)])
    stopifnot(abs(e$p[j]-ref$p.value)<1e-12,abs(e$effect[j]-diff(rev(ref$estimate)))<1e-12)
  }
}
cat('PASS:',checks,'allocation scenarios; full-confounding rank boundary; unchanged artifacts; Welch reference parity\n')
