# Confounding-only design. Source the unchanged canonical and partial helpers first.
figure2_grid <- function() {
  g <- expand.grid(axis=c('batch','run_order','joint'),
    strength=c(.25,.5,.75,.9,.95,1), stringsAsFactors=FALSE)
  rbind(data.frame(axis='none',strength=0),g)
}

figure2_allocate <- function(metadata, priority, plate_seed, axis, strength) {
  a <- partial_confounding_order_metadata(metadata)
  plates <- unique(a$plate)
  stopifnot(length(plates)==5, all(table(a$plate[a$sample_type=='sample'])==90))
  enriched <- .partial_preserve_rng({set.seed(plate_seed); sample(plates,2)})
  target <- ifelse(plates %in% enriched,90,0)
  ideal <- if(axis %in% c('batch','joint')) 36+(target-36)*strength else rep(36,5)
  counts <- floor(ideal+1e-10)
  left <- 180-sum(counts)
  if(left>0) {
    ix <- order(-(ideal-counts),seq_along(plates))[seq_len(left)]
    counts[ix] <- counts[ix]+1L
  }
  labels <- rep(NA_integer_,nrow(a))
  priority <- priority[a$sample_id]
  for(k in seq_along(plates)) {
    ix <- which(a$plate==plates[k] & a$sample_type=='sample')
    lab <- integer(length(ix))
    if(counts[k]>0) lab[order(priority[ix])[seq_len(counts[k])]] <- 1L
    if(axis %in% c('run_order','joint'))
      lab <- .partial_progress_to_early_cases(lab,a$study_order_rank[ix],strength)
    labels[ix] <- lab
  }
  a$phenotype_numeric <- labels
  a$phenotype <- ifelse(a$sample_type=='control','QC',ifelse(labels==1,'case','control'))
  a$confounding_axis <- axis; a$nominal_strength <- strength
  stopifnot(sum(labels==1,na.rm=TRUE)==180,sum(labels==0,na.rm=TRUE)==270)
  a
}

figure2_construct <- function(base,allocation,truth) {
  z <- base$artifact_matrices$clean_log
  beta <- truth$phenotype_effect_log[match(rownames(z),truth$metabolite)]
  code <- allocation$phenotype_numeric-.4
  code[allocation$sample_type=='control'] <- 0
  biological <- tcrossprod(beta,code)
  dimnames(biological) <- dimnames(z)
  technical <- base$artifact_matrices$batch_shift_by_injection + base$artifact_matrices$drift_by_injection
  clean <- z+biological
  raw <- clean+technical
  stopifnot(max(abs((raw-clean)-technical))<1e-12,
    all(biological[,allocation$sample_type=='control']==0),
    sum(truth$responsive)==100,all(abs(beta[truth$responsive])==.35))
  list(raw=exp(raw),clean=exp(clean),technical=technical,biological=biological)
}

figure2_marginal <- function(x,allocation,truth) {
  a <- allocation[allocation$sample_type=='sample',]
  z <- log1p(x[,a$sample_id,drop=FALSE]); y <- a$phenotype_numeric
  z1 <- z[,y==1,drop=FALSE];z0 <- z[,y==0,drop=FALSE]
  n1 <- ncol(z1);n0 <- ncol(z0)
  mu1 <- rowMeans(z1);mu0 <- rowMeans(z0)
  v1 <- rowSums((z1-mu1)^2)/(n1-1);v0 <- rowSums((z0-mu0)^2)/(n0-1)
  se2 <- v1/n1+v0/n0
  df <- se2^2/((v1/n1)^2/(n1-1)+(v0/n0)^2/(n0-1))
  statistic <- (mu1-mu0)/sqrt(se2)
  p <- 2*pt(-abs(statistic),df=df)
  p[se2==0 & mu1==mu0] <- 1
  p[se2==0 & mu1!=mu0] <- 0
  stopifnot(all(is.finite(p)))
  padj <- p.adjust(p,'BH')
  data.frame(feature_id=rownames(x),responsive=truth$responsive[match(rownames(x),truth$metabolite)],
    effect=mu1-mu0,se=sqrt(se2),df=df,p=p,padj=padj,selected=padj<.05)
}

figure2_summarize_effect <- function(e,clean,estimable=TRUE) {
  if(!estimable) return(data.frame(metric=c('effect_recovery','effect_rmse','sign_agreement',
    'true_positive_rate','empirical_fdr','discoveries','false_discoveries','zero_discoveries'),value=NA_real_))
  r <- e$responsive; sel <- e$selected; n <- sum(sel)
  stopifnot(all(is.finite(clean$effect[r])),all(abs(clean$effect[r])>1e-8))
  data.frame(metric=c('effect_recovery','effect_rmse','sign_agreement','true_positive_rate',
    'empirical_fdr','discoveries','false_discoveries','zero_discoveries'),
    value=c(median(e$effect[r]/clean$effect[r]),sqrt(mean((e$effect[r]-clean$effect[r])^2)),
      mean(sign(e$effect[r])==sign(clean$effect[r])),mean(sel[r]),sum(sel & !r)/max(1,n),
      n,sum(sel & !r),as.numeric(n==0)))
}
