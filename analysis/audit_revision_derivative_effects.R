#!/usr/bin/env Rscript
# Evaluate updated frozen derivative exports only; correction and controls are not refitted.
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
frozen<-normalizePath(args[1]);out<-args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Same-feature biological effects compared with raw; raw is a reference, not
# known biological truth. Feature-level summaries are descriptive, without
# pretending that features are independent biological replication units.
association_datasets <- c("mtbls79", "waveica", "sacurine")
assoc_paths <- file.path(frozen, association_datasets,
  c("association_features.csv", "association_features.csv", "sacurine_covariate_associations.csv"))
assoc <- do.call(rbind, lapply(seq_along(assoc_paths), function(i) {
  x <- read.csv(assoc_paths[i], stringsAsFactors = FALSE)
  endpoint <- if (association_datasets[i] == "sacurine") x$endpoint else
    if (association_datasets[i] == "mtbls79") "cow_vs_sheep" else "biological_group"
  effect <- if (association_datasets[i] == "sacurine") x$coefficient else x$effect
  data.frame(dataset_key = association_datasets[i], endpoint = endpoint,
    method_id = x$method_id, feature_id = x$feature_id, effect = effect,
    coefficient = effect, stringsAsFactors = FALSE)
}))
effect_rows <- list(); k <- 0L
groups <- unique(assoc[c("dataset_key", "endpoint")])
for (i in seq_len(nrow(groups))) {
  part <- assoc[assoc$dataset_key == groups$dataset_key[i] & assoc$endpoint == groups$endpoint[i], ]
  raw <- part[part$method_id == "raw", ]
  for (m in unique(part$method_id)) {
    corrected <- part[part$method_id == m, ]
    matched <- match(raw$feature_id, corrected$feature_id)
    reference <- ifelse(is.na(raw$effect), raw$coefficient, raw$effect)
    estimate <- ifelse(is.na(corrected$effect[matched]), corrected$coefficient[matched], corrected$effect[matched])
    ok <- is.finite(reference) & is.finite(estimate)
    if (sum(ok) < 3L) next
    a <- reference[ok]; b <- estimate[ok]
    k <- k + 1L
    effect_rows[[k]] <- data.frame(dataset_key = groups$dataset_key[i], endpoint = groups$endpoint[i],
      method_id = m, common_features = length(a), pearson_vs_raw = cor(a, b),
      spearman_vs_raw = cor(a, b, method = "spearman"),
      sign_concordance_vs_raw = mean(sign(a) == sign(b)),
      slope_through_origin_vs_raw = sum(a * b) / sum(a^2),
      median_absolute_difference = median(abs(b - a)),
      interpretation = "descriptive same-feature effects; raw is not biological truth")
  }
}
write.csv(do.call(rbind, effect_rows), file.path(out, "paired_feature_effects_descriptive.csv"), row.names = FALSE)
write.csv(data.frame(path = assoc_paths, md5 = tools::md5sum(assoc_paths)),
  file.path(out, "frozen_derivative_input_manifest.csv"), row.names = FALSE)

# Exact bootstrap source files used by the query-weighted retrieval update.
sources<-unlist(lapply(c("batchcorr_set1","clinical_fiams"),function(d)
 file.path(frozen,d,c("bootstrap_estimates.csv","bootstrap_replicates.csv","bootstrap_paired_differences.csv","matching_metrics_by_identity.csv"))))
stopifnot(all(file.exists(sources)))
write.csv(data.frame(path=sources,md5=tools::md5sum(sources)),file.path(out,"frozen_identity_retrieval_input_manifest.csv"),row.names=FALSE)
a<-read.csv(file.path(frozen,"clinical_fiams/bootstrap_estimates.csv"));a<-a[a$metric=="batch_weighted_pc_r2_categorical",]
b<-read.csv(file.path(frozen,"clinical_fiams/bootstrap_replicates.csv"));b<-b[b$metric=="batch_weighted_pc_r2_categorical",]
a$recomputed_lower<-vapply(a$method_id,function(id)quantile(b$value[b$method_id==id],.025,names=FALSE),numeric(1))
a$recomputed_upper<-vapply(a$method_id,function(id)quantile(b$value[b$method_id==id],.975,names=FALSE),numeric(1))
a$bootstrap_shift<-a$bootstrap_mean-a$point_estimate
a$point_below_interval<-a$point_estimate<a$ci_lower
stopifnot(max(abs(a$recomputed_lower-a$ci_lower))<1e-12,max(abs(a$recomputed_upper-a$ci_upper))<1e-12)
write.csv(a,file.path(out,"clinical_batch_r2_bootstrap_bias_audit.csv"),row.names=FALSE)
cat("Updated descriptive feature-effect and clinical bootstrap audits.\n")
