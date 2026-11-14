#!/usr/bin/env Rscript
# Reproduce the optional identity-balanced sensitivity endpoint with its frozen
# 27 October 2026 seed and 5,000 paired draws. This is distinct from the main
# query-weighted 2,000-draw endpoint; no correction or selection is refitted.
a <- commandArgs(TRUE)
stopifnot(length(a) == 2L)
p <- file.path(a[1], "batchcorr_set1", "matching_metrics_by_identity.csv")
if (!file.exists(p)) p <- paste0(p, ".gz")
retrieval <- read.csv(p, stringsAsFactors = FALSE)
identities <- sort(unique(retrieval$identity))
methods <- unique(retrieval$method_id)
stopifnot(all(table(retrieval$method_id) == length(identities)))
lookup <- lapply(methods, function(m) {
  part <- retrieval[retrieval$method_id == m, ]
  part[match(identities, part$identity), ]
})
names(lookup) <- methods
set.seed(20261027)
draws <- replicate(5000, sample.int(length(identities), replace = TRUE))
rows <- list()
for (m in methods) for (endpoint in c("top1", "reciprocal_rank", "similarity_margin")) {
  value <- if (endpoint == "similarity_margin") lookup[[m]]$best_match_similarity - lookup[[m]]$best_unmatched_similarity else lookup[[m]][[endpoint]]
  raw <- if (endpoint == "similarity_margin") lookup$raw$best_match_similarity - lookup$raw$best_unmatched_similarity else lookup$raw[[endpoint]]
  boots <- colMeans(matrix(value[draws], nrow(draws)))
  deltas <- colMeans(matrix((value - raw)[draws], nrow(draws)))
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = "BatchCorr Set 1", method_id = m, method = lookup[[m]]$method[1],
    endpoint = endpoint, identities = length(identities), estimate = mean(value),
    lower = unname(quantile(boots, .025)), upper = unname(quantile(boots, .975)),
    paired_difference_vs_raw = mean(value - raw),
    delta_lower = unname(quantile(deltas, .025)), delta_upper = unname(quantile(deltas, .975)),
    uncertainty = "paired identity bootstrap conditional on frozen correction; 5000 draws"
  )
}
dir.create(a[2], recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, rows), file.path(a[2], "batchcorr_identity_retrieval_paired.csv"), row.names = FALSE)
cat("Reproduced identity-balanced sensitivity on", length(methods), "methods with 5000 paired draws.\n")
