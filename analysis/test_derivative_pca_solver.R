#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
source(file.path(repo_root, "analysis", "derivative_biology_helpers.R"), local = TRUE)

if (!requireNamespace("irlba", quietly = TRUE)) {
  stop("The derivative PCA equivalence test requires irlba.", call. = FALSE)
}

set.seed(20260726)
test_matrix <- matrix(stats::rnorm(380L * 340L), nrow = 380L, ncol = 340L)
test_matrix[1L, ] <- 1
target <- factor(rep(letters[1:4], each = 85L))

exact <- pca_basis_from_log_exact(test_matrix)
truncated <- pca_basis_from_log(test_matrix)
exact_r2 <- pca_r2_from_basis(exact, target, "categorical")
truncated_r2 <- pca_r2_from_basis(truncated, target, "categorical")

if (!isTRUE(all.equal(exact_r2, truncated_r2, tolerance = 1e-10))) {
  stop("Truncated and exact weighted-PC R-squared do not agree.", call. = FALSE)
}

permutation <- sample(seq_len(ncol(test_matrix)))
permuted_r2 <- pca_r2_from_log(test_matrix[, permutation, drop = FALSE], target[permutation], "categorical")
if (!isTRUE(all.equal(truncated_r2, permuted_r2, tolerance = 1e-10))) {
  stop("Weighted-PC R-squared is not invariant to a joint sample permutation.", call. = FALSE)
}

cat("Derivative PCA solver equivalence tests passed.\n")
