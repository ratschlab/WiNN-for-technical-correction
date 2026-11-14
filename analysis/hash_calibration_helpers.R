#!/usr/bin/env Rscript
# Compare algorithm bodies across validation and installed/package-source snapshots.
a <- commandArgs(TRUE)
stopifnot(length(a)==2L)
env <- new.env(parent=globalenv());source(a[1],local=env)
functions <- c('.compute_permutation_ljung_box','.ljung_box_statistic_centered')
rows <- lapply(functions,function(name) {
  fun <- get(name,envir=env,inherits=FALSE)
  code <- paste(c(deparse(formals(fun),width.cutoff=500L),deparse(body(fun),width.cutoff=500L)),collapse='\n')
  data.frame(function_name=name,sha256=digest::digest(code,algo='sha256',serialize=FALSE))
})
write.csv(do.call(rbind,rows),a[2],row.names=FALSE)
