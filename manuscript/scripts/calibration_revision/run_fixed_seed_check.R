#!/usr/bin/env Rscript
# Independent check using the production default permutation seed in every fit.
a<-commandArgs(TRUE);stopifnot(length(a)==2L)
source(file.path(normalizePath(a[1]),'R','winn.R'))
rows<-vector('list',1000L)
for(s in seq_len(1000L)) {
 set.seed(8100000L+s)
 x<-adjust_outliers_mad(expm1(3+matrix(rnorm(100L*50L,sd=.15),100L)))
 p<-.compute_permutation_ljung_box(x,10L)$adjusted_p_value
 rows[[s]]<-data.frame(seed=s,n=50L,features=100L,permutations=999L,
   permutation_seed=104729L,rejected_features=sum(p<.05),any_rejection=any(p<.05))
}
d<-do.call(rbind,rows);write.csv(d,a[2],row.names=FALSE)
print(binom.test(sum(d$any_rejection),nrow(d),p=.05))
