#!/usr/bin/env Rscript
# Publication plots from complete validated summaries; no fitting or selection.
a <- commandArgs(TRUE); stopifnot(length(a)==2L)
d <- read.csv(a[1],stringsAsFactors=FALSE); out <- a[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source_out <- file.path(dirname(out),'source_data')
dir.create(source_out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(ggplot2))
labels <- c(asymptotic_BH='WiNN 0.1.6: BH',permutation_maxT='Alternative: maxT')
d$algorithm_label <- factor(labels[d$algorithm],levels=unname(labels))
script_arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
script_root <- dirname(dirname(normalizePath(sub('^--file=','',script_arg))))
source(file.path(script_root,'method_palette.R'))
colours <- setNames(rep(winn_colour,2),unname(labels))
algorithm_scales <- list(scale_colour_manual(values=colours),
  scale_shape_manual(values=c('WiNN 0.1.6: BH'=1,'Alternative: maxT'=16)),
  scale_linetype_manual(values=c('WiNN 0.1.6: BH'=2,'Alternative: maxT'=1)))
base <- theme_classic(base_size=10) + theme(legend.position='bottom',legend.title=element_blank(),
  plot.title=element_text(face='bold',size=10.5),
  axis.title=element_text(size=9.5),axis.text=element_text(colour='#263247',size=9),
  legend.text=element_text(size=9),strip.text=element_text(size=9),
  strip.background=element_rect(fill='#F0F2F5',colour=NA),plot.margin=margin(4,5,4,4))
save <- function(p,stem,data,w=174,h=58) {
  ggsave(file.path(out,paste0(stem,'.pdf')),p,width=w,height=h,units='mm',useDingbats=FALSE)
  write.csv(data,file.path(source_out,paste0(stem,'.csv')),row.names=FALSE)
}
x <- d[d$experiment=='original_null' & d$stage=='mad_prepared',]
p <- ggplot(x,aes(factor(n),100*fwer,colour=algorithm_label,shape=algorithm_label,linetype=algorithm_label,group=algorithm_label)) +
  geom_hline(yintercept=5,linetype=2,colour='#555555',linewidth=.4) +
  geom_line(linewidth=.5) + geom_errorbar(aes(ymin=100*fwer_lower,ymax=100*fwer_upper),width=.09,linewidth=.4) +
  geom_point(size=2) + algorithm_scales +
  labs(title='A  Calibration under independent null profiles',x='Injections per segment',
       y='Experiments with any\nfalse selection (%)') + base
save(p,'S11A_null_calibration',x)
x <- d[d$experiment=='robustness' & d$scenario=='global_null' & d$rho==.6,]
x$distribution <- factor(x$distribution,levels=c('gaussian','student_t3','skewed'),
  labels=c('Gaussian','Heavy-tailed (t3)','Skewed'))
p <- ggplot(x,aes(factor(n),100*fwer,colour=algorithm_label,shape=algorithm_label,linetype=algorithm_label,group=algorithm_label)) +
  geom_hline(yintercept=5,linetype=2,colour='#555555',linewidth=.4) +
  geom_line(linewidth=.5) + geom_point(size=1.7) +
  geom_errorbar(aes(ymin=100*fwer_lower,ymax=100*fwer_upper),width=.09,linewidth=.4) +
  facet_wrap(~distribution,nrow=1) + algorithm_scales +
  labs(title='B  Correlated null profiles (shared-component weight 0.6)',x='Injections per segment',
       y='Experiments with any\nfalse selection (%)') + base
save(p,'S11B_correlated_calibration',x)
x <- d[d$experiment=='robustness' & d$scenario=='mixed_linear_drift' & d$distribution=='gaussian',]
x$correlation <- factor(x$rho,levels=c(0,.6),labels=c('Independent features','Shared-component weight 0.6'))
p <- ggplot(x,aes(factor(n),100*sensitivity,colour=algorithm_label,shape=algorithm_label,linetype=algorithm_label,group=algorithm_label)) +
  geom_line(linewidth=.5) +
  geom_errorbar(aes(ymin=100*pmax(0,sensitivity-sensitivity_se),ymax=100*pmin(1,sensitivity+sensitivity_se)),width=.09,linewidth=.4) +
  geom_point(size=2) + facet_wrap(~correlation,nrow=1) +
  algorithm_scales + scale_y_continuous(limits=c(0,100)) +
  labs(title='C  Detection of known linear drift in 20 of 100 features',x='Injections per segment',
       y='Mean sensitivity (%)') + base
save(p,'S11C_drift_sensitivity',x)
