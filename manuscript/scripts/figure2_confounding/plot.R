# Only Figure 2E-F are written. All numerical inputs must pass summarize.py first.
plot_figure2_confounding <- function(root) {
  suppressPackageStartupMessages(library(ggplot2))
  x <- read.csv(file.path(root,'generated/figure2_confounding/summary.csv'))
  x <- x[x$method=='WiNN' & x$endpoint=='marginal',]
  x$axis <- factor(x$axis,levels=c('batch','run_order','joint'),labels=c('Batch','Run order','Both'))
  x$alignment_percent <- 100*x$strength
  x$lower <- x$mean-x$se;x$upper <- x$mean+x$se
  stopifnot(all(x$n_seeds==10),all(is.finite(x$mean)),all(is.finite(x$se)))
  base <- theme_classic(base_size=9.5,base_family='sans')+
    theme(axis.line=element_line(colour='#596273',linewidth=.16),
      axis.ticks=element_line(colour='#596273',linewidth=.16),
      panel.grid.major=element_line(colour='#E2E7ED',linewidth=.16),
      plot.title=element_text(colour='#172033',size=9.5),plot.title.position='plot',
      axis.title=element_text(size=9.3),axis.text=element_text(size=9),
      legend.text=element_text(size=9),legend.title=element_blank(),
      legend.position='bottom',legend.key.width=grid::unit(5,'mm'),
      legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0),
      plot.margin=margin(3,4,3,4),strip.background=element_rect(fill='#F0F2F5',colour=NA),
      strip.text=element_text(size=9),panel.spacing.x=grid::unit(3,'mm'))
  points <- function(d,x_expand=c(.02,.04)) ggplot(d,aes(alignment_percent,mean,shape=axis,linetype=axis,group=axis))+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=1.5,linewidth=.16,colour='#D7191C',linetype='solid')+
    geom_line(linewidth=.22,colour='#D7191C')+geom_point(size=1.5,colour='#D7191C')+
    scale_shape_manual(values=c(Batch=16,'Run order'=17,Both=15))+
    scale_linetype_manual(values=c(Batch='solid','Run order'='dashed',Both='dotdash'))+
    scale_x_continuous(breaks=c(0,50,100),expand=expansion(mult=x_expand))+
    labs(x='Alignment strength (%)',shape=NULL,linetype=NULL)+base
  e <- x[x$metric=='effect_recovery',]
  pe <- points(e)+scale_y_continuous(limits=c(min(0,min(e$lower)),max(1,max(e$upper))),breaks=c(0,.25,.5,.75,1),expand=expansion(mult=c(.02,.03)))+labs(title='E WiNN: retained group contrast',y='Effect recovery')+
    geom_hline(yintercept=1,colour='grey65',linewidth=.16,linetype='dotted')
  f <- x[x$metric %in% c('true_positive_rate','empirical_fdr'),]
  f$endpoint_label <- factor(f$metric,levels=c('true_positive_rate','empirical_fdr'),labels=c('Sensitivity','FDR'))
  pf <- points(f,c(.05,.15))+facet_wrap(~endpoint_label,nrow=1,scales="free_y")+labs(title='F WiNN: detection and false discoveries',y='Rate')+
    scale_y_continuous(expand=expansion(mult=c(.04,.07)))+
    geom_blank(data=data.frame(endpoint_label=factor(c('Sensitivity','Sensitivity','FDR'),levels=c('Sensitivity','FDR')),y=c(0,1,0)),aes(y=y),inherit.aes=FALSE)+
    geom_hline(data=data.frame(endpoint_label=factor('FDR',levels=c('Sensitivity','FDR')),y=.05),
      aes(yintercept=y),colour='grey55',linewidth=.16)
  specs <- list(list('E','F2E_median_attenuation_ratio_clean_responsive',pe,e),
    list('F','F2F_empirical_fdr',pf,f))
  for(s in specs) {
    ggsave(file.path(root,'figures/main_panels',paste0(s[[2]],'.pdf')),s[[3]],
      width=85,height=50,units='mm',device=grDevices::pdf,useDingbats=FALSE,bg='white')
    write.csv(s[[4]],file.path(root,'figures/source_data',paste0(s[[2]],'.csv')),row.names=FALSE)
  }
  invisible(specs)
}
if(sys.nframe()==0L) {
  arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
  root <- dirname(dirname(dirname(normalizePath(sub('^--file=','',arg)))))
  plot_figure2_confounding(root)
}
