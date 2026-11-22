# Paired gate-policy display; only existing simulation draws are used.
plot_paired_gate <- function(root) {
  source(file.path(root, 'scripts/method_palette.R'), local=TRUE)
  x <- read.csv(file.path(root, 'generated/revision_panel_inputs/selective_forced_truth_rmse.csv'))
  expected_seeds <- sprintf('SIM%02d', seq_len(10L))
  expected_policies <- c('G_SS','G_AS','G_SA','G_AA')
  stopifnot(nrow(x)==40L, !anyDuplicated(x[c('seed_id','method_id')]),
    setequal(unique(x$method_id),expected_policies),
    all(x$metric=='truth_log1p_rmse'), all(is.finite(x$value)), all(x$value>=0),
    all(vapply(expected_policies,function(policy) {
      identical(sort(x$seed_id[x$method_id==policy]),expected_seeds)
    },logical(1))))
  ref <- x[x$method_id=='G_SS', c('seed_id','value')]
  names(ref)[2] <- 'selective_rmse'
  z <- merge(x[x$method_id!='G_SS',],ref,by='seed_id')
  stopifnot(nrow(z)==30L, all(is.finite(z$selective_rmse)),
    all(table(z$method_id)==10L))
  z$paired_difference <- z$value-z$selective_rmse
  z$policy <- factor(z$method_id,levels=c('G_AA','G_SA','G_AS'),labels=c('Force both','Force batch','Force drift'))
  summary <- do.call(rbind,lapply(split(z,z$method_id),function(d) {
    n<-nrow(d); m<-mean(d$paired_difference); se<-sd(d$paired_difference)/sqrt(n)
    data.frame(method_id=d$method_id[1],policy=d$policy[1],n=n,mean_difference=m,
      se=se,ci_lower=m-qt(.975,n-1)*se,ci_upper=m+qt(.975,n-1)*se,
      all_positive=all(d$paired_difference>0),relative_improvement_percent=100*m/mean(d$value))
  }))
  summary$policy<-factor(summary$policy,levels=levels(z$policy))
  suppressPackageStartupMessages(library(ggplot2))
  p<-ggplot(z,aes(1000*paired_difference,policy))+
    geom_vline(xintercept=0,linetype='dashed',colour='#657080',linewidth=.3)+
    geom_point(position=position_jitter(width=0,height=.08,seed=2026),colour=winn_colour,size=1.3,alpha=.7)+
    geom_segment(data=summary,aes(x=1000*ci_lower,xend=1000*ci_upper,y=policy,yend=policy),inherit.aes=FALSE,linewidth=.45,colour='#172033')+
    geom_point(data=summary,aes(x=1000*mean_difference,y=policy),inherit.aes=FALSE,shape=23,fill=winn_colour,colour='#172033',size=2.3,stroke=.3)+
    scale_x_continuous(expand=expansion(mult=c(.06,.12)))+
    labs(title='A Paired benefit of selection',x='RMSE difference: forced - selective (x 1000)',y=NULL)+
    theme_classic(base_size=9.5,base_family='sans')+
    theme(plot.title=element_text(face='bold',size=9.5,colour='#172033'),axis.title=element_text(size=9),
      axis.text=element_text(size=9,colour='#263247'),plot.margin=margin(3,5,3,4))
  list(plot=p,data=z,summary=summary)
}
if(sys.nframe()==0L) {
  a<-commandArgs(FALSE);path<-normalizePath(sub('^--file=','',grep('^--file=',a,value=TRUE)))
  root<-dirname(dirname(dirname(path))); result<-plot_paired_gate(root)
  out<-file.path(root,'generated/final_review_gating_matching');dir.create(out,recursive=TRUE,showWarnings=FALSE)
  write.csv(result$data,file.path(out,'paired_gate_rmse.csv'),row.names=FALSE)
  write.csv(result$summary,file.path(out,'paired_gate_rmse_summary.csv'),row.names=FALSE)
  write.csv(result$data,file.path(root,'figures/source_data/F2A_selective_forced_truth_rmse.csv'),row.names=FALSE)
  ggplot2::ggsave(file.path(root,'figures/main_panels/F2A_selective_forced_truth_rmse.pdf'),result$plot,width=85,height=50,units='mm',device=grDevices::pdf,useDingbats=FALSE)
  print(result$summary)
}
