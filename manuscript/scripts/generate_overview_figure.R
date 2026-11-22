#!/usr/bin/env Rscript
# Reproducible schematic; traces illustrate four technical-effect patterns,
# not observed data or fitted benchmark results.
args_full <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub('^--file=', '', grep('^--file=', args_full, value=TRUE)))
root <- dirname(dirname(script))
pdf(file.path(root, 'figures/F1_winn_overview.pdf'), width=175/25.4,
    height=53/25.4, pointsize=8, family='Helvetica', useDingbats=FALSE)
par(mar=c(0,0,0,0), xaxs='i', yaxs='i')
plot.new(); plot.window(xlim=c(0,175), ylim=c(0,65))
ink <- '#263247'; gray <- '#A5ADB8'; blue <- '#3568AE'; red <- '#BA5255'
rect(31,1,175,64,col='#F6F7F9',border=NA)
text(103,61.5,'WiNN workflow',font=2,cex=1.12,col=ink)
segments(32,58.8,174,58.8,col='#BBC2CC',lwd=.5)
text(15.5,55,'Randomized\ndesign',font=2,cex=1,col=ink)
text(45,55,'Metabolite\nprofiles',font=2,cex=1,col=ink)
text(64.5,55,'Run-order\ntest',cex=.94,col=ink)
text(84,55,'GAM\ndetrending',font=2,cex=1,col=ink)
text(104,55,'Batch\ntest',cex=.94,col=ink)
text(124,55,'Batch\ncorrection',font=2,cex=1,col=ink)
text(158,55,'Corrected\nprofiles',font=2,cex=1,col=ink)
text(15.5,43,'Balanced groups\nacross batches',cex=.94,col=ink)
pal <- c('#7992BE','#8CB69B','#BE8F9E','#63666A')
set.seed(20261122)
for (b in 1:2) {
  x0 <- if(b==1) 2 else 16.4
  rect(x0,21,x0+12.3,34,col='white',border=gray,lwd=.6)
  arrows(x0+1,32.5,x0+11.4,32.5,length=.035,col=ink,lwd=.6)
  groups <- sample(rep(1:4, each=5))
  for(k in 1:20) symbols(x0+1.7+((k-1)%%4)*2.95,
      30.4-floor((k-1)/4)*1.95,circles=.65,inches=FALSE,add=TRUE,
      bg=pal[groups[k]],fg='white',lwd=.4)
}
text(15.5,17,'Randomized within\neach batch',cex=.94,col=ink)
for(k in 1:3) symbols(4.4+(k-1)*2.7,8,circles=.8,inches=FALSE,add=TRUE,bg=pal[k],fg=NA)
text(12,8,'Groups',adj=c(0,.5),cex=.88,col=ink)
symbols(4.4,4.4,circles=.8,inches=FALSE,add=TRUE,bg=pal[4],fg=NA)
text(6.4,4.4,'QC',adj=c(0,.5),cex=.88,col=ink)
yc <- c(43,32,21,10)
t <- seq(0,1,length.out=100); seg <- pmin(4, floor(t*4)+1)
noise <- rnorm(length(t),sd=.13)
batch <- c(-.35,.45,-.65,.65)[seg]
drift <- .5*sin(8*pi*t)+.14*cos(20*pi*t)
raw <- list(noise, noise+drift, noise+batch,noise+drift+batch)
profile <- function(x0,y,v,fit=NULL,fitcol=blue) {
  rect(x0,y-3.6,x0+26,y+3.6,border='#89929E',col='white',lwd=.45)
  for(b in 1:3) segments(x0+26*b/4,y-3.6,x0+26*b/4,y+3.6,col='#C4CAD2',lty=2,lwd=.45)
  xx <- x0+.7+24.6*t
  lines(xx,y+2.2*v,col=ink,lwd=.5)
  if(!is.null(fit)) lines(xx,y+2.2*fit,col=fitcol,lwd=1.0)
}
mark <- function(x,y,pass) {
  if(pass) {lines(c(x-1.1,x-.35,x+1.25), c(y,y-.8,y+1),col='#447A58',lwd=1.1)}
  else {segments(x-.8,y-.8,x+.8,y+.8,col='#B45C62',lwd=1.0);segments(x-.8,y+.8,x+.8,y-.8,col='#B45C62',lwd=1.0)}
}
for(i in 1:4) {
  profile(32,yc[i],raw[[i]])
  mark(64.5,yc[i],i%%2==0)
  if(i%%2==0) profile(71,yc[i],raw[[i]],drift+if(i==4)batch else 0)
  mark(104,yc[i],i>=3)
  if(i>=3) profile(111,yc[i],noise+batch,batch,red)
  profile(145,yc[i],noise)
}
text(45,48.3,'Batches 1-4',cex=.88,col=ink)
# Stage arrows connect the process without claiming every row is corrected.
for(x in c(29,68.5,108,140.5)) arrows(x,26,x+2,26,length=.035,col='#6581A5',lwd=.7)
arrows(32,4.8,58,4.8,length=.035,col=ink,lwd=.5)
text(45,2.5,'Run order',cex=.88,col=ink)
dev.off()
