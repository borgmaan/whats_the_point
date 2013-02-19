#!/usr/bin/env Rscript
# Andrew Borgman
# 2/13/2013
# Exploring results of pairwise Fst analysis
library(ggplot2)

#setwd("/media/3C52A21052A1CF48/pointer_vs_non/")
setwd("/media/3C52A21052A1CF48/pointer_vs_pointer/")

files = dir()
for (f in files){
  dat = read.table(f, stringsAsFactors=F)
  names(dat) = c("SNP", "CHR", "BP", "d_i")
  
  #Plotting!
  dat$BP <- dat$BP / 100
  lastbase <- 0
  newbase <- 0
  ticks <- NULL
  
  for(i in 1:length(unique(dat$CHR))){
    cat(i)
    if(i == 1){
      newbase <- range(dat[dat$CHR == i,]$BP)[2]
      ticks <- c(ticks, (newbase - lastbase) / 2)
      lastbase <- newbase
    }
    else{
      dat[dat$CHR == i,]$BP <- dat[dat$CHR == i,]$BP + lastbase
      newbase <- range(dat[dat$CHR == i,]$BP)[2]
      ticks <- c(ticks, ((newbase - lastbase) / 2)+lastbase)
      lastbase <- newbase
    }
  }
  
  #Specifying colors
  cols = rep(c('#000000','#585858'),40)
  
  #Plotting points
  #outfile_name = paste("/media/3C52A21052A1CF48/images_pointer_vs_non/", f,sep="")
  outfile_name = paste("/media/3C52A21052A1CF48/images_pointer_vs_pointer/", f,sep="")
  
  png(outfile_name, height=600, width=1800)
  par(mar=c(5.1,6,4.3,2.1))
  plot(d_i ~ BP, xlab='Chromosome',ylab=expression(d(italic(i))),main=f, axes=F,frame.plot=T, cex.lab = 1.8, cex.main=1.8, data=dat, col="#FFFFFF")
  axis(1, at=ticks, cex.axis = 1.8, lab=seq(1,length(unique(dat$CHR))))
  axis(2, cex.axis = 1.8)
  for (i in 1:39) {
    with(dat[dat$CHR==i, ],points(BP, d_i, col=cols[i],pch=20))  
  }
  dev.off()
}

