#!/usr/bin/env Rscript
# Andrew Borgman
# 2/15/2013
# Comparing d(i) values from within and between runs
library(ggplot2)

pn_files = paste("/media/3C52A21052A1CF48/pointer_vs_non/", dir("/media/3C52A21052A1CF48/pointer_vs_non/"), sep="")
pp_files = paste("/media/3C52A21052A1CF48/pointer_vs_pointer/", dir("/media/3C52A21052A1CF48/pointer_vs_pointer/"), sep="")
files = data.frame(cbind(pn_files,pp_files))
files[,1] = as.character(files[,1])
files[,2] = as.character(files[,2])

for (i in 1:length(files[,1])){

  # Grab within and against data for a breed
  pn = read.table(files[i,1], stringsAsFactors=F)
  pp = read.table(files[i,2], stringsAsFactors=F)
  
  # Grab breed name
  f = strsplit(files[i,1], split="/")[[1]][length(strsplit(files[i,1], split="/")[[1]])]
  
  # Subtract away pointer/pointer d(i) from pointer/non d(i)
  dat = pn[,1:3]
  dat$diff = pn[,4] - pp[,4]
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
  outfile_name = paste("/media/3C52A21052A1CF48/diff_images/", f,sep="")
  
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