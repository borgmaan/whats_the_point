rm#!/usr/bin/env Rscript
# Andrew Borgman
# 2/16/2013
# Looking for extended regions of elevated Fst to 
# infer haplotypes over 
library(ggplot2)
library(reshape)

#setwd("/home/andrew/Dropbox/pointing/scripts/")
setwd("/home/andrew/projects/whats_the_point/visualization/")
source('fst_functions.R')

pp = df_from_dir("/media/3C52A21052A1CF48/pointer_vs_pointer/")
pn = df_from_dir("/media/3C52A21052A1CF48/pointer_vs_non/")

# Take differences and compute averages across all combinations
diffs = pn[, 4:length(pn)] - pp[, 4:length(pp)] 
all_avgs = apply(diffs, 1, function(x){return(mean(x, na.rm=T))})

dat = cbind(pp[,1:3], all_avgs)

# Plot this shit
names(dat) = c("SNP", "CHR", "BP", "d_i")

sig = dat[dat$d_i > quantile(dat$d_i, .999), ]
write.csv(sig, "/home/andrew/Dropbox/pointing/hits/adj_di_sig.csv", row.names=F)


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
#outfile_name = paste("/media/3C52A21052A1CF48/images_pointer_vs_pointer/", f,sep="")

png("/home/andrew/Dropbox/pointing/hits/d_i_no_labs.png", height=700, width=1800)

par(mar=c(5.1,6,4.3,2.1))
plot(d_i ~ BP, xlab='Chromosome',ylab=expression(d(italic(i))),main="Adjusted D(i) Values Averaged\nAcross All Pointing Breeds", axes=F,frame.plot=T, cex.lab = 1.8, cex.main=1.8, data=dat, col="#FFFFFF")
axis(1, at=ticks, cex.axis = 1.8, lab=c(seq(1,(length(unique(dat$CHR)) - 1)), "X"), lwd=4)
axis(2, cex.axis = 1.8, lwd=4)
box(lwd=4)
for (i in 1:39) {
  with(dat[dat$CHR==i, ],points(BP, d_i, col=cols[i],pch=20))  
}
abline(h=quantile(dat$d_i, .999), lwd=5, col='red')


dev.off()


