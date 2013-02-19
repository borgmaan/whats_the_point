#!/usr/bin/env Rscript
# Andrew Borgman
# 2/16/2013
# Analysis of variance simulation

setwd("/home/andrew/Dropbox/pointing/var_sims/")

# In case we need margins
defaults = par()$mar

dat = read.csv('sim_results.csv', stringsAsFactors=F, check.names=F)

#### Difference analysis

mean_cols = dat[, grep("*_mean", names(dat))]

# Lots of histograms
png("mean_differences.png", height=600, width=600)
par(mar=defaults)
par(mfrow=c(5, 3))
for (i in 1:length(mean_cols)){
  samp_size = strsplit(names(mean_cols)[i], "_")[[1]][1]
  hist(mean_cols[,i], breaks=50, xlim = c(-0.1, 0.1), main=paste("Sample Size of", samp_size), xlab="Fst Deviance", col='grey', cex.lab=1.3, cex.main=1.3, cex.axis=1.3)
}
dev.off()

# Average differences
par(mfrow=c(1, 1))
averages = apply(mean_cols, 2, function(x){return(mean(x, na.rm=T))})
names(averages) = paste("Samples of", seq(6,90,6))
new_mars = defaults
new_mars[1] = 6.5
new_mars[2] = 6.5
par(mar=new_mars)
barplot(-averages, las=2, ylab="Avg. Difference From True", main="Deviance From True Fst\nfor Random Samples")

# SD analysis
sd_cols = dat[, grep("*_sd", names(dat))]

# Histocrams
png("sd_differences.png", height=600, width=600)
par(mar=defaults)
par(mfrow=c(5, 3))
for (i in 1:length(sd_cols)){
  samp_size = strsplit(names(sd_cols)[i], "_")[[1]][1]
  hist(sd_cols[,i], breaks=50, xlim = c(0, 0.3), main=paste("Sample Size of", samp_size), xlab="STEV of Fst Deviance", col='grey', cex.lab=1.3, cex.main=1.3, cex.axis=1.3)
}
dev.off()

# Far clots
par(mfrow=c(1, 1))
averages = apply(sd_cols, 2, function(x){return(mean(x, na.rm=T))})
names(averages) = paste("Samples of", seq(6,90,6))
new_mars = defaults
new_mars[1] = 6.5
new_mars[2] = 6.5
par(mar=new_mars)
barplot(averages, las=2, ylab="Avg. STDEV from Sampling", main="STDEV of Fst Estimates\nfor Random Samples")



