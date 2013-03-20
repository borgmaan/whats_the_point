#!/usr/bin/env Rscript
# Andrew Borgman
# 3/11/2013
# Targeted interrogation of XP-EHH & Fst calcs
library(ggplot2)
library(stringr)
source("~/multiplot.R")

# Daterrrr
setwd("/home/andrew/Dropbox/pointing/")
ehh_vals = read.table("/home/andrew/projects/whats_the_point/data/xpehh_results.tsv", header=T, sep="\t", stringsAsFactors=F)
di_vals = read.csv("/home/andrew/projects/whats_the_point/data/adjusted_di_avgs.csv", stringsAsFactors=F)

# Merge'em
mm = merge(ehh_vals, di_vals, by.x=1, by.y=1)

# Renamer
mm = mm[, -c(7,8)]
names(mm) = c("SNP", "CHR", "BP", "iHs_Pop_1", "iHs_Pop_2", "XP_EHH", "D_i")

# Grab some percentile lines
di_pct = quantile(mm$D_i, .999)
ehh_pct = quantile(mm$XP_EHH, .999)

# Region Plotting Function
region_plotter <- function(chrom, start_pos, stop_pos){

  # Grab subsetted DF for chromosome of interest
  tmp = mm[mm$CHR == chrom, ]
  
  # Zooming to region
  zoomer = tmp[tmp$BP > start_pos, ]
  zoomer = zoomer[zoomer$BP < stop_pos, ]
  ranger = range(zoomer$BP)
  
  rect_left = seq(start_pos, stop_pos, 2000000)
  rectangles_ehh <- data.frame(
    xmin = rect_left,
    xmax = rect_left + 1000000,
    ymin =  min(zoomer$XP_EHH),
    ymax = max(zoomer$XP_EHH) + 0.2
  )
  
  rect_left = seq(0, 18, 2) * 1000000
  rectangles_di <- data.frame(
    xmin = rect_left,
    xmax = rect_left + 1000000,
    ymin =  min(zoomer$D_i),
    ymax = max(zoomer$D_i) + 15
  )
  
  + scale_x_continuous(breaks=seq(0,20000000,5000000), labels=c("0", "5,000,000", "10,000,000", "15,000,000", "20,000,000"))
  
  p1 = ggplot() + geom_rect(data=rectangles_ehh, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.5) +  
    geom_line(data=zoomer,aes(x=BP, y=XP_EHH), size=1.3, col="#117733") +  theme_bw(22) + ylab("XP-EHH") + xlab(str_c("Chr. ", chrom, " Base Position")) +
    ggtitle("Cross-Population Extended Haplotype Homozygosity")  +  theme(panel.border = element_rect(colour = "black", size=2), axis.ticks = element_line(size=1.25),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_abline(slope=0, intercept=ehh_pct, col="red", size=1.5, alpha=0.8)
  
  p2 = ggplot()  + geom_rect(data=rectangles_di, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.5) +   
    geom_line(data=zoomer,aes(x=BP, y=D_i), size=1, col="#377EB9") +  theme_bw(22) + ylab("D(i)") + xlab("Chr. 10 Base Position") + 
    ggtitle(expression(Adjusted~D[(i)]~Statistics)) + scale_x_continuous(breaks=seq(0,20000000,5000000), labels=c("0", "5,000,000", "10,000,000", "15,000,000", "20,000,000")) +
    theme(panel.border = element_rect(colour = "black", size=2), axis.ticks = element_line(size=1.25),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_abline(slope=0, intercept=di_pct, col="red", size=1.5, alpha=0.8)
  
  multiplot(p1,p2)
  
  
}

# Tweaking for pub_quality
# Region 10
chrom = 10
start_pos = 0
stop_pos = 21000000

# Region 11
chrom = 11
start_pos = 25000000
stop_pos = 36000000


# Grab subsetted DF for chromosome of interest
tmp = mm[mm$CHR == chrom, ]

# Zooming to region
zoomer = tmp[tmp$BP > start_pos, ]
zoomer = zoomer[zoomer$BP < stop_pos, ]
ranger = range(zoomer$BP)

rect_left = seq(start_pos, stop_pos - 1000000, 2000000)
rectangles_ehh <- data.frame(
  xmin = rect_left,
  xmax = rect_left + 1000000,
  ymin =  min(zoomer$XP_EHH),
  ymax = max(zoomer$XP_EHH) + 0.2
)

rect_left = seq(start_pos, stop_pos - 1000000, 2000000)
rectangles_di <- data.frame(
  xmin = rect_left,
  xmax = rect_left + 1000000,
  ymin =  min(zoomer$D_i),
  ymax = max(zoomer$D_i) + 15
)


p1 = ggplot() + geom_rect(data=rectangles_ehh, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.5) +  
  geom_line(data=zoomer,aes(x=BP, y=XP_EHH), size=1.3, col="#117733") +  theme_bw(22) + ylab("XP-EHH") + xlab(str_c("Chr. ", chrom, " Base Position")) +
  ggtitle("Cross-Population Extended Haplotype Homozygosity") + scale_x_continuous(breaks=seq(25000000, 35000000, 2500000), labels=format(seq(25000000, 35000000, 2500000),scientific=F) ) + 
  theme(panel.border = element_rect(colour = "black", size=2), axis.ticks = element_line(size=1.25),
                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_abline(slope=0, intercept=ehh_pct, col="red", size=1.5, alpha=0.8)

p2 = ggplot()  + geom_rect(data=rectangles_di, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.5) +   
  geom_line(data=zoomer,aes(x=BP, y=D_i), size=1, col="#377EB9") +  theme_bw(22) + ylab("D(i)") + xlab("Chr. 10 Base Position") + 
  ggtitle(expression(Adjusted~D[(i)]~Statistics)) + scale_x_continuous(breaks=seq(25000000, 35000000, 2500000), labels=format(seq(25000000, 35000000, 2500000),scientific=F) ) + 
  theme(panel.border = element_rect(colour = "black", size=2), axis.ticks = element_line(size=1.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_abline(slope=0, intercept=di_pct, col="red", size=1.5, alpha=0.8)

multiplot(p1,p2)




# All Chromosome 11
chrom = 11
start_pos = 0
stop_pos = 75000000

# Grab subsetted DF for chromosome of interest
tmp = mm[mm$CHR == chrom, ]

# Zooming to region
zoomer = tmp[tmp$BP > start_pos, ]
zoomer = zoomer[zoomer$BP < stop_pos, ]
ranger = range(zoomer$BP)

rect_left = seq(start_pos, stop_pos - 5000000, 10000000)
rectangles_ehh <- data.frame(
  xmin = rect_left,
  xmax = rect_left + 5000000,
  ymin =  min(zoomer$XP_EHH),
  ymax = max(zoomer$XP_EHH) + 0.2
)

rect_left = seq(start_pos, stop_pos - 5000000, 10000000)
rectangles_di <- data.frame(
  xmin = rect_left,
  xmax = rect_left + 5000000,
  ymin =  min(zoomer$D_i),
  ymax = max(zoomer$D_i) + 15
)


p1 = ggplot() + geom_rect(data=rectangles_ehh, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.5) +  
  geom_line(data=zoomer,aes(x=BP, y=XP_EHH), size=1.3, col="#117733") +  theme_bw(25) + ylab("XP-EHH") + xlab(str_c("Chr. ", chrom, " Base Position")) +
  ggtitle("Cross-Population Extended Haplotype Homozygosity") + scale_x_continuous(breaks=seq(start_pos, stop_pos, 10000000), labels=format(seq(start_pos, stop_pos, 10000000),scientific=F) ) + 
  theme(panel.border = element_rect(colour = "black", size=2), axis.ticks = element_line(size=1.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_abline(slope=0, intercept=ehh_pct, col="red", size=1.5, alpha=0.8)

p2 = ggplot()  + geom_rect(data=rectangles_di, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.5) +   
  geom_line(data=zoomer,aes(x=BP, y=D_i), size=1, col="#377EB9") +  theme_bw(25) + ylab(expression(D[(i)])) + xlab(str_c("Chr. ", chrom, " Base Position")) + 
  ggtitle(expression(Adjusted~D[(i)]~Statistics)) + scale_x_continuous(breaks=seq(start_pos, stop_pos, 10000000), labels=format(seq(start_pos, stop_pos, 10000000),scientific=F) ) + 
  theme(panel.border = element_rect(colour = "black", size=2), axis.ticks = element_line(size=1.25), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope=0, intercept=di_pct, col="red", size=1.5, alpha=0.8)

multiplot(p1,p2)

