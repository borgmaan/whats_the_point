#!/usr/bin/env Rscript
# Andrew Borgman
# 3/11/2013
# Targeted interrogation of XP-EHH & Fst calcs
setwd("/home/andrew/Dropbox/pointing/")
source("~/multiplot.R")


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

# Tenner
tenner = mm[mm$CHR == 10, ]
range(tenner$BP)

# Plotter
p1 = ggplot() + geom_line(data=tenner,aes(x=BP, y=XP_EHH), size=1, col="#117733") +  theme_bw(22) + ylab("XP-EHH") + xlab("Chr. 10 Base Position") + ggtitle("XP-EHH Values") +
  theme(panel.border = element_rect(colour = "black", size=2), axis.ticks = element_line(size=1.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_abline(slope=0, intercept=ehh_pct, col="red", size=1.5, alpha=0.8)


p2 = ggplot() + geom_line(data=tenner,aes(x=BP, y=D_i), size=1, col="#377EB9") +  theme_bw(22) + ylab("D(i)") + xlab("Chr. 10 Base Position") + ggtitle("Adjusted D(i) Values") +
  theme(panel.border = element_rect(colour = "black", size=2), axis.ticks = element_line(size=1.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_abline(slope=0, intercept=di_pct, col="red", size=1.5, alpha=0.8)

multiplot(p1,p2)


# Zooming even farther (CHR 10  +)
zoomer = tenner[tenner$BP < 20000000, ]

rect_left = seq(0, 18, 2) * 1000000
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

p1 = ggplot() + geom_rect(data=rectangles_ehh, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.5) +  
  geom_line(data=zoomer,aes(x=BP, y=XP_EHH), size=1.3, col="#117733") +  theme_bw(22) + ylab("XP-EHH") + xlab("Chr. 10 Base Position") +
  ggtitle("XP-EHH Values") + scale_x_continuous(breaks=seq(0,20000000,5000000), labels=c("0", "5,000,000", "10,000,000", "15,000,000", "20,000,000")) +
  theme(panel.border = element_rect(colour = "black", size=2), axis.ticks = element_line(size=1.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_abline(slope=0, intercept=ehh_pct, col="red", size=1.5, alpha=0.8)

p2 = ggplot()  + geom_rect(data=rectangles_di, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.5) +   
  geom_line(data=zoomer,aes(x=BP, y=D_i), size=1, col="#377EB9") +  theme_bw(22) + ylab("D(i)") + xlab("Chr. 10 Base Position") + 
  ggtitle("Adjusted D(i) Values") + scale_x_continuous(breaks=seq(0,20000000,5000000), labels=c("0", "5,000,000", "10,000,000", "15,000,000", "20,000,000")) +
  theme(panel.border = element_rect(colour = "black", size=2), axis.ticks = element_line(size=1.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_abline(slope=0, intercept=di_pct, col="red", size=1.5, alpha=0.8)

multiplot(p1,p2)
