#!/usr/bin/env Rscript
# Andrew Borgman
# Making figures for Mark's talk

library(ggplot2)

#setwd("/home/andrew/Dropbox/pointing/rip/")
setwd("/Users/andrewborgman/Dropbox/pointing/rip/")

#call_rates = read.csv("/home/andrew/Dropbox/pointing/qc_data/breed_call_rates.csv")
call_rates = read.csv("/Users/andrewborgman/Dropbox/pointing/qc_data/breed_call_rates.csv", check.names=F))
call_rates$Mean_Call = call_rates$Mean_Call * 100
call_rates$SE_Call = call_rates$SE_Call * 100
call_rates$dummy_order = as.numeric(call_rates$Type)
call_rates = call_rates[order(call_rates$Type, decreasing=T), ]
#call_rates$Breed <- transform(call_rates$Breed, variable=reorder(call_rates$Breed, call_rates$Mean_Call) ) 

ggplot(call_rates, aes(x=reorder(Breed, -dummy_order), y=Mean_Call, fill=Type)) +
  geom_bar(colour="black", stat="identity") +
  geom_errorbar(width=.25, aes(ymin=Mean_Call-SE_Call, ymax=Mean_Call+SE_Call)) +
  ylim(0,100) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  theme_bw(17) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Call Rate Distribution by Breed") + ylab("Mean Call Rate") + xlab("Breed") +
  theme(panel.border = element_rect(colour = "black", size=3), axis.ticks = element_line(size=1.25))
  

# Quality score plots
quality_scores = read.csv("/media/3C52A21052A1CF48/data_backup/Dropbox/Desktop/snpQC/genoQual.csv", stringsAsFactors=F)
bad_markers = read.table("/home/andrew/ultimateExclusionList.txt", stringsAsFactors=F)


p1 <- ggplot(quality_scores, aes(x=AvgGenoQual)) + 
        geom_histogram(aes(y=..density..), binwidth=.05,colour="black", fill="white") +
        geom_density(alpha=.2, fill="#FF6666") +
        theme_bw(18) + ggtitle("Distribution of Marker Quality\nScores for Canine HD Chip") +
        ylab("") + xlab("Average Quality Scores from 1699 Dogs") +
        xlim(0,1) + opts(panel.border = element_rect(colour = "black", size=2))

# Remove problematic plots
quality_scores = quality_scores[-which(quality_scores[,2] %in% bad_markers[,1]), ]
quality_scores = quality_scores[quality_scores$AvgGenoQual > 0.60, ]


p2 <- ggplot(quality_scores, aes(x=AvgGenoQual)) + 
        geom_histogram(aes(y=..density..), binwidth=.05,colour="black", fill="white") +
        geom_density(alpha=.2, fill="#FF6666") +
        theme_bw(18) + ggtitle("Distribution of Marker Quality Scores for Canine HD Chip\nProblematic Markers Removed") +
        ylab("") + xlab("Average Quality Scores from 1699 Dogs") +
        xlim(0,1) + opts(panel.border = element_rect(colour = "black", size=2))

multiplot(p1,p2)

# Informative SNPs across some breeds 
dat = read.csv('breed_allele_freqs.csv',stringsAsFactors=F)
dat[,2] = dat[,2] - dat[,3]
dat$tot = dat[,2] + dat[,3]
dat = transform(dat, Breed=reorder(Breed, tot))

dat_stacked = data.frame(matrix(NA,ncol=3,nrow=(2*length(dat[,1]))))
dat_stacked[,1] = rep(dat[,1],2)
dat_stacked[1:length(dat[,2]),2] = dat[,2]
dat_stacked[(length(dat[,2])+1):length(dat_stacked[,2]),2] = dat[,3]
dat_stacked[,3] = c(rep("MAF > 5%",25),rep("MAF > 20%",25))
names(dat_stacked) = c("Breed","Pct_Of_SNPs","Cutoff")

ggplot(dat_stacked, aes(x = Breed, y = Pct_Of_SNPs)) +
    geom_bar(aes(fill = Cutoff)) +  scale_fill_manual(values=c("#000000", "#FFEF00")) +
    theme_bw(16) +
    ggtitle('Informative SNPs Across Breeds') + xlab('Breed') + ylab('Percent Of SNPs') +
    opts(axis.text.x=theme_text(angle=90, hjust=1)) + 
    opts(panel.border = element_rect(colour = "black", size=2)) 
    






multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

ggplot(df, aes(x=rating)) + geom_histogram(binwidth=.5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(rating, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)