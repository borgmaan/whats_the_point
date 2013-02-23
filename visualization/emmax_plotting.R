#Script to generate plots from GWAS study
#################################
#	 Making Manhattan Plot		#
#################################
#png("emmax_manhattan.png",width=1500,height=500)

setwd("/Users/andrewborgman/Dropbox/pointing/rip/")
#Reading in association, map, and bad_snp files
assoc_raw <- read.table("emmax_results.ps",stringsAsFactors=F,header=T)
map <- read.table("canFam3.map",stringsAsFactors=F)
bad <- read.table("ultimateExclusionList.txt",stringsAsFactors=F)
#bad_2 <- read.table("exclude_from_emmax",stringsAsFactors=F)
#bad <- rbind(bad,bad_2)
names(assoc_raw) = c("SNP","Beta","SE(Beta)","P")
names(map) = c("CHR","SNP","Stupid","BP")

# Merging and excluding bad snps
assoc = merge(assoc_raw,map,by.x=1,by.y=2)
assoc = assoc[-which(assoc$SNP %in% bad[,1]),]

#Significance!
assoc$logP <- -log10(assoc$P)
#Writing out significant snps with correct base pair position
sel <- which(assoc$logP > quantile(assoc$logP,.999 ,na.rm=T))
sig <- assoc[sel,]
write.csv(sig,'sig.snps.emmax.csv',row.names=F)
 
#Making chromosomes numeric X->39 Y->40
sel <- which(assoc$CHR == 'X' | assoc$CHR == 'Y')
assoc$CHR[sel] <- 39
assoc$CHR <- as.numeric(assoc$CHR)
#Plotting!
assoc$BP <- assoc$BP / 100
lastbase <- 0
newbase <- 0
ticks <- NULL
for(i in 1:length(unique(assoc$CHR))){
  cat(i)
  if(i == 1){
    newbase <- range(assoc[assoc$CHR == i,]$BP)[2]
    ticks <- c(ticks, (newbase - lastbase) / 2)
    lastbase <- newbase
  }
  else{
    assoc[assoc$CHR == i,]$BP <- assoc[assoc$CHR == i,]$BP + lastbase
    newbase <- range(assoc[assoc$CHR == i,]$BP)[2]
    ticks <- c(ticks, ((newbase - lastbase) / 2)+lastbase)
    lastbase <- newbase
  }
}

#Specifying colors
cols = rep(c('#000000','#585858'),40)

par(mar=c(5.1,6,4.3,2.1))


#Plotting points
tit = "EMMAX GWAS: 132 Pointers vs. 246 Non-Pointers\n6 Randomly Sampled Dogs per Breed"

# Hack Job
for(i in 1){
	if (ceiling((range(assoc$logP,na.rm=T)[2])+2) > 8){
		plot(logP~BP,ylim=c(0,9), xlab='Chromosome',ylab=expression(-log[10](italic(p))),main=tit,axes=F,frame.plot=T, cex.lab = 1.8, cex.main=1.8, col="#FFFFFF", data=assoc)
	}
	else { 
		plot(logP~BP,ylim=c(0,8), xlab='Chromosome',ylab=expression(-log[10](italic(p))),main=tit,axes=F,frame.plot=T, cex.lab = 1.8, cex.main=1.8, data=assoc)
	}
}
axis(1, at=ticks, cex.axis = 1.5, lab=c(seq(1,38),"X"))
yticks <- 0:9
axis(2,at=yticks,cex.axis = 1.5,lab=yticks, lwd=4)
for (i in 1:39) {
	with(assoc[assoc$CHR==i, ],points(BP, logP, col=cols[i],pch=20))  
}
abline(h = -log10(.05 / length(assoc[,1])),col='red',lwd=5)
box(lwd=4)

#Tagging hot loci
hot = assoc[assoc$logP > -log10(.05 / length(assoc[,1])), ]
with(hot,points(BP,logP,col='#006400',pch=20))
dev.off()



