#!/usr/bin/env Rscript
# Andrew Borgman
# Phylogenetic analysis of SNP data
library(adegenet)
library(multicore)
library(ape)
library(ctc)
library(gplots)

# Read in data
dat = read.PLINK("/home/andrew/projects/whats_the_point/data/all_breeds.raw", chunkSize=1000, multicore=TRUE, n.cores=6)
dat = read.PLINK("/home/andrew/projects/whats_the_point/data/chr_10_breeds.raw", chunkSize=1000, multicore=TRUE, n.cores=6)
dat = read.PLINK("/home/andrew/Dropbox/pointing/rip/chr_24_breeds_top.raw", chunkSize=1000, multicore=TRUE, n.cores=6)
#dat = read.PLINK("/home/andrew/projects/whats_the_point/data/chr_10_breeds_top.raw", chunkSize=1000, multicore=TRUE, n.cores=6)

# Getting phylogenetic trees
head(glNA(dat))
D = dist(as.matrix(dat))
fit = hclust(D)

phy_tree = as.phylo(fit)
plot(phy_tree,type="c",font=1,no.margin=T)

# Writing out tree in Newick topology format
write.table(hc2Newick(fit), file="/home/andrew/Dropbox/pointing/rip/newick_tree_24.txt",row.names=FALSE,col.names=FALSE,quote=F)

datmat = as.matrix(dat)
heatmap.2(datmat, Colv=F,Rowv=F, trace='none', col=c('red','orange','yellow'))