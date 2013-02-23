library(ggplot2)
library(ape)
library(gplots)
library(ctc)

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

p = ggplot(dat_stacked, aes(x = Breed, y = Pct_Of_SNPs)) +
  geom_bar(aes(fill = Cutoff)) +  scale_fill_manual(values=c("#000000", "#FFEF00"))
p = p + opts(title='Informative SNPs in Each Breed') + opts(axis.text.x=theme_text(angle=-90, hjust=0)) + xlab('Breed') + ylab('Percent Of SNPs')

# Phylogenetic tree attempt
#phy = read.table("kinship_matrix",sep="\t",header=T)
phy = read.table("plink.mibs",sep=" ",header=F)
blast = read.csv('/home/andrew/new_master_blaster.csv',stringsAsFactors=F)

# Getting breeds for each of the dogs
nms = data.frame(phy[,1])
mm = merge(nms,blast,by.x=1,by.y=2,sort=F)
mm$row_name = paste(mm[,1],mm[,3],sep="_")

rownames(phy) = mm$row_name  #phy[,1]
phy = phy[,-1]
phy.mat = data.matrix(phy)
phy.mat = -phy.mat

# Heatmap for fun
#heatmap.2(phy.mat)
#test = phy.mat * lower.tri(phy.mat)

# Cluster dendrogram
testes = as.dist(phy.mat)
fit = hclust(testes)
plot(fit)

phy_tree = as.phylo(fit)
plot(phy_tree,type="c",font=1,no.margin=T)


# Writing out tree in Newick topology format
write.table(hc2Newick(fit), file="newick_tree_all_hot_region.txt",row.names=FALSE,col.names=FALSE,quote=F)

# Fixing labels to make nicer looking images
for(z in 1:length(fit$labels)){
  cat(z,"\n")
  lab_str = fit$labels[z]
  new_lab = strsplit(lab_str,"_")[[1]][2]
  fit$labels[z] = new_lab
  
}

#########################################
#         CHR 25 Only                   #
#########################################
system("")