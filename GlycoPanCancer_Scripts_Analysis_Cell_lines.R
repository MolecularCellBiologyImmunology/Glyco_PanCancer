library(plyr)
library(cluster)

#For each cell line dataset, we changed the names of the cell lines according to the Supplementary Table 4.

###CLUSTERING###

d <- #Load dataset matrix
Glyco <- read.delim("Glycosylation_related_genes.txt", header = FALSE, sep = "\t", check.names=FALSE)
genes <- intersect(rownames(d), Glyco[,1])
x <- d[genes,]
mads <- apply(x,1,IQR)
x <- x[rev(order(mads))[1:100],]
x <- sweep(x,1, apply(x,1,median,na.rm=T))
x <- as.matrix(x)

res.klijn <-
  ConsensusClusterPlus(x, maxK=10, reps=1000, pItem=0.85, pFeature=1,
                       clusterAlg="hc", title = "Consensus Clustering - Klijn",
                       verbose = TRUE, distance="pearson", seed = 123,
                       innerLinkage="ward.D2", finalLinkage="ward.D2",
                       writeTable = FALSE, plot="pdf")

#Evaluate Clustering using silhouette width
xx <- as.dist(1-cor(x, method = "pearson"))
sil.res <- silhouette(res[[5]]$consensusClass, xx)
plot(sil.res)

#Finding Consensus Clusters
#Performing a Cluster of Clusters

n <- 5 #Defined 5 clusters per dataset
res.klijn <- as.matrix(res.klijn[[n]]$consensusClass)
res36133 <- as.matrix(res36133[[n]]$consensusClass)
res3610 <- as.matrix(res3610[[n]]$consensusClass)
res783 <- as.matrix(res783[[n]]$consensusClass)

res_combined <- rbind.fill(as.data.frame(t(res.3)), as.data.frame(t(res36133.3)), as.data.frame(t(res783.3)), as.data.frame(t(res3610.3)))

rownames(res_combined) <- c("Klijn", "GSE36133", "MTAB783", "MTAB3610")
res_combined <- t(res_combined)
mod1 <- model.matrix(~0+factor(res_combined[,1]))
rownames(mod1) <- names(na.omit(res_combined[,1]))

mod2 <- model.matrix(~0+factor(res_combined[,2]))
rownames(mod2) <- names(na.omit(res_combined[,2]))

mod3 <- model.matrix(~0+factor(res_combined[,3]))
rownames(mod3) <- names(na.omit(res_combined[,3]))

mod4 <- model.matrix(~0+factor(res_combined[,4]))
rownames(mod4) <- names(na.omit(res_combined[,4]))

res_combined_split <- rbind.fill(as.data.frame(t(mod1)), as.data.frame(t(mod2)), as.data.frame(t(mod3)), as.data.frame(t(mod4)))
rownames(res_combined_split) <- c("Klijn-A", "Klijn-B", "Klijn-C", "Klijn-D", "Klijn-E",
                   "GSE36133-A", "GSE36133-B", "GSE36133-C", "GSE36133-D", "GSE36133-E",
                   "MTAB783-A", "MTAB783-B", "MTAB783-C", "MTAB783-D", "MTAB783-E",
                   "MTAB3610-A", "MTAB3610-B", "MTAB3610-C", "MTAB3610-D", "MTAB3610-E")

res_combined_split[is.na(res_combined_split)] <- 0
res_combined_split <- as.matrix(res_combined_split)

#Performing Clustering
res.res<- ConsensusClusterPlus(res_combined_split, maxK=10, reps=1000, pItem=0.9, pFeature=1, clusterAlg="hc", title = "Consensus Clustering - All", verbose = TRUE,
                               distance="pearson", seed = 123, innerLinkage="ward.D2", finalLinkage="ward.D2",  writeTable = FALSE, plot="pdf")



res_combined <- as.data.frame(res_combined)
res_combined$Cluster <- res.res[[n]]$consensusClass

###Differential Expression
setwd("C:/Users/ernes/Dropbox/Data R/Cell lines/Clustering/Limma")
#Limma
library(edgeR)
library(limma)

##Klijn
res_combined <- subset(res_combined, res_combined$Klijn > 0)
res_combined <- res_combined[order(rownames(res_combined)),]
res_combined$Cluster2 <- chartr("1234567890", "ABCDEFGHI2", res_combined$Cluster)
d.limma <- d[ ,rownames(res_combined)]


#For each cluster, perform DGA against all the other.
design <- model.matrix(~0+factor(ifelse(res_combined$Cluster == 1, 1, 2)))
colnames(design) <- c("GlycoA","GlycoB")
fit <- lmFit(d.limma, design)
contrasts <- makeContrasts(GlycoA-GlycoB, levels=design)
contr.fit <- eBayes(contrasts.fit(fit,contrasts))
Topsgenes.A <- topTable(contr.fit,coef=1, number=60000, adjust="fdr")
