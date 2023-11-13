######### The transcriptional landscape of glycosylation related genes ########
## Scripts used for the Analysis of TCGA dataset
## Ernesto Rodriguez

##Pre-processing
#EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv can be downloaded from https://gdc.cancer.gov/about-data/publications/pancanatlas
data <- read.delim("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", header = TRUE, row.names = 1, sep = "\t", check.names=FALSE)
data <- as.matrix(data)
rownames(data) <- gsub("\\|.*","", rownames(data))
data <- data[!duplicated(rownames(data)),]
data <- data[-1,]

#Selecting only Tumor samples
data <- data[,as.numeric(substr(colnames(data), 14, 15)) < 10]

#Selecting Glycosylation-related genes
Glycosylation_genes <- read.delim("GlycoPanCancer_Glycogenes.txt", header = 0)
a <- intersect(rownames(data), Glycosylation_genes$V1)
data_glyco_genes <- data[a,]

### Gene set enrichment analysis
#Generate a score for each Glycosylation-related pathway
library(GSVA)
library(qusage)

genesets <- read.gmt("GlycoPanCancer_GlycosPathways.gmt")
ssgsea_res <- gsva(data, genesets, verbose=TRUE, method="ssgsea")
data_glyco_pathways <- as.matrix(ssgsea_res)

### Clustering
library(cluster)
library(ConsensusClusterPlus)

#Selecting 100 most variable genes using Inter Quantile Range (IQR)
mads <- apply(data_glyco_genes,1,IQR)
data_glyco_clust <- data_glyco_genes[rev(order(mads))[1:100],]
data_glyco_clust <- sweep(data_glyco_clust,1, apply(x,1,median,na.rm=T))
data_glyco_clust <- as.matrix(data_glyco_clust)

results_ccp <-
  ConsensusClusterPlus(data_glyco_clust, maxK=25, reps=1000, pItem=0.9,
                       pFeature=1, clusterAlg="hc", title = "Consensus Clustering - TCGA",
                       verbose = TRUE, distance="pearson", seed = 123,
                       innerLinkage="ward.D2", finalLinkage="ward.D2", writeTable = FALSE, plot="pdf")

#Evaluate results of Clustering by calculating silhouette width
#Calculate distance as 1 - pearson correlation
data_dist <- as.dist(1-cor(data_glyco_clust, method = "pearson"))
n <- #Number of Clusters evaluating according to output of ConsensusClusterPlus
sil.res <- silhouette(results_ccp[[n]]$consensusClass, data_dist)
plot(sil.res)


## Run tSNE
library(Rtsne)
tSNE_data <- Rtsne(t(data_glyco_clust), perplexity = 100, max_iter = 2000, verbose = T)
tSNE_data <- as.data.frame(tSNE_data$Y, row.names = colnames(data_glyco_clust))

## Survival Analysis
#The data variable showld be replaced by:
              # data_glyco for the Survival analysis based on glycosylation-related genes
              # data_glyco_pathways for the analysis of glycosylation-related pathways

#Results Clustering 16 Clusters
results_ccp <- as.data.frame(results_ccp[[16]]$consensusClass)
colnames(results_ccp) <- "Cluster"

#Selecting one sample per patient
data <- data[,order(colnames(data))]
data <- data[,!duplicated(substr(colnames(data), 1, 12))]
results_ccp <- results_ccp[colnames(data),]
colnames(data) <- substr(colnames(data), 1, 12) #Patient identifier
rownames(results_ccp) <- substr(rownames(results_ccp), 1, 12)


#TCGA_Clinical variable is the TCGA Clinical data obtained from Liu et al. (https://doi.org/10.1016/j.cell.2018.02.052, Supplementary Table 1)
#Discard samples from TCGA Projects DLBC, PCPG, TGCT and THYM as recommended by Liu et al. (Table 3)
TCGA_Clinical <- TCGA_Clinical [!clin$type %in% c("DLBC", "PCPG", "TGCT", "THYM"),]

#Selecting samples with clinical data and cluster information
samp <- intersect(rownames(TCGA_Clinical), colnames(data))

TCGA_Clinical <- TCGA_Clinical[samp,]
results_ccp <- results_ccp[samp,]
data <- data[,samp]

#Add a 0 before 1 digit numbers
TCGA_Clinical$GlycoCluster <-
  ifelse(substr(results_ccp$Cluster, 2, 2) == "",
         paste(0, results_ccp$Cluster, sep = ""), results_ccp$Cluster)

#Combine Clustering and TCGA Project in a single Variable
TCGA_Clinical$GlycoCluster_Project <-
  paste(TCGA_Clinical$GlycoCluster, clin2$type, sep = "_")
TCGA_Clinical <- subset(TCGA_Clinical, TCGA_Clinical$OS.time != "#N/A")

TCGA_Clinical$OS.time <- as.numeric(TCGA_Clinical$OS.time)/365 #Change from days to years
TCGA_Clinical$OS <- as.numeric(TCGA_Clinical$OS) + 1

#Select condition with 60 or more samples
a <- table(TCGA_Clinical$GlycoCluster_Project)
a <- subset(a, a >= 60)

library(survival)
library(survminer)
library(plyr)

#Loop for Survival Analysis

options(warn=1)
for (n in names(a)) {

#Select data for a given combination of GlycoCluster and TCGA Project
    clin_loop <- subset(TCGA_Clinical, TCGA_Clinical$GlycoCluster_Project == n)
    data_loop <- data[,rownames(clin_loop)]
    Survival.data <- data.frame(OS = TCGA_Clinical$OS, TimeOS =  TCGA_Clinical$OS.time, t(data_loop),
                                row.names = rownames(TCGA_Clinical))
   
#For each gene/score, classify samples between high and low expression
    for (n2 in colnames(Survival.data)[-1:-2]) {
      p2 <- quantile(Survival.data[,n2], probs = seq(0, 1, 1/3), na.rm = T)
      Survival.data[,n2] <- ifelse(Survival.data[,n2] <= p2[2], 1,
                                   ifelse(Survival.data[,n2] > p2[length(names(p2))-1], 2, NA))
      
    }
    
    Survival.data <- Survival.data[,colSums(Survival.data, na.rm = T) > 0]
    
    covariates <- colnames(Survival.data)[-1:-2]
    univ_formulas <- sapply(covariates,
                            function(x) as.formula(paste('Surv(TimeOS, OS)~', x)))
    
    univ_models <- lapply(univ_formulas, function(x){coxph(x, data = Survival.data)})
    
#Extract Beta Coefficient and p value 
    univ_results <- lapply(univ_models,
                           function(x){ 
                             x <- summary(x)
                             p.value<-signif(x$wald["pvalue"], digits=2)
                             wald.test<-signif(x$wald["test"], digits=2)
                             beta<-signif(x$coef[1], digits=3)
                             res<-c(beta, wald.test, p.value)
                             names(res)<-c("beta", "wald.test", "p.value")
                             return(res)
                           })
    

#Test for proportional hazards assumption of a Cox regression model
    univ_zph <- lapply(univ_models, function(x){cox.zph(x)})
    
    univ_results_zph <- lapply(univ_zph,
                               function(x){
                                 p.value<-signif(x$table[1,3], digits=3)
                                 return(p.value)
                               })
    
}



