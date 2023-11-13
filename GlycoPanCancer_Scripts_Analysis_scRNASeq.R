######### The transcriptional landscape of glycosylation related genes ########
## Scripts used for the analysis of scRNA-Seq Data
## Ernesto Rodriguez

##### Integration ####

#Datasets used in this manuscript are depicted in the "Data availability" section of Materials and Methods and
#can be downloaded from the sources described.
#Before integration, each dataset was downsampled to 10000 cells.
#The list of cells used per dataset can be found in https://github.com/MolecularCellBiologyImmunology/Glyco_PanCancer
#To reproduce the results depicted in the manuscript, we recommend to generate individual expression matrices per dataset
#using the correspondent cells (e.g., saved as *.cvs)

#Loading data and generating Seurat Object
library(Seurat)
dir <- #Directory where the .cvs files with the scRNA-Seq data are

setwd(dir)
files <- list.files(pattern = ".cvs$")

seurat_list <- list()
for (n in files) {
  d <- read.csv(n)
  d <- CreateSeuratObject(counts = d, min.cells = 3, min.features = 200)
  nam <- gsub(".cvs", "", n)
  print(nam)
  seurat_list[[nam]] <- d
  rm("d")
}

#Normalization and selecting of Variable Features per dataset
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

#Selecting genes to be used during integration and generate PCA for each dataset
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Identify Integration Anchors using the reciprocal PCA (RPCA) method
#As a reference, we used the Breast Cancer dataset.
seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                         anchor.features = features,
                                         reduction = "rpca", 
                                         k.anchor = 20,
                                         reference = 2)

#Select genes to be integrated (all genes across the datasets)
all.features <- Reduce(f=intersect, x = lapply(seurat_list, rownames))

#Integration of the data
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, features.to.integrate = all.features)

#Adding Original Dataset in metadata
a <- NULL
for(n in names(seurat_list)){
  a <- c(a, rep(n, length(colnames(seurat_list[[n]]))))
}
a <- as.factor(a)
names(a) <- colnames(seurat_integrated)
seurat_integrated$Dataset <- a

#Dimensional Reduction and Clustering of Integrated Object
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:50)
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:50)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 2)
DimPlot(seurat_integrated, reduction = "umap", label = T)

#Assign Cluster to Populations
Populations <- read.table("GlycoPanCancer_Cluster_to_Celltype.txt", header = T)
seurat_integrated$Clusters <- Populations$Cell_type[match(tumor$integrated_snn_res.2, Populations$Clusters)]


#### Clean up "Tumor" Cluster
#Selecting "Tumor Cluster" and discard CD45+ datsets
Tumor <- subset(seurat_integrated, idents = "Tumor")
Idents(Tumor) <- Tumor$Class
Tumor <- subset(Tumor, idents = as.character(unique(Idents(Tumor)))[c(-6, -12)])

#Renormalize and reclustering individual datasets
Tumor_list <- SplitObject(Tumor, split.by = "ident")

Tumor_list <- lapply(X = Tumor_list, FUN = function(x) {
  DefaultAssay(x) <- 'RNA'
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, verbose = FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:50)
  x <- FindNeighbors(x, dims = 1:50)
  x <- FindClusters(x, resolution = 2)
})



#Discard clusters that correspond to tissue resident cell types, wrongly clustered cells.
Tumor_list$Astro <- subset(Tumor_list$Astro, idents = c(12, 18), invert = TRUE)
Tumor_list$Breast <- subset(Tumor_list$Breast, idents = c(22, 30, 7, 27, 20), invert = TRUE)
Tumor_list$CRC <- subset(Tumor_list$CRC, idents = c(12), invert = TRUE)
Tumor_list$cSCC <- subset(Tumor_list$cSCC, idents = c(8, 23, 24, 27), invert = TRUE)
Tumor_list$ESCC_CD45neg <- subset(Tumor_list$ESCC_CD45neg, idents = c(32,12,4), invert = TRUE)
Tumor_list$GBM <- subset(Tumor_list$GBM, idents = c(33, 34, 27, 23, 8, 28), invert = TRUE)
Tumor_list$HCC <- subset(Tumor_list$HCC, idents = c(12,22,0,9), invert = TRUE)
Tumor_list$HNSCC <- subset(Tumor_list$HNSCC, idents = c(22, 14, 15, 12), invert = TRUE)
Tumor_list$Lung <- subset(Tumor_list$Lung, idents = c(5,12,15,16,20,22), invert = TRUE)
Tumor_list$Melanoma <- subset(Tumor_list$Melanoma, idents = c(2,3,4,15,19,20,23), invert = TRUE)
Tumor_list$OV <- subset(Tumor_list$OV, idents = c(0:2, 4,6,7,10,12:14,17), invert = TRUE)
Tumor_list$PDAC <- subset(Tumor_list$PDAC, idents = c(21, 11, 27), invert = TRUE)
Tumor_list$UVM <- subset(Tumor_list$UVM, idents = c(16), invert = TRUE)

#Selecting cells in the original dataset
A <- NULL
for (n in names(Tumor_list)){
  A <- c(A, colnames(Tumor_list[[n]]))
}

Tumor_Final <- subset(Tumor, cells = A)

#Separate Breast cancer samples between Triple Negative Breast cancer and Others
Tumor_Final$Clusters <- ifelse(Tumor_Final$Class == "Breast" & Tumor_Final$Condition == "TNBC", "TNBC",
                               ifelse(Tumor_Final$Class == "Breast" & Tumor_Final$Condition != "TNBC",
                                      "Breast", as.character(Tumor_Final$Class)))

Idents(Tumor_Final) <- factor(Tumor_Final$Clusters,
                              levels = rev(c( "cSCC", "ESCC_CD45neg", "HNSCC", "CRC", "PDAC",
                                              "Lung", "OV", "Melanoma", 'UVM', "TNBC", "Breast", "HCC",
                                              "Astro", "Oligo", "GBM")))