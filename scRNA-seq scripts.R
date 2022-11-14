######### The transcriptional landscape of glycosylation related genes ########
## Scripts used for scRNA-Seq Data
## Ernesto Rodriguez

##### Integration ####

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

seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})


features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                         anchor.features = features,
                                         reduction = "rpca", 
                                         k.anchor = 20,
                                         reference = 2)

a <- NULL
for(n in names(seurat_list)){
  a <- c(a, rep(n, length(colnames(seurat_list[[n]]))))
}
a <- as.factor(a)

all.features <- Reduce(f=intersect, x = lapply(seurat_list, rownames))
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, features.to.integrate = all.features)

names(a) <- colnames(seurat_integrated)
seurat_integrated$Class <- a

tumor <- seurat_integrated
tumor <- ScaleData(tumor, verbose = FALSE)
tumor <- RunPCA(tumor, verbose = FALSE)
tumor <- RunUMAP(tumor, reduction = "pca", dims = 1:30)
tumor <- FindNeighbors(tumor, dims = 1:30)
tumor <- FindClusters(tumor, resolution = 1)
DimPlot(tumor, reduction = "umap", label = T)
