library(SingleCellExperiment)
library(scater)
library(scran)
library(Seurat)

###########################Male########################################

#cont'd from the object from https://github.com/CB-postech/Chicken_GermCell/
seurat_male <- readRDS("seurat_male_QC_CCcorrected.rds")

#subset E2.5 to E8
dr_male <- seurat_male[, seurat_male$condition %in% c("E2.5_Male", "E6_Male", "E8_Male")]
dr_male <- dr_male[rowSums(dr_male@assays$RNA@counts) != 0, ]

sce_male <- as.SingleCellExperiment(dr_male)

#Normalization
clusters <- quickCluster(sce_male, method="igraph")
sce_male <- computeSumFactors(sce_male, clusters=clusters)
sce_male <- logNormCounts(sce_male)

#HVG selection
dec <- modelGeneVar(sce_male)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)
top.hvgs <- getTopHVGs(dec, n = 750)

dr_male <- as.Seurat(sce_male)
dr_male@assays$RNA@var.features <- top.hvgs

dr_male@reductions$PCA <- NULL
dr_male@reductions$UMAP <- NULL

dr_male <- ScaleData(dr_male, features = top.hvgs)

#dimensionality reduction
dr_male <- RunPCA(dr_male, npcs = 50, weight.by.var = F)
plot(dr_male@reductions$pca@stdev)
PCA = 15

set.seed(10)
dr_male <- FindNeighbors(dr_male, reduction = "pca", dims = 1:PCA) %>%
            FindClusters(resolution = 0.8) %>%
            RunUMAP(reduction = "pca", dims = 1:PCA)



#####################Female###############################

#cont'd from the object from https://github.com/CB-postech/Chicken_GermCell/
seurat_female <- readRDS("seurat_female_QC_CCcorrected.rds")

#subset E2.5 to E8
dr_female <- seurat_female[, seurat_female$condition %in% c("E2.5_Female", "E6_Female", "E8_Female")]
dr_female <- dr_female[rowSums(dr_female@assays$RNA@counts) != 0, ]

sce_female <- as.SingleCellExperiment(dr_female)

#Normalization
clusters <- quickCluster(sce_female, method="igraph")
sce_female <- computeSumFactors(sce_female, clusters=clusters)
sce_female <- logNormCounts(sce_female)

#HVG selection
dec <- modelGeneVar(sce_female)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)
top.hvgs <- getTopHVGs(dec, n = 750)

dr_female <- as.Seurat(sce_female)
dr_female@assays$RNA@var.features <- top.hvgs

dr_female@reductions$PCA <- NULL
dr_female@reductions$UMAP <- NULL

dr_female <- ScaleData(dr_female, features = top.hvgs)

#dimensionality reduction
dr_female <- RunPCA(dr_female, npcs = 50, weight.by.var = F)
plot(dr_female@reductions$pca@stdev)
PCA = 15

set.seed(10)
dr_female <- FindNeighbors(dr_female, reduction = "pca", dims = 1:PCA) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(reduction = "pca", dims = 1:PCA)

saveRDS(dr_male, "Male_2.5-8.rds")
saveRDS(dr_female, "Female_2.5-8.rds")