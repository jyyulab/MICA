#/usr/bin/env Rscript

library(SingleCellExperiment)
library(SC3)
library(scater)
library(dplyr)
library(Matrix)
library(data.table)
library(pdfCluster)


author <- 'PBMC_20k'
level <- 'SilverStd'
mat <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_MICA_input.txt"),
                  sep="\t", header=TRUE, row.names=1)
mat_t <- transpose(mat)

s_obj <- CreateSeuratObject(counts = mat_t, project = "Yan", min.cells = 3, min.features = 200)

# Create Seurat object
s_obj <- FindVariableFeatures(s_obj, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(s_obj), 10)
# plot variable features with and without labels
# plot1 <- VariableFeaturePlot(s_obj)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2


# PCA
all.genes <- rownames(s_obj)
s_obj <- ScaleData(s_obj, features = all.genes)
s_obj <- RunPCA(s_obj, features = VariableFeatures(object = s_obj))
# VizDimLoadings(s_obj, dims = 1:2, reduction = "pca")

# Clustering
s_obj <- FindNeighbors(s_obj, dims = 1:10)
s_obj <- FindClusters(s_obj, resolution = 0.4)
s_obj <- RunUMAP(s_obj, dims = 1:10)
DimPlot(s_obj, reduction = "umap")

write.table(s_obj@reductions$umap@cell.embeddings, file='/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/PBMC_20k/PBMC_20k_Seurat_UMAP.txt', sep='\t')
write.table(s_obj@reductions$pca@cell.embeddings, file='/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/PBMC_20k/PBMC_20k_Seurat_PCA.txt', sep='\t')


saveRDS(s_obj, file = paste0("/Users/lding/Documents/MICA/Datasets/HPC/", level ,"/", author,"/", author,"_seurat.rds"))


# Calculate ARI
true_label_file <- paste0('/Users/lding/Documents/MICA/Datasets/HPC/', level, '/', author, '/', author, '_true_label.txt')
true_labels <- read.table(file=true_label_file, sep="\t", header=TRUE, row.names=1)
adj.rand.index(true_labels$label, as.numeric(s_obj$seurat_clusters))


library(aricode)
AMI(true_labels$label, as.numeric(s_obj$seurat_clusters))



# Calculate silhouette
library(scclusteval)
silhouette <- CalculateSilhouette(s_obj, dims=1:50)
mean(silhouette$width)
# 0.2059292
