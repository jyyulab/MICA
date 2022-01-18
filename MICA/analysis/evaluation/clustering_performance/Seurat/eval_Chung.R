#/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(pdfCluster)


# Read raw TPM matrix
raw_mat <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/raw/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix_no_pool_T.txt',
                      sep=",", header=TRUE, row.names=1)
raw_mat_t <- transpose(raw_mat)
rownames(raw_mat_t) <- colnames(raw_mat)
colnames(raw_mat_t) <- rownames(raw_mat)



s_obj <- CreateSeuratObject(counts = raw_mat_t, project = "Chung", min.cells = 3, min.features = 50)

s_obj[["percent.mt"]] <- PercentageFeatureSet(s_obj, pattern = "^MT-")
VlnPlot(s_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(s_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(s_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

s_obj <- subset(s_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
s_obj <- NormalizeData(s_obj, normalization.method = "LogNormalize", scale.factor = 10000)

s_obj <- FindVariableFeatures(s_obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(s_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(s_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(s_obj)
s_obj <- ScaleData(s_obj, features = all.genes)
s_obj <- RunPCA(s_obj, features = VariableFeatures(object = s_obj))

s_obj <- FindNeighbors(s_obj, dims = 1:10)
s_obj <- FindClusters(s_obj, resolution = 0.08)
s_obj <- RunUMAP(s_obj, dims = 1:10)
DimPlot(s_obj, reduction = "umap")

# Calculate ARI
author <- 'Chung'
level <- 'SilverStd'
true_label_file <- paste0('/Users/lding/Documents/MICA/Datasets/HPC/', level, '/', author, '/', author, '_true_label.txt')
true_labels <- read.table(file=true_label_file, sep="\t", header=TRUE, row.names=1)
seurat_clusters <- s_obj$seurat_clusters[rownames(true_labels)]
adj.rand.index(true_labels$label, as.numeric(seurat_clusters))

write.csv(s_obj@reductions$umap@cell.embeddings, file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/Chung_k5_seurat_umap.csv')









# Read input matrix
author <- 'Chung'
level <- 'SilverStd'
pp_mat <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/", level ,"/", author, "/", author, "_MICA_input.txt"),
                     sep="\t", header=TRUE, row.names=1)
pp_mat_t <- transpose(pp_mat)
colnames(pp_mat_t) <- rownames(pp_mat)
rownames(pp_mat_t) <- colnames(pp_mat)

scgnn_true_label_file <- '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/scGNN/Chung_cell_label.csv'
scgnn_true_label <- read.table(file=scgnn_true_label_file, sep=",", header=TRUE, row.names=1)

pp_mat_t <- pp_mat_t[,rownames(scgnn_true_label)]

s_obj <- CreateSeuratObject(counts = pp_mat_t, project = "Chung", min.cells = 3, min.features = 50)

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
s_obj <- FindClusters(s_obj, resolution = 0.15)
s_obj <- RunUMAP(s_obj, dims = 1:10)
DimPlot(s_obj, reduction = "umap")

write.table(s_obj@reductions$umap@cell.embeddings, file='/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Chung/Chung_Seurat_UMAP.txt', sep='\t')
write.table(s_obj@reductions$pca@cell.embeddings, file='/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Chung/Chung_Seurat_PCA.txt', sep='\t')


adj.rand.index(scgnn_true_label$cell_type, as.numeric(s_obj$seurat_clusters))


saveRDS(s_obj, file = paste0("/Users/lding/Documents/MICA/Datasets/HPC/", level ,"/", author,"/", author,"_seurat.rds"))

# Calculate ARI
true_label_file <- paste0('/Users/lding/Documents/MICA/Datasets/HPC/', level, '/', author, '/', author, '_true_label.txt')
true_labels <- read.table(file=true_label_file, sep="\t", header=TRUE, row.names=1)
adj.rand.index(true_labels$label, as.numeric(s_obj$seurat_clusters))


# Calculate silhouette
library(scclusteval)
silhouette <- CalculateSilhouette(s_obj, dims=1:50)
mean(silhouette$width)
# 4.19E-05


