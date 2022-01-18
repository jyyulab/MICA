#/usr/bin/env Rscript

library(SingleCellExperiment)
library(SC3)
library(scater)
library(dplyr)
library(Matrix)
library(data.table)
library(pdfCluster)


zeisel <- readRDS('/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Zeisel/Zeisel_seurat.rds')

# log1p based?
mat <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Zeisel/Zeisel_MICA_input.txt"),
                  sep="\t", header=TRUE, row.names=1)
mat_t <- transpose(mat)
mat_pow2 <- exp(1)^mat_t - 1
ann <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Zeisel/Zeisel_true_label.txt"),
                  sep="\t", header=TRUE, row.names=1)


# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = mat_pow2,
    logcounts =  as.matrix(mat_t)
  ), 
  colData = ann
)


# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
sce <- runPCA(sce)


plotPCA(sce, colour_by = "label")
sce <- sc3(sce, ks = 7, biology = TRUE)


sc3_plot_silhouette(sce, k = 7)
col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
adj.rand.index(col_data$label, col_data$sc3_7_clusters)



library(umap)
euclidean_laplacian_umap <- umap(sce@metadata$sc3$transformations$euclidean_laplacian)
write.table(euclidean_laplacian_umap$layout, file='/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Zeisel/SC3/Zeisel_SC3_euclidean_laplacian_UMAP.txt', sep='\t')
