#/usr/bin/env Rscript

library(SingleCellExperiment)
library(SC3)
library(scater)
library(dplyr)
library(Matrix)
library(data.table)
library(pdfCluster)


klein <- readRDS('/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Klein/klein.rds')

# log1p based?
# mat <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Klein/Klein_MICA_input.txt"),
#                   sep="\t", header=TRUE, row.names=1)
# mat_t <- transpose(mat)
# mat_pow2 <- exp(1)^mat_t - 1
ann <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Klein/Klein_true_label.txt"),
                  sep="\t", header=TRUE, row.names=1)


# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = klein@assays[["data"]]@listData[["counts"]],
    logcounts = klein@assays[["data"]]@listData[["logcounts"]]
  ), 
  colData = ann
)


# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
sce <- runPCA(sce)


plotPCA(sce, colour_by = "label")
sce <- sc3(sce, ks = 4, biology = TRUE)

col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
adj.rand.index(col_data$label, col_data$sc3_4_clusters)

library(aricode)
AMI(col_data$label, col_data$sc3_4_clusters)



library(umap)
euclidean_laplacian_umap <- umap(sce@metadata$sc3$transformations$euclidean_laplacian)
write.table(euclidean_laplacian_umap$layout, file='/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/SC3/Klein_SC3_euclidean_laplacian_UMAP.txt', sep='\t')

