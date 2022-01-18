#/usr/bin/env Rscript

library(SingleCellExperiment)
library(SC3)
library(scater)
library(dplyr)
library(Matrix)
library(data.table)
library(pdfCluster)


# log1p based
mat <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/Chung_MICA_input.txt"),
                  sep="\t", header=TRUE, row.names=1)
mat_t <- transpose(mat)
<<<<<<< HEAD
colnames(mat_t) <- rownames(mat)
rownames(mat_t) <- colnames(mat)
# ann <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/Chung_true_label.txt"),
#                   sep="\t", header=TRUE, row.names=1)

scgnn_true_label_file <- '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/scGNN/Chung_cell_label.csv'
scgnn_true_label <- read.table(file=scgnn_true_label_file, sep=",", header=TRUE, row.names=1)

mat_t <- mat_t[,rownames(scgnn_true_label)]
mat_pow2 <- exp(1)^mat_t - 1
=======
mat_pow2 <- exp(1)^mat_t - 1
ann <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/Chung_true_label.txt"),
                  sep="\t", header=TRUE, row.names=1)
>>>>>>> b77f5c014b1a3e01faa6131792bee9d39e39b760


# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = mat_pow2,
    logcounts = as.matrix(mat_t)
  ), 
<<<<<<< HEAD
  colData = scgnn_true_label
=======
  colData = ann
>>>>>>> b77f5c014b1a3e01faa6131792bee9d39e39b760
)


# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
sce <- runPCA(sce)


<<<<<<< HEAD
plotPCA(sce, colour_by = "cell_type")
sce <- sc3(sce, ks = 4, biology = TRUE)


col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
adj.rand.index(col_data$cell_type, col_data$sc3_4_clusters)


library(umap)
euclidean_laplacian_umap <- umap(sce@metadata$sc3$transformations$euclidean_laplacian)
write.table(euclidean_laplacian_umap$layout, file='/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Chung/SC3/Chung_SC3_euclidean_laplacian_UMAP.txt', sep='\t')

=======
plotPCA(sce, colour_by = "label")
sce <- sc3(sce, ks = 5, biology = TRUE)

col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
adj.rand.index(col_data$label, col_data$sc3_5_clusters)
>>>>>>> b77f5c014b1a3e01faa6131792bee9d39e39b760
