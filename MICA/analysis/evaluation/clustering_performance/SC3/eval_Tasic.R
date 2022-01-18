#/usr/bin/env Rscript

library(SingleCellExperiment)
library(SC3)
library(scater)
library(dplyr)
library(Matrix)
library(data.table)
library(pdfCluster)

# log1p based?
mat <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Tasic/Tasic_MICA_input.txt"),
                  sep="\t", header=TRUE, row.names=1)
mat_t <- transpose(mat)
mat_pow2 <- exp(1)^mat_t - 1
ann <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Tasic/Tasic_true_label.txt"),
                  sep="\t", header=TRUE, row.names=1)

label <- ann[rownames(ann) %in% rownames(mat),]
rownames2 <- rownames(ann)[rownames(ann) %in% rownames(mat)]
ann2 <- data.frame(label)
rownames(ann2) <- rownames2
ann2 <- ann[rownames(mat),]


# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = mat_pow2,
    logcounts = as.matrix(mat_t)
  ), 
  colData = ann2
)


# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
sce <- runPCA(sce)


plotPCA(sce, colour_by = "label")
sce <- sc3(sce, ks = 8, biology = TRUE)

col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
adj.rand.index(col_data$label, col_data$sc3_8_clusters)
