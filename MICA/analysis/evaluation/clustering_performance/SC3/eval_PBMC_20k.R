#!/usr/bin/env Rscript
library(SingleCellExperiment)
library(SC3)
library(scater)
library(dplyr)
library(Matrix)
library(data.table)
# library(pdfCluster)


# log1p based?
mat <- read.table(file=paste0("/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/datasets/SilverStd/PBMC_20k/PBMC_20k_MICA_input.txt"),
                  sep="\t", header=TRUE, row.names=1)
mat_t <- transpose(mat)
mat_pow2 <- exp(1)^mat_t - 1
ann <- read.table(file=paste0("/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/datasets/SilverStd/PBMC_20k/PBMC_20k_true_label.txt"),
                  sep="\t", header=TRUE, row.names=1)


# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = mat_pow2,
    logcounts = as.matrix(mat_t)
  ), 
  colData = ann
)


# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
sce <- runPCA(sce)


# plotPCA(sce, colour_by = "label")
sce <- sc3(sce, ks = 10, biology = TRUE, n_cores = 10)

col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
print(adj.rand.index(col_data$label, col_data$sc3_10_clusters))