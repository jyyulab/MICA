#!/usr/bin/env Rscript
library(SingleCellExperiment)
library(SC3)
library(scater)
library(dplyr)
library(Matrix)
library(data.table)
library(pdfCluster)


# log1p based?
mat <- read.table(file=paste0("/Users/lding/Documents/scMINER/PBMC14k_input/PBMC_20k_MICA_input_filter_14k.txt"),
                  sep="\t", header=TRUE, row.names=1)
mat_t <- t(as.matrix(mat))
mat_pow2 <- exp(1)^mat_t - 1
ann <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt"),
                  sep="\t", header=TRUE, row.names=1)
ann <- ann[rownames(mat),]


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


plotPCA(sce, colour_by = "X")
sce <- sc3(sce, ks = 10, biology = TRUE)

col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
adj.rand.index(col_data$label, col_data$sc3_5_clusters)


save.image(file='/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/tests/SilverStd/PBMC_20k/SC3/PBMC_20k_SC3.RData')
# plotPCA(sce, colour_by = "label")
sce <- sc3(sce, ks = 10, biology = TRUE, n_cores = 10)

col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
print(adj.rand.index(col_data$label, col_data$sc3_10_clusters))
