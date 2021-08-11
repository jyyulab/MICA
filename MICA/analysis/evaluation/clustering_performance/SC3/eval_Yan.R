#/usr/bin/env Rscript

library(SingleCellExperiment)
library(SC3)
library(scater)
library(dplyr)
library(Matrix)
library(data.table)
library(pdfCluster)


yan <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/GoldernStd/Yan/Yan_MICA_input.txt"),
                  sep="\t", header=TRUE, row.names=1)
yan_t <- transpose(yan)
yan_pow2 <- 2^yan_t - 1
ann <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/GoldernStd/Yan/Yan_true_label.txt"),
                  sep="\t", header=TRUE, row.names=1)


# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = yan_pow2,
    logcounts = as.matrix(yan_t)
  ), 
  colData = ann
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