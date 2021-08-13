#/usr/bin/env Rscript

library(SingleCellExperiment)
library(SC3)
library(scater)
library(dplyr)
library(Matrix)
library(data.table)
library(pdfCluster)


mat <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/pollen/Pollen_MICA_input.txt"),
                  sep="\t", header=TRUE, row.names=1)
mat_t <- transpose(mat)
mat_pow2 <- exp(1)^mat_t - 1


# mat <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/pollen/pollen_log2.txt"),
#                   sep="\t", header=TRUE, row.names=1)
# mat_t <- transpose(mat)
# mat_pow2 <- 2^mat_t - 1
ann <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/GoldernStd/pollen/Pollen_true_label.txt"),
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


plotPCA(sce, colour_by = "label")
sce <- sc3(sce, ks = 11, biology = TRUE)

col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
adj.rand.index(col_data$label, col_data$sc3_11_clusters)











pollen = readRDS(file = '/Users/lding/Documents/MICA/Datasets/HPC/GoldenStd/Pollen/pollen.rds')
mat <- pollen@assays[["data"]]@listData[["logcounts"]]
mat_pow <- 2^mat - 1

obs <- read.table(file = '/Users/lding/Documents/MICA/Datasets/HPC/GoldenStd/Pollen/Pollen_obs.txt', sep = '\t', 
                  header = TRUE, row.names = 1)
mat <- mat[,rownames(obs)]
mat_pow <- mat_pow[,rownames(obs)]
cell_types <- pollen@colData[rownames(obs),]


# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = mat_pow,
    logcounts = as.matrix(mat)
  ), 
  colData = cell_types@listData[["cell_type1"]]
)
# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)

sce <- runPCA(sce)

plotPCA(sce, colour_by = "X")
sce <- sc3(sce, ks = 11, biology = TRUE)

col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
adj.rand.index(col_data$X, col_data$sc3_11_clusters)


cells <- pollen@colData@rownames
labels <- pollen@colData@listData[["cell_type1"]]
true_label <- data.frame(cells, labels)
write.table(true_label, file='/Users/lding/Documents/MICA/Datasets/HPC/GoldenStd/Pollen/Pollen_true_label.txt', sep="\t")
