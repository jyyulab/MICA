#!/usr/bin/env Rscript


library(scMINER)
library(NetBID2)


mat.20k <- read.table(file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_MICA_input.txt', 
                      sep = '\t',
                      header = TRUE,
                      row.names = 1)
mat.20k <- as.matrix(mat.20k)

cluster.res <- read.table(file = '/Users/lding/Git/scMINER/docs/docs/images/clustering_UMAP_euclidean_24_1.82212.txt',
                          sep = '\t', header = TRUE, stringsAsFactors = F)

mat.14k <- mat.20k[cluster.res$ID,]

eset.14k <- generate.eset(exp_mat=t(mat.14k),
                          phenotype_info=rownames(mat.14k),
                          feature_info=colnames(mat.14k),
                          annotation_info='exp.log2')


pData(eset.14k)$cellType <- as.factor(cluster.res$label)
pData(eset.14k)$X <- cluster.res$X
pData(eset.14k)$Y <- cluster.res$Y
exp.log2 <- t(mat.14k)


# clustering_UMAP_euclidean_24_3.0
colors <- c('#00BFC4', '#7CAE00',  '#00BF7D', '#00B0F6',  '#CD9600', '#F8766D', '#C77CFF')
MICAplot(input_eset = eset.14k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "cellType", colors = colors, pct = 0.5)
