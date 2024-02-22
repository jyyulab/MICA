# LCA vs MICA node2vec on PBMC14k

library(scLCA)
library(pdfCluster)

data(myscExampleData)
names(myscExampleData)
dim(myscExampleData$datamatrix)
table(myscExampleData$truelabel)


PBMC14k <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_MICA_input_filter_14k.txt', sep = '\t', header = TRUE, row.names = 1)
PBMC14k_t <- t(PBMC14k)
myclust.res <- myscLCA(PBMC14k_t, clust.max=7)


true.label.20k <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt', 
                             sep = '\t', header = T)
true.label.14k <- true.label.20k[true.label.20k$cell %in% names(myclust.res[[1]]),]


# Calculate ARI
adj.rand.index(true.label.14k$label, myclust.res[[1]])
