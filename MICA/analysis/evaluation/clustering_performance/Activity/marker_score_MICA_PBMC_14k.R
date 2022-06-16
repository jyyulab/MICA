# MICA activity-based clustering and marker score

library(scMINER)
library(NetBID2)


mat.20k <- read.table(file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_MICA_input.txt', 
                      sep = '\t', header = TRUE, row.names = 1)
mat.20k <- as.matrix(mat.20k)

filtered_cells <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/filtered_cells_14k.txt', 
                             sep='\t', header = TRUE, row.names = 1)
mat.14k <- mat.20k[filtered_cells$cell,]


eset.14k <- generate.eset(exp_mat=t(mat.14k),
                          phenotype_info=rownames(mat.14k),
                          feature_info=colnames(mat.14k),
                          annotation_info='exp.log2')


cluster.res <- read.table(file = '/Users/lding/Documents/scMINER/PBMC14k_input/activity_6863_drivers/clustering_UMAP_euclidean_24_2.01375.txt', 
                          sep = '\t', header = TRUE)
pData(eset.14k)$cellType <- as.factor(cluster.res$label)
pData(eset.14k)$X <- cluster.res$X
pData(eset.14k)$Y <- cluster.res$Y
exp.log2 <- t(mat.14k)





# clustering_UMAP_euclidean_24_1.49
colors <- c('#00BFC4', '#9590FF', '#00B0F6', '#CD9600', '#F8766D', '#C77CFF', '#FF62BC', '#00BF7D')
MICAplot(input_eset = eset.14k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "cellType", colors = colors, pct = 0.5)


genes_of_interest <-c("FOXP3", "IL2RA", "TIGIT")
genes_of_interest <-c("CD3D", "CD27", "IL7R", "SELL", "CCR7", "IL32", "GZMA", "GZMK", "DUSP2", 
                      "CD8A", "GZMH", "GZMB", "CD79A", "CD79B", "CD86", "CD14")
feature_highlighting(input_eset = eset.14k, target = genes_of_interest, 
                     feature = "gene", ylabel = "log2Exp", x = "X", y = "Y", pct.size = 0.5)




true.label.20k <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt', 
                             sep = '\t', header = T)
true.label.14k <- true.label.20k[true.label.20k$cell %in% cluster.res$ID,]
pData(eset.14k)$trueLabel <- as.factor(true.label.14k$label)


colors <- c('#C3C3C3', '#C3C3C3', '#CD9600', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3')
MICAplot(input_eset = eset.14k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "trueLabel", colors = colors, pct = 0.5)


colors <- c('#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#F8766D', '#C3C3C3', '#C3C3C3')
MICAplot(input_eset = eset.14k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "trueLabel", colors = colors, pct = 0.5)






exp.count <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/Filtered_DownSampled_SortedPBMC_data.csv', 
                        sep=',', row.names = 1, header = TRUE)
exp.count <- exp.count[filtered_cells$cell,]


cluster.res.act <- read.table(file = '/Users/lding/Documents/scMINER/PBMC14k_input/activity_6863_drivers/clustering_UMAP_euclidean_24_2.01375.txt', 
                              sep = '\t', header = TRUE)
cluster.res.exp <- read.table(file = '/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3_1/PBMC_14k/MICA/filter_CTL/clustering_UMAP_euclidean_24_2.72.txt',
                              sep = '\t', header = TRUE)




# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells.
FOXP3.freq <- as.data.frame( table(exp.count[cluster.res.act[cluster.res.act$label == 3,]$ID, 'FOXP3']) )
colnames(FOXP3.freq) <- c('category', 'count.activity')
FOXP3.freq$count.expression <- c(table(exp.count[cluster.res.exp[cluster.res.exp$label == 4,]$ID, 'FOXP3']))

FOXP3.freq$Freq.activity <- FOXP3.freq$count.activity / sum(FOXP3.freq$count.activity)
FOXP3.freq$Freq.expression <- FOXP3.freq$count.expression / sum(FOXP3.freq$count.expression)

mFOXP3.freq <- melt(FOXP3.freq, id=c('category'))
mFOXP3.freq <- mFOXP3.freq[c(8,9,11,12),]
mFOXP3.freq$label_ypos <- c(0.05, 0.006, 0.009, 0.0003)


# Create the barplot
ggplot(data=mFOXP3.freq, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=value), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_brewer(palette='Set2')+
  theme_minimal()




# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells.
IL2RA.freq <- as.data.frame( table(exp.count[cluster.res.act[cluster.res.act$label == 3,]$ID, 'IL2RA']) )
colnames(IL2RA.freq) <- c('category', 'count.activity')
IL2RA.freq$count.expression <- c(table(exp.count[cluster.res.exp[cluster.res.exp$label == 4,]$ID, 'IL2RA']))

IL2RA.freq$Freq.activity <- IL2RA.freq$count.activity / sum(IL2RA.freq$count.activity)
IL2RA.freq$Freq.expression <- IL2RA.freq$count.expression / sum(IL2RA.freq$count.expression)

mIL2RA.freq <- melt(IL2RA.freq, id=c('category'))
mIL2RA.freq <- mIL2RA.freq[c(12,13,14,15,17,18,19,20),]
mIL2RA.freq$label_ypos <- c(0.05, 0.006, 0.009, 0.0003)


# Create the barplot
ggplot(data=mIL2RA.freq, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=value), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_brewer(palette='Set2')+
  theme_minimal()




# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells.
TIGIT.freq <- as.data.frame( table(exp.count[cluster.res.act[cluster.res.act$label == 3,]$ID, 'TIGIT']) )
colnames(TIGIT.freq) <- c('category', 'count.activity')
TIGIT.freq$category <- as.vector(TIGIT.freq$category)
TIGIT.freq[nrow(TIGIT.freq) + 1, ] <- c("4", 1)
TIGIT.freq <- TIGIT.freq[c(1, 2 , 3, 4, 6, 5),]
TIGIT.freq$count.expression <- c(table(exp.count[cluster.res.exp[cluster.res.exp$label == 4,]$ID, 'TIGIT']))

TIGIT.freq$count.activity <- as.numeric(TIGIT.freq$count.activity)
TIGIT.freq$Freq.activity <- TIGIT.freq$count.activity / sum(TIGIT.freq$count.activity)
TIGIT.freq$Freq.expression <- TIGIT.freq$count.expression / sum(TIGIT.freq$count.expression)

mTIGIT.freq <- melt(TIGIT.freq, id=c('category'))
mTIGIT.freq <- mTIGIT.freq[c(14,15,16,17,18,20,21,22,23,24),]
mTIGIT.freq$label_ypos <- c(0.06, 0.05, 0.006, 0.009, 0.0003)
mTIGIT.freq$category <- as.factor(mTIGIT.freq$category)
mTIGIT.freq$value <- as.numeric(mTIGIT.freq$value)


# Create the barplot
ggplot(data=mTIGIT.freq, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=value), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_brewer(palette='Set2')+
  theme_minimal()








set.seed(20190708)
genes <- paste("gene",1:1000,sep="")
x <- list(
  A = sample(genes,300), 
  B = sample(genes,525), 
  C = sample(genes,440),
  D = sample(genes,350)
)


if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")


if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library("ggvenn")
names(x) <- c("Stage 1","Stage 2","Stage 3", "Stage 4")
ggvenn(
  x, columns = c("Stage 1", "Stage 2", "Stage 3"),
  stroke_size = 0.5
)



cluster.res.act3 <- cluster.res.act[cluster.res.act$label == 3,]
cluster.res.exp4 <- cluster.res.exp[cluster.res.exp$label == 4,]
true.label.14k.Treg <- true.label.14k[true.label.14k$label == 'CD4+/CD25 T Reg',]

x <- list(
  activity = cluster.res.act3$ID,
  expression = cluster.res.exp4$ID,
  true.label = true.label.14k.Treg$cell
)

ggvenn(
  x, columns = c("activity", "expression", "true.label"),
  stroke_size = 0.5
)




y <- data.frame(abundance = c('activity', 'expression'), overlap_percentage = c(63, 60))

# Change the colors manually
p <- ggplot(data=y, aes(x=abundance, y=overlap_percentage, fill=abundance)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()
# Use custom colors
p + scale_fill_manual(values=c('#999999','#E69F00'))
# Use brewer color palettes
# p + scale_fill_brewer(palette="Blues")

