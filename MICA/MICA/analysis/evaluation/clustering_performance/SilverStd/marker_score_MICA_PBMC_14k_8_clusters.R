# MICA clustering and marker score

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




cluster.res <- read.table(file = '/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3_1/PBMC_14k/MICA/filter_CTL/clustering_UMAP_euclidean_24_2.72.txt', 
                          sep = '\t', header = TRUE)
# cluster.res <- read.table(file = '/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3_1/PBMC_14k/MICA/filter_CTL/clustering_UMAP_euclidean_24_1.49.txt', 
#                           sep = '\t', header = TRUE)
# cluster.res <- read.table(file = '/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3_1/PBMC_14k/MICA/filter_CTL/clustering_UMAP_euclidean_24_1.35.txt', 
#                           sep = '\t', header = TRUE)
# cluster.res <- read.table(file = '/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3_1/PBMC_14k/MICA/filter_CTL/clustering_UMAP_euclidean_24_3.0.txt',
#                           sep = '\t', header = TRUE)
# cluster.res <- read.table(file = '/Users/lding/Git/scMINER/docs/docs/images/clustering_UMAP_euclidean_24_1.82212.txt',
#                           sep = '\t', header = TRUE)
# mat.14k <- mat.14k[cluster.res$ID,]
# eset.14k <- generate.eset(exp_mat=t(mat.14k),
#                           phenotype_info=rownames(mat.14k),
#                           feature_info=colnames(mat.14k),
#                           annotation_info='exp.log2')



pData(eset.14k)$cellType <- as.factor(cluster.res$label)
pData(eset.14k)$X <- cluster.res$X
pData(eset.14k)$Y <- cluster.res$Y
exp.log2 <- t(mat.14k)


# clustering_UMAP_euclidean_24_3.0
colors <- c('#00BFC4', '#7CAE00',  '#00BF7D', '#00B0F6',  '#CD9600', '#F8766D', '#C77CFF')
MICAplot(input_eset = eset.14k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "cellType", colors = colors, pct = 0.5)


# clustering_UMAP_euclidean_24_1.49
colors <- c('#00BFC4', '#9590FF', '#00B0F6', '#CD9600', '#F8766D', '#C77CFF', '#FF62BC', '#00BF7D', '#7CAE00')
MICAplot(input_eset = eset.14k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "cellType", colors = colors, pct = 0.5)


# clustering_UMAP_euclidean_24_1.35
colors <- c('#A3A500', '#9590FF', '#00BFC4', '#00B0F6', '#CD9600', '#F8766D', '#C77CFF', '#FF62BC', '#00BF7D', '#7CAE00')
MICAplot(input_eset = eset.14k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "cellType", colors = colors, pct = 0.5)



true.label.20k <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt', 
                             sep = '\t', header = T)
true.label.14k <- true.label.20k[true.label.20k$cell %in% cluster.res$ID,]
pData(eset.14k)$trueLabel <- as.factor(true.label.14k$label)






colors <- c( '#00BFC4', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3')
MICAplot(input_eset = eset.14k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "trueLabel", colors = colors, pct = 0.5)

colors <- c( '#C3C3C3', '#C3C3C3', '#CD9600', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3')
MICAplot(input_eset = eset.14k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "trueLabel", colors = colors, pct = 0.5)

library(pdfCluster)
adj.rand.index(pData(eset.14k)$trueLabel, pData(eset.14k)$cellType)






# Create data
data <- data.frame(
  name=c("MICA_CP10K", "MICA_CPM"),
  ARI=c(0.8130037, 0.8401542)
)

library(ggplot2)
library(scales)
# 3: Using RColorBrewer
ggplot(data, aes(x=as.factor(name), y=ARI, fill=as.factor(name))) + 
  geom_bar(stat = "identity", width=0.7) + 
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position="none") +
  theme_bw() + scale_y_continuous(limits=c(0.5, 0.9), oob = rescale_none)






# Heatmap
library(openxlsx)
ref<-read.xlsx("/Users/lding/Git/scMINER/docs/tests/Ref/PBMC_markers.xlsx", sheet="Markers_10_cell_types_PBMC_20k")
head(ref)

celltypes <- c("CD4TCM", "CD4TReg", "CD4NaiveT", "CD8NaiveCTL", "Bcell", "NK", "Monocyte")

ac <- matrix(NA, nrow=ncol(exp.log2), ncol=length(celltypes), dimnames = list(colnames(exp.log2), celltypes))

for(i in 1:length(celltypes)){
  cat(i,"\n")
  ref.sel<-dplyr::filter(ref,celltype==celltypes[i])
  n <- length(unique(ref.sel$markers))
  print(celltypes[i])
  print(n)
  
  if(n>1){
    mat<-t(exp.log2[ref.sel$markers,])%*%as.vector(ref.sel$weight)
    ac[,i]<-mat[,1]/n
  }else if (n==1){
    ac[,i]<-exp.log2[ref.sel$markers,]
  }
}

colnames(ac) <- celltypes

ac_norm<-apply(ac,2,scale) #column normalization

n_mtx<-(ac>0.5)
group_name = 'cellType'
df_n<-data.frame(label=pData(eset.14k)[,group_name],n_mtx)
df_n<-aggregate(.~label,df_n,mean)
library(reshape2)
df_n_melt<-melt(df_n,id.vars = "label")





df <- data.frame(label=pData(eset.14k)[,group_name],ac_norm);
df <- df[,colSums(is.na(df))<nrow(df)]; # remove NA columns
df <- aggregate(.~label,df,mean)
input<-t(apply(df[,-1],1,scale))  # row normalization
input<-as.data.frame(cbind(df[,1],input))
rownames(input)<-rownames(df)
colnames(input)<-colnames(df)
df_melt<-melt(input,id.vars = "label")

d <- cbind(df_melt,df_n_melt[,3])
colnames(d)<-c("Cluster", "CellType", "MarkerScore","ExpressionPercentage")
d$Cluster<-as.factor(d$Cluster)





pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_14k/MICA/Heatmap_marker_score_Adam_markers_2.72_raw.pdf", 
    height=10, width=10)

h.mat <- t(input[,-1])
# h.mat <- h.mat[,c(8, 4, 7, 2, 3, 6, 5, 1)]

Heatmap(as.matrix(h.mat), col = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)), name = 'MarkerScore',
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,)
dev.off()









# Draw Sankey plot
# Read true label file
true.label.20k <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt', 
                             sep = '\t', header = T)
true.label.14k <- true.label.20k[true.label.20k$cell %in% cluster.res$ID,]

true.label.14k$MICA.cluster <- as.factor(cluster.res$label)
true.label.14k <- true.label.14k[true.label.14k$cell %in% pbmc@assays$RNA@counts@Dimnames[[2]],]
true.label.14k$Seurat.cluster <- pbmc@meta.data$seurat_clusters

# test.factor.MICA <- factor(true.label.20k$MICA.cluster, levels=c(1,2,3,4,5,6,7,8,9,10), labels=c(9, 2, 4, 8, 6, 3, 7, 1, 5, 10))
# test.factor.Seurat <- factor(true.label.20k$Seurat.cluster, levels=c(0,1,2,3,4,5,6,7,8,9), labels=c(1, 8, 2, 3, 4, 5, 6, 7, 9, 10))

test.factor.MICA <- factor(true.label.14k$MICA.cluster, levels=c(1,2,3,4,5,6,7,8), 
                           labels=c('M.Monocyte', 'M.CD8NaiveCTL', 'M.Bcell1', 'M.CD4TReg', 'M.NK', 'M.Bcell2',
                                    'M.CD4NaiveT', 'M.CD4TCM'))
test.factor.Seurat <- factor(true.label.14k$Seurat.cluster, levels=c(0,1,2,3,4,5,6,7), 
                             labels=c('S.CD4TCM', 'S.CD4NaiveT', 'S.NK', 'S.Bcell1',  'S.Monocyte1',  'S.CD8NaiveCTL',
                                      'S.Bcell2', 'S.Monocyte2'))
test.factor.true.label <- as.factor(true.label.14k$label)
test.factor.true.label <- factor(test.factor.true.label, levels=levels(test.factor.true.label), 
                                 labels=c('Monocyte', 'Bcell', 'CD4TReg', 'CD4NaiveT', 'CD4TCM', 'NK', 'CD8NaiveCTL'))

true.label.14k.reorder <- true.label.14k
true.label.14k.reorder$MICA.cluster <- test.factor.MICA
true.label.14k.reorder$Seurat.cluster <- test.factor.Seurat
true.label.14k.reorder$label <- test.factor.true.label



library(ggsankey)

df <- true.label.14k.reorder %>%
  make_long(MICA.cluster, label, Seurat.cluster)

# install.packages("ggplot2")
library(ggplot2)
# install.packages("dplyr")
library(dplyr) # Also needed

ggplot(df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey(flow.alpha = 0.9, node.color = 'gray60') +
  theme_sankey(base_size = 16)




true.label.20k <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt', 
                             sep = '\t', header = T)
true.label.14k <- true.label.20k[true.label.20k$cell %in% cluster.res$ID,]
pData(eset.14k)$trueLabel <- as.factor(true.label.14k$label)


colors <- c('#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#F8766D', '#C3C3C3', '#C3C3C3')
MICAplot(input_eset = eset.14k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "trueLabel", colors = colors, pct = 0.5)


colors <- c('#C3C3C3', '#C3C3C3', '#CD9600', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3')
MICAplot(input_eset = eset.14k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "trueLabel", colors = colors, pct = 0.5)



true.label.14k$MICA.cluster <- as.factor(cluster.res$label)
true.label.14k$Seurat.cluster <- pbmc@meta.data$seurat_clusters

test.factor.MICA <- factor(true.label.14k$MICA.cluster, levels=c(1,2,3,4,5,6,7,8), 
                           labels=c('M.Monocyte', 'M.CD8NaiveCTL', 'M.Bcell1', 'M.CD4TReg', 'M.NK', 'M.Bcell2',
                                    'M.CD4NaiveT', 'M.CD4TCM'))
test.factor.Seurat <- factor(true.label.14k$Seurat.cluster, levels=c(0,1,2,3,4,5,6,7), 
                             labels=c('S.CD4TCM', 'S.CD4NaiveT', 'S.NK', 'S.Bcell1',  'S.Monocyte1',  'S.CD8NaiveCTL',
                                      'S.Bcell2', 'S.Monocyte2'))
test.factor.true.label <- as.factor(true.label.14k$label)
test.factor.true.label <- factor(test.factor.true.label, levels=levels(test.factor.true.label), 
                                 labels=c('Monocyte', 'Bcell', 'CD4TReg', 'CD4NaiveT', 'CD4TCM', 'NK', 'CD8NaiveCTL'))

true.label.14k.reorder <- true.label.14k
true.label.14k.reorder$MICA.cluster <- test.factor.MICA
true.label.14k.reorder$Seurat.cluster <- test.factor.Seurat
true.label.14k.reorder$label <- test.factor.true.label



# Refined annotations using true labels
true.label.14k.refined <- true.label.14k.reorder
true.label.14k.refined$MICA.cluster <- factor(true.label.14k.refined$MICA.cluster, 
                                              levels=c('M.Monocyte', 'M.CD8NaiveCTL', 'M.Bcell1', 'M.CD4TReg', 'M.NK', 'M.Bcell2',
                                                       'M.CD4NaiveT', 'M.CD4TCM'), 
                                              labels=c('Monocyte', 'CD8NaiveCTL', 'Bcell1', 'CD4TReg',
                                                       'NK', 'Bcell2', 'CD4NaiveT', 'CD4TCM'))
true.label.14k.refined$Seurat.cluster <- factor(true.label.14k.refined$Seurat.cluster, 
                                                levels = c('S.CD4TCM', 'S.CD4NaiveT', 'S.NK', 'S.Bcell1',  'S.Monocyte1',  'S.CD8NaiveCTL',
                                                           'S.Bcell2', 'S.Monocyte2'),
                                                labels = c('CD4TCM', 'CD4NaiveT', 'NK', 'Bcell1', 'Monocyte1', 'CD8NaiveCTL',
                                                           'Bcell2', 'Monocyte2'))



# Calculate cluster cell type percentages using the true labels
Seurat.subcluster <- true.label.14k.refined[true.label.14k.refined$Seurat.cluster == 'Monocyte2',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(Seurat.subcluster$label)),
  count=c(table(Seurat.subcluster$label))
)

# Compute percentages
data$fraction <- data$count / sum(data$count)




# Calculate cluster cell type percentages using the true labels
MICA.subcluster <- true.label.14k.refined[true.label.14k.refined$MICA.cluster == 'CD4TCM',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(MICA.subcluster$label)),
  count=c(table(MICA.subcluster$label))
)

# Compute percentages
data$fraction <- data$count / sum(data$count)






library(ggplot2)

# Calculate cluster cell type percentages using the true labels
subcluster <- true.label.14k.refined[true.label.14k.refined$Seurat.cluster == 'CD4TCM',]

# Donut charts 1
# Create test data.
data <- data.frame(
  category=names(table(subcluster$label)),
  count=c(table(subcluster$label))
)

# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category, "\n fraction: ", data$fraction)

# Make the plot
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=8) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")
