# MICA clustering and marker score
library(scMINER)
library(NetBID2)


mat.20k <- read.table(file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_MICA_input.txt', 
           sep = '\t', header = TRUE, row.names = 1)
mat.20k <- as.matrix(mat.20k)
eset.20k <- generate.eset(exp_mat=t(mat.20k),
                          phenotype_info=rownames(mat.20k),
                          feature_info=colnames(mat.20k),
                          annotation_info='exp.log2')



cluster.res <- read.table(file = '/Users/lding/Documents/MICA/outputs/PBMC_20k/MICA/clustering_UMAP_euclidean_24_1.8.txt', 
                          sep = '\t', header = TRUE)

pData(eset.20k)$cellType <- as.factor(cluster.res$label)
pData(eset.20k)$X <- cluster.res$X
pData(eset.20k)$Y <- cluster.res$Y
exp.log2 <- t(mat.20k)



true.label.20k <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt', 
                             sep = '\t', header = T)
pData(eset.20k)$trueLabel <- as.factor(true.label.20k$label)


# c(9, 2, 4, 8, 6, 3, 7, 1, 5, 10)
colors <- c('#00B0F6', '#9590FF', '#00BF7D', '#D89000', '#E76BF3', '#39B600', '#00BFC4', '#A3A500', '#F8766D', '#FF62BC')
MICAplot(input_eset = eset.20k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "cellType", colors = colors, pct = 0.5)


colors <- c('#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#F8766D', '#C3C3C3', '#C3C3C3', '#C3C3C3')
MICAplot(input_eset = eset.20k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "trueLabel", colors = colors, pct = 0.5)


colors <- c('#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#CD9600', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3')
MICAplot(input_eset = eset.20k,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "trueLabel", colors = colors, pct = 0.5)





# Read true label file
true.label.20k <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt', 
                             sep = '\t', header = T)
true.label.20k$MICA.cluster <- as.factor(cluster.res$label)
true.label.20k <- true.label.20k[true.label.20k$cell %in% pbmc@assays$RNA@counts@Dimnames[[2]],]
true.label.20k$Seurat.cluster <- pbmc@meta.data$seurat_clusters

# test.factor.MICA <- factor(true.label.20k$MICA.cluster, levels=c(1,2,3,4,5,6,7,8,9,10), labels=c(9, 2, 4, 8, 6, 3, 7, 1, 5, 10))
# test.factor.Seurat <- factor(true.label.20k$Seurat.cluster, levels=c(0,1,2,3,4,5,6,7,8,9), labels=c(1, 8, 2, 3, 4, 5, 6, 7, 9, 10))

test.factor.MICA <- factor(true.label.20k$MICA.cluster, levels=c(1,2,3,4,5,6,7,8,9,10), 
                           labels=c('M.Monocyte', 'M.CD8CTL', 'M.Bcell', 'M.CD4TCM', 'M.HSPC', 'M.CD8NaiveCTL',
                                    'M.NK', 'M.CD4NaiveT', 'M.CD4HeloperT', 'M.DC'))
test.factor.Seurat <- factor(true.label.20k$Seurat.cluster, levels=c(0,1,2,3,4,5,6,7,8,9), 
                             labels=c('S.CD4HelperT/CD4NaiveT', 'S.CD8NaiveCTL', 'S.CD4TCM', 'S.CD4HelperT/CD4TCM',  'S.NK', 'S.Bcell',
                                      'S.Monocyte1', 'S.HSPC1', 'S.HSPC2', 'S.Monocyte2'))
test.factor.true.label <- as.factor(true.label.20k$label)
test.factor.true.label <- factor(test.factor.true.label, levels=levels(test.factor.true.label), 
                                 labels=c('Monocyte', 'Bcell', 'HSPC', 'CD4HelperT', 'CD4TReg', 'CD4NaiveT', 'CD4TCM', 'NK', 'CD8CTL', 'CD8NaiveCTL'))

true.label.20k.reorder <- true.label.20k
true.label.20k.reorder$MICA.cluster <- test.factor.MICA
true.label.20k.reorder$Seurat.cluster <- test.factor.Seurat
true.label.20k.reorder$label <- test.factor.true.label





library(ggsankey)

df <- true.label.20k.reorder %>%
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








# Refined annotations using true labels
true.label.20k.refined <- true.label.20k.reorder
true.label.20k.refined$MICA.cluster <- factor(true.label.20k.refined$MICA.cluster, 
                                              levels=c('M.Monocyte', 'M.CD8CTL', 'M.Bcell', 'M.CD4TCM', 'M.HSPC', 'M.CD8NaiveCTL',
                                                       'M.NK', 'M.CD4NaiveT', 'M.CD4HeloperT', 'M.DC'), 
                                              labels=c('Monocyte', 'CD8CTL', 'Bcell', 'CD4TReg.CD4Helper',
                                                       'HSPC1', 'CD8NaiveCTL.CD8CTL', 'NK', 'CD4NaiveT.CD4HelperT.CD4TReg', 
                                                       'CD4TCM', 'HSPC2'))

true.label.20k.refined$Seurat.cluster <- factor(true.label.20k.refined$Seurat.cluster, 
                                                levels = c('S.CD4HelperT/CD4NaiveT', 'S.CD8NaiveCTL', 'S.CD4TCM', 'S.CD4HelperT/CD4TCM',  
                                                           'S.NK', 'S.Bcell', 'S.Monocyte1', 'S.HSPC1', 'S.HSPC2', 'S.Monocyte2'),
                                                labels = c('CD4HelperT.CD4NaiveT', 'CD8NaiveCTL', 
                                                           'CD4TCM', 'CD4HelperT.CD4TCM', 'NK', 'Bcell', 'Monocyte1', 
                                                           'HSPC1', 'HSPC2', 'Monocyte2'))




library(ggplot2)

# Calculate cluster cell type percentages using the true labels
Seurat.cluster <- true.label.20k.refined[true.label.20k.refined$Seurat.cluster == 'Monocyte2',]

# Donut charts 1
# Create test data.
data <- data.frame(
  category=names(table(Seurat.cluster$label)),
  count=c(table(Seurat.cluster$label))
)

# Compute percentages
data$fraction <- data$count / sum(data$count)







# Calculate cluster cell type percentages using the true labels
MICA.cluster.CD4TCM <- true.label.20k.refined[true.label.20k.refined$MICA.cluster == 'CD4TCM',]

# Donut charts 1
# Create test data.
data <- data.frame(
  category=names(table(MICA.cluster.CD4TCM$label)),
  count=c(table(MICA.cluster.CD4TCM$label))
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
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")




# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(MICA.cluster.CD4TCM$label)),
  count=c(table(MICA.cluster.CD4TCM$label))
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
  geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
  scale_fill_brewer(palette=5) +
  scale_color_brewer(palette=5) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")





# Calculate cluster cell type percentages using the true labels
Seurat.cluster.CD4TCM <- true.label.20k.refined[true.label.20k.refined$Seurat.cluster == 'CD4HelperT.CD4TCM',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(Seurat.cluster.CD4TCM$label)),
  count=c(table(Seurat.cluster.CD4TCM$label))
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
  geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
  scale_fill_brewer(palette=3) +
  scale_color_brewer(palette=3) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")





# Calculate cluster cell type percentages using the true labels
MICA.cluster.CD4TReg <- true.label.20k.refined[true.label.20k.refined$MICA.cluster == 'CD4TReg.CD4Helper',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(MICA.cluster.CD4TReg$label)),
  count=c(table(MICA.cluster.CD4TReg$label))
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
  geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
  scale_fill_brewer(palette=2) +
  scale_color_brewer(palette=2) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")






# Calculate cluster cell type percentages using the true labels
MICA.cluster.HSPC <- true.label.20k.refined[true.label.20k.refined$MICA.cluster == 'HSPC',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(MICA.cluster.HSPC$label)),
  count=c(table(MICA.cluster.HSPC$label))
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
  geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
  scale_fill_brewer(palette=4) +
  scale_color_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")





# Calculate cluster cell type percentages using the true labels
Seurat.cluster.HSPC1 <- true.label.20k.refined[true.label.20k.refined$Seurat.cluster == 'HSPC1',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(Seurat.cluster.HSPC1$label)),
  count=c(table(Seurat.cluster.HSPC1$label))
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
  geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
  scale_fill_brewer(palette=11) +
  scale_color_brewer(palette=11) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")





# Calculate cluster cell type percentages using the true labels
Seurat.cluster.HSPC2 <- true.label.20k.refined[true.label.20k.refined$Seurat.cluster == 'HSPC2',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(Seurat.cluster.HSPC2$label)),
  count=c(table(Seurat.cluster.HSPC2$label))
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
  geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
  scale_fill_brewer(palette=13) +
  scale_color_brewer(palette=13) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")





# Calculate cluster cell type percentages using the true labels
Seurat.cluster.Monocyte1 <- true.label.20k.refined[true.label.20k.refined$Seurat.cluster == 'Monocyte1',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(Seurat.cluster.Monocyte1$label)),
  count=c(table(Seurat.cluster.Monocyte1$label))
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
  geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
  scale_fill_brewer(palette=11) +
  scale_color_brewer(palette=11) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")






# Calculate cluster cell type percentages using the true labels
Seurat.cluster.Monocyte2 <- true.label.20k.refined[true.label.20k.refined$Seurat.cluster == 'Monocyte2',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(Seurat.cluster.Monocyte2$label)),
  count=c(table(Seurat.cluster.Monocyte2$label))
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
  geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
  scale_fill_brewer(palette=13) +
  scale_color_brewer(palette=13) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")




# Calculate cluster cell type percentages using the true labels
Seurat.cluster.HSPC2 <- true.label.20k.refined[true.label.20k.refined$Seurat.cluster == 'HSPC2',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(Seurat.cluster.HSPC2$label)),
  count=c(table(Seurat.cluster.HSPC2$label))
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
  geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
  scale_fill_brewer(palette=13) +
  scale_color_brewer(palette=13) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")






# Calculate cluster cell type percentages using the true labels
MICA.subcluster<- true.label.20k.refined[true.label.20k.refined$MICA.cluster == 'HSPC2',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(MICA.subcluster$label)),
  count=c(table(MICA.subcluster$label))
)

# Compute percentages
data$fraction <- data$count / sum(data$count)






# Calculate cluster cell type percentages using the true labels
MICA.subcluster<- true.label.20k.refined[true.label.20k.refined$MICA.cluster == 'CD4NaiveT.CD4HelperT.CD4TReg',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(MICA.subcluster$label)),
  count=c(table(MICA.subcluster$label))
)

# Compute percentages
data$fraction <- data$count / sum(data$count)






# Calculate cluster cell type percentages using the true labels
MICA.subcluster<- true.label.20k.refined[true.label.20k.refined$MICA.cluster == 'HSPC',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(MICA.subcluster$label)),
  count=c(table(MICA.subcluster$label))
)

# Compute percentages
data$fraction <- data$count / sum(data$count)






# Calculate cluster cell type percentages using the true labels
Seurat.subcluster <- true.label.20k.refined[true.label.20k.refined$Seurat.cluster == 'CD8CTL',]


# Donut charts 2
# Create test data.
data <- data.frame(
  category=names(table(Seurat.subcluster$label)),
  count=c(table(Seurat.subcluster$label))
)

# Compute percentages
data$fraction <- data$count / sum(data$count)









library(openxlsx)
ref<-read.xlsx("/Users/lding/Git/scMINER/docs/tests/Ref/PBMC_markers.xlsx", sheet="Markers_10_cell_types_PBMC_20k")
head(ref)

celltypes <- c( "CD4HelperT",  "CD4TCM", "CD4TReg", "CD4NaiveT","CD8CTL", "CD8NaiveCTL", "Bcell", "NK", "Monocyte", "HSPC")

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
df_n<-data.frame(label=pData(eset.20k)[,group_name],n_mtx)
df_n<-aggregate(.~label,df_n,mean)
library(reshape2)
df_n_melt<-melt(df_n,id.vars = "label")


df <- data.frame(label=pData(eset.20k)[,group_name],ac_norm);
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

# p<-draw.bubblePlot2(df=d, xlab="Cluster",ylab="CellType",
#                     clab="MarkerScore",slab="ExpressionPercentage",
#                     plot.title="Cell type annotation for each cluster")






pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_3/PBMC_20k/MICA/Heatmap_marker_score_Adam_markers_24_1.8_2.pdf", 
    height = 10, width=10)

h.mat <- t(input[,-1])
h.mat <- h.mat[,c(9, 4, 8, 2, 6, 3, 7, 1, 5, 10)]

Heatmap(as.matrix(h.mat), col = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)), name = 'MarkerScore',
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,)
dev.off()



pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_20k/MICA/Heatmap_marker_score_Adam_markers_24_2.2.pdf", 
    height = 10, width=10)

h.mat <- t(input[,-1])
h.mat <- h.mat[,c(7, 5, 4, 6, 10, 3, 8, 1, 9, 2)]

Heatmap(as.matrix(h.mat), col = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)), name = 'MarkerScore',
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,)
dev.off()



pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_20k/MICA/Heatmap_marker_score_Adam_markers_36_1.8.pdf", 
    height = 10, width=10)

h.mat <- t(input[,-1])
h.mat <- h.mat[,c(8, 5, 4, 7, 10, 3, 6, 1, 2, 9)]

Heatmap(as.matrix(h.mat), col = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)), name = 'MarkerScore',
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,)
dev.off()




pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_20k/MICA/Heatmap_marker_score_Adam_markers_44_1.8.pdf", 
    height = 10, width=10)

h.mat <- t(input[,-1])
h.mat <- h.mat[,c(5, 6, 4, 7, 8, 3, 9, 1, 10, 2)]

Heatmap(as.matrix(h.mat), col = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)), name = 'MarkerScore',
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,)
dev.off()




feature_highlighting<-function(input_eset,target=NULL,
                               feature="geneSymbol",
                               x="X",y="Y",
                               wrap_by=NULL,
                               ylabel="Expression",pct.size=0.8,
                               title.size=15,ncol=4, alpha=0.8,
                               colors=colorRampPalette(c("#E3E3E3", "#BCA2FC","#4900FE"),interpolate="linear")(8)){
  
  # change it to expr is ok
  input<-as.matrix(exprs(input_eset))
  indx<-which(fData(input_eset)[,feature]%in%target)
  if(length(indx)==0) stop("Target feature not found.")
  
  gn<-fData(input_eset)[,feature][indx]
  id.vars<-c(x,y,wrap_by)
  projection<-pData(input_eset)[colnames(input),id.vars]
  
  #gene expression visualized as columns
  if (length(indx)!=1) {
    target_values <- t(as.matrix(input[indx,]))
    colnames(target_values)<-gn
    
    proj_target <- cbind(projection,target_values)
    proj_target_melt <- reshape2::melt(proj_target, id.vars=id.vars)
    
    p<- ggplot(proj_target_melt, aes_string(x, y)) +
      theme_classic()+
      facet_wrap(c("variable",wrap_by),scales = "free",ncol = ncol)
    labs(title="")
    
  }else{
    target_values <- input[indx,]
    proj_target <- cbind(projection,target=target_values)
    proj_target_melt <- reshape2::melt(proj_target, id.vars=id.vars)
    
    p<- ggplot(proj_target_melt, aes_string(x, y)) +
      theme_classic()+
      labs(title=target,scales = "free")
    
    if(!is.null(wrap_by)) p <- p + facet_wrap(c(wrap_by),scales = "free",ncol = ncol)
  }#indx = 1
  
  p<- p + geom_point(aes(colour=value), size = pct.size, stroke = 0, shape = 16, alpha=alpha) +
    scale_colour_gradientn(colors=colors)   +
    theme(plot.title = element_text(size = title.size, face = "bold"),
          axis.title = element_text(size = 10),
          legend.title = element_text(size = 10))+
    labs(color=ylabel)
  
  return(p)
}



eset <- readMICAoutput(eset = eset.20k, load_ClusterRes = TRUE, 
                       output_file = "/Users/lding/Documents/MICA/outputs/PBMC_20k/MICA/clustering_UMAP_euclidean_24_1.8.txt")

# UMAP plot of selected markers
library(ggplot2)
gn.sel<-c("CD4", "IL2RA", "AQP3", "IL32")
p <- feature_highlighting(input_eset = eset, target = gn.sel, 
                          feature="gene", ylabel = "log2Exp", 
                          x="X", y="Y", pct.size = 0.25)

filename <- '/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_20k/MICA/feature_highlighting_selected_markers_CD4TCM.pdf'
ggsave(plot = p, filename = filename, units="in",
       width = 8, height = 2.5, dpi = 300, useDingbats = F)




# HSPC
gn.sel<-c("CD34", "MYB", "KIT", "CD164")
p <- feature_highlighting(input_eset = eset, target = gn.sel, 
                          feature="gene", ylabel = "log2Exp", 
                          x="X", y="Y", pct.size = 0.25)

filename <- '/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_20k/MICA/feature_highlighting_selected_markers_HSPC.pdf'
ggsave(plot = p, filename = filename, units="in",
       width = 8, height = 2.5, dpi = 300, useDingbats = F)





# Heatmap of selected markers
gn.sel<-c("CD3D","CD27","IL7R","SELL","CCR7","IL32","GZMA",
          "GZMK","DUSP2","CD8A","GZMH","GZMB","CD79A","CD79B","CD86","CD14")



#' @title Visualize gene expression level on scRNA-seq data via heatmap
#' @description This plot will visualiz feature info in scatter plot by outputing a ggplot object
#' @param input_eset Input expression set
#' @param feature a character, which feature to visualize
#' @param target a character or a character vector indicating feature names
#' @param group_name a character, label to visualize on the top of heatmap
#' @param name character, name of value visualized in color scale
#' @param cluster_rows logical, if or not cluster rows
#' @param colors color palette
#' @param plot_name character, name of heamap
#' @param save_plot logical, whether to save plots or not
#' @param width numerical
#' @param height numerical
#' @param ... parameter to be passed to ComplexHeatmap::Heatmap
#' @return a ggplot object
#'
#' @export
feature_heatmap <- function(input_eset,target,feature="geneSymbol",
                            group_name="label",name="log2Exp",
                            save_plot=TRUE,width=4,height=8,
                            cluster_rows = FALSE,
                            colors = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)),
                            plot_name="GeneHeatmap.png",
                            ...){
  
  input <- exprs(input_eset)
  gn<-intersect(target,fData(input_eset)[,feature])
  indx<-match(gn,fData(input_eset)[,feature])
  
  exp<-exprs(input_eset)[indx,]
  rownames(exp)<-gn
  lab<-pData(input_eset)[,group_name];names(lab) <- sampleNames(input_eset)
  
  #re-order expressionmatrix and label
  order <- c(5,2,3,8,4,7,1,6)
  ranks <- lab[order(match(lab, order))]
  ranks <- names(ranks)
  exp.ordered<-as.matrix(exp[,ranks])
  lab.ordered<-lab[ranks]
  df<-data.frame(scMINER=lab.ordered)
  
  #Define color annotations
  n<-length(unique(lab.ordered))
  ncols <- scales::hue_pal()(n)
  names(ncols) <- unique(lab.ordered)
  myanndf = HeatmapAnnotation(df = df,col=list(scMINER = ncols))
  mycolors = colors
  
  hmp <- Heatmap(exp.ordered, col = mycolors, name = name,
                 show_row_names = TRUE,
                 show_column_names = FALSE,
                 cluster_rows = cluster_rows,
                 cluster_columns = FALSE,
                 top_annotation = myanndf,
                 ...)
  
  if(save_plot){
    pdf(file = plot_name, width=width, height=height)
    ComplexHeatmap::draw(hmp)
    dev.off()
  }
  
  return(hmp)
}



# Heatmap
feature_heatmap(input_eset = eset.log2.v3, feature="geneSymbol", target = gn.sel, group_name = "cellType",
                save_plot = TRUE, width = 10, height = 5,
                name = "log2Exp", plot_name="/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_12k/MICA/Heatmap_Tracy_markers_log2.pdf")





gn.sel <- c('CTSS',
            'FCN1',
            'NEAT1',
            'LYZ',
            'PSAP',
            'S100A9',
            'AIF1',
            'MNDA',
            'SERPINA1',
            'TYROBP',
            'NKG7',
            'KLRD1',
            'TYROBP',
            'GNLY',
            'FCER1G',
            'PRF1',
            'CD247',
            'KLRF1',
            'CST7',
            'GZMB',
            'CD79A',
            'RALGPS2',
            'CD79B',
            'MS4A1',
            'BANK1',
            'CD74',
            'HLA-DQA1',
            'MEF2C',
            'CD74',
            'HLA-DPA1',
            'HLA-DPB1',
            'HLA-DQA1',
            'CCDC88A',
            'HLA-DRA',
            'HLA-DMA',
            'CST3',
            'HLA-DQB1',
            'HLA-DRB1',
            'TCF7',
            'CD4',
            'CCR7',
            'IL7R',
            'FHIT',
            'LEF1',
            'MAL',
            'NOSIP',
            'LDHB',
            'PIK3IP1',
            'CD8B',
            'S100B',
            'CCR7',
            'RGS10',
            'NOSIP',
            'LEF1',
            'CRTAM',
            'CD8A',
            'OXNAD1',
            'IL7R',
            'TMSB10',
            'CD4',
            'ITGB1',
            'LTB',
            'AQP3',
            'LDHB',
            'IL32',
            'MAL',
            'CD8B',
            'ANXA1',
            'CD8A',
            'KRT1',
            'YBX3',
            'IL7R',
            'NELL2',
            'LDHB',
            'IL7R',
            'CCL5',
            'GZMK',
            'IL32',
            'GZMA',
            'KLRB1',
            'LTB',
            'AQP3',
            'CD8B',
            'ANXA1',
            'CD8A',
            'KRT1',
            'YBX3',
            'IL7R',
            'NELL2',
            'LDHB')

# c("NaiveT", "CD4TCM", "CD8TCM", "TEM", "Monocyte", "NK", "Bcell", "DC")
# c( 5,        2,        3,        8,     4,          7,    1,       6)
gn.sel <- c('CCR7', 'MAL', 'IL7R', 'AQP3', 'CD8B', 'NELL2', 'IL32', 'GZMK', 'KLRB1', 'CTSS', 'FCN1', 'SERPINA1', 
            'NKG7', 'GNLY', 'PRF1', 'GZMB', 'CD79A', 'CD79B', 'MS4A1', 'CST3', 'HLA-DRB1')

pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_12k/MICA/feature_vlnplot_8_cell_types_21_markers.pdf", 
    height = 50, width=50)
feature_vlnplot(input_eset=eset.log2.v3, target=gn.sel, feature = "geneSymbol",
                group_by = "cellType", ylabel = "log2Exp", ncol = 4)
dev.off()




# Heatmap
feature_heatmap(input_eset = eset.log2.v3, feature="geneSymbol", target = gn.sel, group_name = "cellType",
                save_plot = TRUE, width = 10, height = 5,
                name = "log2Exp", plot_name="/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_12k/MICA/Heatmap_21_markers_log2.pdf")

