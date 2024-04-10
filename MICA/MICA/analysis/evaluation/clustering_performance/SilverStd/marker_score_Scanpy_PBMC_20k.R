# MICA clustering and marker score
library(scMINER)
library(NetBID2)


mat.20k <- read.table(file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_MICA_input.txt', 
                      sep = '\t', header = TRUE, row.names = 1)
mat.20k <- as.matrix(mat.20k)


cluster.res <- read.table(file = '/Users/lding/Documents/MICA/outputs/PBMC_20k/Scanpy/PBMC_20k_predict_label.txt', 
                          sep = '\t', header = TRUE)
mat.20k <- mat.20k[cluster.res$Cell,]


eset.20k <- generate.eset(exp_mat=t(mat.20k),
                          phenotype_info=rownames(mat.20k),
                          feature_info=colnames(mat.20k),
                          annotation_info='exp.log2')


pData(eset.20k)$cellType <- as.factor(cluster.res$leiden)
# pData(eset.20k)$X <- cluster.res$X
# pData(eset.20k)$Y <- cluster.res$Y
exp.log2 <- t(mat.20k)








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






pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_3/PBMC_20k/Scanpy/Heatmap_marker_score_Adam_markers_orig.pdf", 
    height = 10, width=10)

h.mat <- t(input[,-1])
# h.mat <- h.mat[,c(1, 2, 8, 3, 4, 6, 5, 10, 7, 9)]

Heatmap(as.matrix(h.mat), col = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)), name = 'MarkerScore',
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,)
dev.off()









cluster.res.MICA <- read.table(file = '/Users/lding/Documents/MICA/outputs/PBMC_20k/MICA/clustering_UMAP_euclidean_24_1.8.txt', 
                               sep = '\t', header = TRUE)

# Read true label file
true.label.20k <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt', 
                             sep = '\t', header = T)
true.label.20k$MICA.cluster <- as.factor(cluster.res.MICA$label)
true.label.20k <- true.label.20k[true.label.20k$cell %in% cluster.res$Cell,]
true.label.20k$Scanpy.cluster <- as.factor(cluster.res$leiden)






# test.factor.MICA <- factor(true.label.20k$MICA.cluster, levels=c(1,2,3,4,5,6,7,8,9,10), labels=c(9, 2, 4, 8, 6, 3, 7, 1, 5, 10))
# test.factor.Seurat <- factor(true.label.20k$Seurat.cluster, levels=c(0,1,2,3,4,5,6,7,8,9), labels=c(1, 8, 2, 3, 4, 5, 6, 7, 9, 10))

test.factor.MICA <- factor(true.label.20k$MICA.cluster, levels=c(1,2,3,4,5,6,7,8,9,10), 
                           labels=c('M.Monocyte', 'M.CD8CTL', 'M.Bcell', 'M.CD4TCM', 'M.HSPC1', 'M.CD8NaiveCTL',
                                    'M.NK', 'M.CD4NaiveT', 'M.CD4HelperT', 'M.HSPC2'))
# c(1, 2, 8, 3, 4, 6, 5, 10, 7, 9)
test.factor.Scanpy <- factor(true.label.20k$Scanpy.cluster, levels=c(4,7,3,1,6,2,5,0,8,9), 
                             labels=c('S.Monocyte1', 'S.CD8CTL', 'S.Bcell', 'S.CD4TCM', 'S.HSPC1', 'S.CD8NaiveCTL',
                                      'S.NK', 'S.CD4NaiveT', 'S.HSPC2', 'S.Monocyte2'))
test.factor.true.label <- as.factor(true.label.20k$label)
test.factor.true.label <- factor(test.factor.true.label, levels=levels(test.factor.true.label), 
                                 labels=c('Monocyte', 'Bcell', 'HSPC', 'CD4HelperT', 'CD4TReg', 'CD4NaiveT', 'CD4TCM', 'NK', 'CD8CTL', 'CD8NaiveCTL'))

true.label.20k.reorder <- true.label.20k
true.label.20k.reorder$MICA.cluster <- test.factor.MICA
true.label.20k.reorder$Scanpy.cluster <- test.factor.Scanpy
true.label.20k.reorder$label <- test.factor.true.label




library(ggsankey)
df <- true.label.20k.reorder %>%
  make_long(MICA.cluster, label, Scanpy.cluster)

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
