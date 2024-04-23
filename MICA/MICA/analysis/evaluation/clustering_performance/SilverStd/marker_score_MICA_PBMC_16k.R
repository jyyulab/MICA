# MICA clustering and marker score

library(scMINER)
library(NetBID2)


mat.20k <- read.table(file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_MICA_input.txt', 
                      sep = '\t', header = TRUE, row.names = 1)
mat.20k <- as.matrix(mat.20k)


filtered_cells <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/filtered_cells_16k.txt', 
                             sep='\t', header = TRUE, row.names = 1)
mat.16k <- mat.20k[filtered_cells$cell,]


eset.16k <- generate.eset(exp_mat=t(mat.16k),
                          phenotype_info=rownames(mat.16k),
                          feature_info=colnames(mat.16k),
                          annotation_info='exp.log2')



cluster.res <- read.table(file = '/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/MICA/clustering_UMAP_euclidean_24_3.32.txt', 
                          sep = '\t', header = TRUE)

pData(eset.16k)$cellType <- as.factor(cluster.res$label)
pData(eset.16k)$X <- cluster.res$X
pData(eset.16k)$Y <- cluster.res$Y
exp.log2 <- t(mat.16k)





# Heatmap
library(openxlsx)
ref<-read.xlsx("/Users/lding/Git/scMINER/docs/tests/Ref/PBMC_markers.xlsx", sheet="Markers_10_cell_types_PBMC_20k")
head(ref)

celltypes <- c("CD4TCM", "CD4TReg", "CD4NaiveT", "CD8NaiveCTL", "Bcell", "NK", "Monocyte", "HSPC")

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
df_n<-data.frame(label=pData(eset.16k)[,group_name],n_mtx)
df_n<-aggregate(.~label,df_n,mean)
library(reshape2)
df_n_melt<-melt(df_n,id.vars = "label")


df <- data.frame(label=pData(eset.16k)[,group_name],ac_norm);
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



pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/MICA/Heatmap_marker_score_Adam_markers.pdf", 
    height = 8, width=8)

h.mat <- t(input[,-1])
h.mat <- h.mat[,c(4, 7, 5, 2, 6, 1, 3, 8)]

Heatmap(as.matrix(h.mat), col = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)), name = 'MarkerScore',
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,)
dev.off()




