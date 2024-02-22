# Seurat clustering and marker score

library(Seurat)
library(dplyr)
library(patchwork)


pbmc <- readRDS(file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20K_seurat.rds')
filtered_cells <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/filtered_cells_16k.txt', 
                             sep='\t', header = TRUE, row.names = 1)

pbmc <- subset(x = pbmc, cells = names(pbmc$orig.ident)[names(pbmc$orig.ident) %in% filtered_cells$cell])



pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)


pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.2)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")


saveRDS(pbmc, file = "/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/pbmc_20k_seurat_filter_16k.rds")





# Make maker list heatmap plot
# Seurat clustering
exp.count <- as.matrix(pbmc@assays$RNA@counts)

# CPM normalization
norm = 1e5
exp.norm <- sweep(exp.count, 2, norm/unname(Matrix::colSums(exp.count)), '*')

# log transformation
exp.log2 <- log(exp.norm+1, base=2)


library(NetBID2)
pbmc.eset.log2 <- generate.eset(exp_mat=exp.log2, 
                                phenotype_info=pbmc@meta.data,
                                feature_info=pbmc@assays$RNA@counts@Dimnames[[1]],
                                annotation_info='exp.log2')

gene.symbol <- rownames(exprs(pbmc.eset.log2))
fData(pbmc.eset.log2) <- data.frame(gene.symbol)




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
group_name = 'seurat_clusters'
df_n<-data.frame(label=pData(pbmc.eset.log2)[,group_name],n_mtx)
df_n<-aggregate(.~label,df_n,mean)
library(reshape2)
df_n_melt<-melt(df_n,id.vars = "label")


df <- data.frame(label=pData(pbmc.eset.log2)[,group_name],ac_norm);
df <- df[,colSums(is.na(df))<nrow(df)];#remove NA columns
df <- aggregate(.~label,df,mean)
input<-t(apply(df[,-1],1,scale))#row normalization
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



pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_16k/Seurat/Heatmap_marker_score_Adam_markers.pdf", 
    height = 8, width=8)

h.mat <- t(input[,-1])
h.mat <- h.mat[,c(1, 2, 6, 3, 5, 4, 7, 8)]

Heatmap(as.matrix(h.mat), col = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)), name = 'MarkerScore',
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,)
dev.off()


