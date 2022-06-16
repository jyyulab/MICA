# Seurat clustering and marker score

library(Seurat)
library(dplyr)
library(patchwork)




# Create Seurat object
raw_counts<-read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/Filtered_DownSampled_SortedPBMC_data.csv', 
                       sep=",", row.names = 1, header = TRUE)
head(raw_counts)
raw_counts_t <- t(raw_counts)
pbmc <- CreateSeuratObject(counts = raw_counts_t, min.cells = 3, min.genes = 200, project = "pbmc")



# Filter cells
filtered_cells <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/filtered_cells_14k.txt', 
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
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 4000)
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 6000)
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 8000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# pbmc <- RunPCA(pbmc, features = rownames(pbmc@assays$RNA@counts))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)


pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.27)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")


saveRDS(pbmc, file = "/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/pbmc_20k_seurat_filter_to_14k.rds")






# find all markers of Seurat cluster 0
cluster2.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.05, group.by = 'seurat_clusters')
cluster2.markers[c('CD4','IL2RA','FOXP3','IL2RA','TIGIT','TNFRSF18','ENTPD1','CTLA4','LAIR2'),]
head(cluster2.markers, n = 10)
VlnPlot(pbmc, features = c("FOXP3", "IL2RA"))
write.table(cluster2.markers, file='/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3_1/PBMC_14k/Seurat/filter_cytoT/cluster0.markers.txt', sep = '\t')




# find all markers of Seurat cluster 4
cluster.res <- read.table(file = '/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3_1/PBMC_14k/MICA/filter_CTL/clustering_UMAP_euclidean_24_2.72.txt', 
                          sep = '\t', header = TRUE)
cluster.res.factor <- as.factor(cluster.res$label)
pbmc@meta.data$MICA_clusters <- cluster.res.factor
VlnPlot(pbmc, features = c("FOXP3", "IL2RA"), group.by = 'MICA_clusters')
cluster4.markers.MICA <- FindMarkers(pbmc, ident.1 = 4, group.by = 'MICA_clusters', min.pct = 0.05)
cluster4.markers.MICA[c('CD4','IL2RA','FOXP3','IL2RA','TIGIT','TNFRSF18','ENTPD1','CTLA4','LAIR2'),]

cluster4.markers.MICA <- FindMarkers(pbmc, ident.1 = 4, group.by = 'MICA_clusters', min.pct = 0.05)
head(cluster4.markers.MICA, n = 10)
cluster4.markers.MICA[c('CD4','IL2RA','FOXP3','IL2RA','TIGIT','TNFRSF18','ENTPD1','CTLA4','LAIR2'),]
write.table(cluster4.markers.MICA, 
            file='/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3_1/PBMC_14k/MICA/filter_CTL/cluster4.markers.txt', 
            sep = '\t')





library(EnhancedVolcano)
EnhancedVolcano(cluster2.markers,
                lab = rownames(cluster2.markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-50,
                FCcutoff = 0.5,
                selectLab = c('IL2RA', 'FOXP3', 'TIGIT', 'TNFRSF18',
                              'RTKN2', 'AC133644.2', 'CD4', 'CTLA4', 'FCRL3', 'LAIR2', 'IKZF2', # Azimuth Treg 
                              'B2M', 'S100A4', 'FOXP3', 'IKZF2', 'TRAC',  # Azimuth Treg Memory
                              'LEF1', 'C12orf57', 'TOMM7', 'CCR7', 'LDHB'),  # Azimuth Treg Naive
                drawConnectors = TRUE,
                labFace = 'bold')


EnhancedVolcano(cluster4.markers.MICA,
                lab = rownames(cluster4.markers.MICA),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-50,
                FCcutoff = 0.5,
                selectLab = c('IL2RA', 'FOXP3', 'TIGIT', 'TNFRSF18', 'LRRC32',
                              'RTKN2', 'AC133644.2', 'CD4', 'CTLA4', 'FCRL3', 'LAIR2', 'IKZF2', # Azimuth Treg 
                              'B2M', 'S100A4', 'FOXP3', 'IKZF2', 'TRAC',  # Azimuth Treg Memory
                              'LEF1', 'C12orf57', 'TOMM7', 'CCR7', 'LDHB'),  # Azimuth Treg Naive
                drawConnectors = TRUE,
                labFace = 'bold')





# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells.
FOXP3.freq <- as.data.frame( table(exp.count['FOXP3', rownames(pbmc@meta.data[pbmc@meta.data$MICA_clusters == 4,])]) )
colnames(FOXP3.freq) <- c('category', 'Freq.MICA')
FOXP3.freq$Freq.Seurat <- c(table(exp.count['FOXP3', rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters == 1,])]))
FOXP3.freq$Freq.MICA <- FOXP3.freq$Freq.MICA / sum(FOXP3.freq$Freq.MICA)
FOXP3.freq$Freq.Seurat <- FOXP3.freq$Freq.Seurat / sum(FOXP3.freq$Freq.Seurat)

mFOXP3.freq <- melt(FOXP3.freq, id=c('category'))
mFOXP3.freq <- mFOXP3.freq[c(2,3,5,6),]
mFOXP3.freq$label_ypos <- c(0.05, 0.006, 0.009, 0.0003)



# Create the barplot
ggplot(data=mFOXP3.freq, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=value), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_brewer(palette='Set2')+
  theme_minimal()






# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells.
IL2RA.freq <- as.data.frame( table(exp.count['IL2RA', rownames(pbmc@meta.data[pbmc@meta.data$MICA_clusters == 4,])]) )
colnames(IL2RA.freq) <- c('category', 'Freq.MICA')
IL2RA.freq$Freq.Seurat <- c(table(exp.count['IL2RA', rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters == 1,])]), 0, 0, 0)
IL2RA.freq$Freq.MICA <- IL2RA.freq$Freq.MICA / sum(IL2RA.freq$Freq.MICA)
IL2RA.freq$Freq.Seurat <- IL2RA.freq$Freq.Seurat / sum(IL2RA.freq$Freq.Seurat)

mIL2RA.freq <- melt(IL2RA.freq, id=c('category'))
mIL2RA.freq <- mIL2RA.freq[c(2,3,4,5,7),]
mIL2RA.freq$label_ypos <- c(0.06, 0.008, 0.0007, 0.0007, 0.006)

# Create the barplot
ggplot(data=mIL2RA.freq, aes(x=variable, y=value, fill=category)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=value), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_brewer(palette='Set2')+
  theme_minimal()






# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells.

TIGIT <- table(exp.count['TIGIT', rownames(pbmc@meta.data[pbmc@meta.data$MICA_clusters == 4,])])
# TIGIT <- c(TIGIT, "7"=0)
# TIGIT <- c(TIGIT, "8"=0)
TIGIT.freq <- as.data.frame( TIGIT )
# TIGIT.freq$category <- as.integer( rownames(TIGIT.freq) )
colnames(TIGIT.freq) <- c('category', 'Freq.MICA')
TIGIT.freq$Freq.Seurat <- as.vector( table(exp.count['TIGIT', rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters == 0,])]) )
TIGIT.freq$Freq.MICA <- TIGIT.freq$Freq.MICA / sum(TIGIT.freq$Freq.MICA)
TIGIT.freq$Freq.Seurat <- TIGIT.freq$Freq.Seurat / sum(TIGIT.freq$Freq.Seurat)

table(exp.count['TIGIT', rownames(pbmc@meta.data[pbmc@meta.data$MICA_clusters == 4,])])

mTIGIT.freq <- melt(TIGIT.freq, id=c('category'))
mTIGIT.freq <- mTIGIT.freq[c(2,3,4,5,6,8,9,10,11,12),]
mTIGIT.freq$label_ypos <- c(0.05, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001)



# Create the barplot
ggplot(data=mTIGIT.freq, aes(x=variable, y=value, fill=as.factor(category))) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=value), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_brewer(palette='Set2')+
  theme_minimal()









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





pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3_1/PBMC_14k/Seurat/Heatmap_marker_score_Adam_markers4_reorder.pdf", 
    height=10, width=10)

h.mat <- t(input[,-1])
h.mat <- h.mat[,c(1, 2, 6, 4, 7, 3, 5, 8)]


Heatmap(as.matrix(h.mat), col = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)), name = 'MarkerScore',
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,)
dev.off()


  

true.label.20k <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt', 
                             sep = '\t', header = T)
true.label.14k <- true.label.20k[true.label.20k$cell %in% cluster.res$ID,]

DimPlot(pbmc, cells.highlight = true.label.14k[true.label.14k$label == 'CD4+/CD45RO+ Memory',], 
        cols.highlight = c('#F8766D'), reduction = "umap", pt.size = 0.01, sizes.highlight = 0.01)


DimPlot(pbmc, cells.highlight = true.label.14k[true.label.14k$label == 'CD4+/CD25 T Reg',], 
        cols.highlight = c('#CD9600'), reduction = "umap", pt.size = 0.01, sizes.highlight = 0.01)



