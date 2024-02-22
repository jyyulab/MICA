# Seurat clustering and marker score

library(Seurat)
library(dplyr)
library(patchwork)


# Seurat
pbmc.12k <- Read10X(data.dir = "/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_12k/PBMC12k_input")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.12k, project = "pbmc12k", min.cells = 3, min.features = 200)
pbmc


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
pbmc <- FindClusters(pbmc, resolution = 0.12)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")


saveRDS(pbmc, file = "/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_12k/pbmc_12k.rds")



pbmc.markers <- FindAllMarkers(pbmc, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
markers_plus <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 40, order_by = avg_log2FC)

markers_minus <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_min(n = 40, order_by = avg_log2FC)


# Plot markers from reference dataset https://azimuth.hubmapconsortium.org/references/#Human%20-%20PBMC
# CD4+ Naive T Markers
VlnPlot(pbmc, features = c("TCF7", "CD4", "CCR7", "IL7R", "FHIT", "LEF1", "MAL", "NOSIP", "LDHB", "PIK3IP1"))
# CD8+ Naive T Markers
VlnPlot(pbmc, features = c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
# CD4+ Effector Memory T
VlnPlot(pbmc, features = c("IL7R", "CCL5", "FYB1", "GZMK", "IL32", "GZMA", "KLRB1", "TRAC", "LTB", "AQP3"))
# CD8+ Effector Memory T
VlnPlot(pbmc, features = c("CCL5", "GZMH", "CD8A", "TRAC", "KLRD1", "NKG7", "GZMK", "CST7", "CD8B", "TRGC2"))
# CD4+ Central Memory T
VlnPlot(pbmc, features = c("IL7R", "TMSB10", "CD4", "ITGB1", "LTB", "TRAC", "AQP3", "LDHB", "IL32", "MAL"))
# CD8+ Central Memory T
VlnPlot(pbmc, features = c("CD8B", "ANXA1", "CD8A", "KRT1", "LINC02446", "YBX3", "IL7R", "TRAC", "NELL2", "LDHB"))
# NK
VlnPlot(pbmc, features = c("NKG7", "KLRD1", "TYROBP", "GNLY", "FCER1G", "PRF1", "CD247", "KLRF1", "CST7", "GZMB"))
# B cell
VlnPlot(pbmc, features = c("CD79A", "RALGPS2", "CD79B", "MS4A1", "BANK1", "CD74", "TNFRSF13C", "HLA-DQA1", "IGHM", "MEF2C"))
# DC
VlnPlot(pbmc, features = c("CD74", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CCDC88A", "HLA-DRA", "HLA-DMA", "CST3", "HLA-DQB1", 
                           "HLA-DRB1"))
# Monocyte
VlnPlot(pbmc, features = c("CTSS", "FCN1", "NEAT1", "LYZ", "PSAP", "S100A9", "AIF1", "MNDA", "SERPINA1", "TYROBP"))
# pDC
VlnPlot(pbmc, features = c("CCDC50", "UGCG", "TCF4", "LILRA4", "IRF8", "IL3RA", "PLD4", "IRF7", "SERPINF1", "ITM2C"))
# Platelet
VlnPlot(pbmc, features = c("GNG11", "PPBP", "NRGN", "PF4", "CAVIN2", "TUBB1", "HIST1H2AC", "PRKAR2B", "CLU", "F13A1"))




new.cluster.ids <- c("Naive T (CCR7+,LEF1+)", "CD4+ Memory T (AQP3+,IL7R+,IL32+)", "CD8+ Memory T (GZMK+,CD8B+)", 
                     "CD8+ Effector Memory T (CD8A+,GZMH+,CD8B+)", "Monocyte/DC (CTSS+,LYZ+,CST3+)", 
                     "NK (KLRD1+,GNLY+)", "B cell (CD79A+, CD79B+)", "pDC (CCDC50+, IRF7+)", "Platelet (GNG11+, PPBP+)")

new.cluster.ids <- c("Naive T (CCR7+,LEF1+)", "Memory T (AQP3+,IL7R+,IL32+)", 
                     "CD8+ Effector Memory T (CD8A+,GZMH+,CD8B+)", "Monocyte/DC (CTSS+,LYZ+,CST3+)", 
                     "NK (KLRD1+,GNLY+)", "B cell (CD79A+, CD79B+)", "pDC (CCDC50+, IRF7+)", "Platelet (GNG11+, PPBP+)")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()









# Make maker list heatmap plot
# Seurat clustering
exp.count <- as.matrix(pbmc@assays$RNA@counts)

# CPM normalization
norm = 1e5
exp.norm <- sweep(exp.count, 2, norm/unname(Matrix::colSums(exp.count)), '*')

# log transformation
exp.log2 <- log(exp.norm+1, base=2)





# Violin plot, Tracy's marker genes
library(Seurat)
library(Biobase)
library(BisqueRNA)
library(NetBID2)
library(scMINER)

pbmc.eset <- SeuratToExpressionSet(pbmc, delimiter='-', position=2, version="v3")

pbmc.eset.log2 <- generate.eset(exp_mat=exp.log2, 
                                phenotype_info=pData(pbmc.eset),
                                feature_info=fData(pbmc.eset),
                                annotation_info='exp.log2')

gene.symbol <- rownames(exprs(pbmc.eset.log2))
fData(pbmc.eset.log2) <- data.frame(gene.symbol)


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
  ranks<-names(sort(lab,decreasing = FALSE))
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



gn.sel<-c("CD3D","CD27","IL7R","SELL","CCR7","IL32","GZMA",
          "GZMK","DUSP2","CD8A","GZMH","GZMB","CD79A","CD79B","CD86","CD14")

# Heatmap
feature_heatmap(input_eset = pbmc.eset.log2, feature="gene.symbol", target = gn.sel, group_name = "cellType",
                save_plot = TRUE, width = 10, height = 5,
                name = "log2Exp", plot_name="/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_12k/Seurat/Heatmap_Tracy_markers_log2.pdf")



gn.sel <- c('CCR7', 'MAL', 'IL7R', 'AQP3', 'CD8B', 'NELL2', 'IL32', 'GZMK', 'KLRB1', 'CTSS', 'FCN1', 'SERPINA1', 
            'NKG7', 'GNLY', 'PRF1', 'GZMB', 'CD79A', 'CD79B', 'MS4A1', 'CST3', 'HLA-DRB1')

# Heatmap
feature_heatmap(input_eset = pbmc.eset.log2, feature="gene.symbol", target = gn.sel, group_name = "cellType",
                save_plot = TRUE, width = 10, height = 5,
                name = "log2Exp", plot_name="/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_12k/Seurat/Heatmap_21_markers_log2.pdf")




# Cell type heatmap
celltypes <- c("NaiveT", "CD4Mem", "CD8Mem", 
               "CD8EffMem", "Mono/DC", 
               "NK", "Bcell", "pDC", "Platelet")
ac <- matrix(NA, nrow=ncol(exp.log2), ncol=length(celltypes), dimnames = list(colnames(exp.log2), celltypes))


# Adam's markers
celltypes.dup <- c('NaiveT', 'NaiveT', 'CD4Mem', 'CD4Mem', 'CD4Mem', 'CD8Mem', 'CD8Mem', 'CD8EffMem', 'CD8EffMem', 'CD8EffMem',
                   'Mono/DC', 'Mono/DC', 'Mono/DC', 'NK', 'NK', 'Bcell', 'Bcell', 'pDC', 'pDC', 'Platelet', 'Platelet')
markers       <- c('CCR7', 'LEF1', 'AQP3', 'IL7R', 'IL32', 'GZMK', 'CD8B', 'CD8A', 'GZMH', 'CD8B', 
                   'CTSS', 'LYZ', 'CST3', 'KLRD1', 'GNLY', 'CD79A', 'CD79B', 'CCDC50', 'IRF7', 'GNG11', 'PPBP')
ref=data.frame(celltype=celltypes.dup, markers=markers, weight=rep(c(1),times=21))


# Tracy's markers
celltypes.dup <- c('NaiveT', 'NaiveT', 'NaiveT', 'NaiveT', 'CD4Mem', 'CD4Mem', 'CD8Mem', 'CD8Mem', 'CD8Mem', 'CD8Mem',
                          'CD8EffMem', 'CD8EffMem', 'CD8EffMem', 'CD8EffMem', 'CD8EffMem', 'NK', 'NK', 'Bell', 'Bell', 
                          'Mono/DC', 'Mono/DC', 'pDC', 'pDC', 'Platelet', 'Platelet')
markers <- c('CD3D', 'CD27', 'SELL', 'CCR7', 'CD3D', 'IL32', 'IL32', 'GZMA', 'GZMK', 'DUSP2',
                    'CD3D', 'IL32', 'GZMA', 'GZMH', 'GZMB', 'GZMA', 'GZMB', 'CD79A', 'CD79B',
                    'CD86', 'CD14', 'CCDC50', 'IRF7', 'GNG11', 'PPBP')
ref=data.frame(celltype=celltypes.dup, markers=markers, weight=rep(c(1),times=25))





library(openxlsx)
ref<-read.xlsx("/Users/lding/Git/scMINER/docs/tests/Ref/Azimuth_PBMC_markers.xlsx",sheet="Markers_8_cell_types")
head(ref)

celltypes <- c("CD4Naive", "CD8Naive", "CD4TCM", "CD8TCM", "CD4TEM", "CD8TEM", "Monocyte", "NK", 
               "Bcell", "DC")

celltypes <- c("NaiveT", "CD4TCM", "CD8TCM", "TEM", "Monocyte", "NK", "Bcell", "DC")

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



pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_12k/Seurat/Heatmap_marker_score_Azimuth_markers_8_cell_types.pdf", 
    height = 10, width=10)
Heatmap(as.matrix(t(input[,-1])), col = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)), name = 'MarkerScore',
               show_row_names = TRUE,
               show_column_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,)
dev.off()
