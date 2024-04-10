library(scMINER)
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
pbmc <- FindClusters(pbmc, resolution = 0.16)

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



new.cluster.ids <- c("Naive T (CCR7+,LEF1+)", "CD4+ Memory T (AQP3+,IL7R+,IL32+)", "CD8+ Memory T (GZMK+,CD8B+)", "CD8+ Effector Memory T (CD8A+,GZMH+,CD8B+)", "Monocyte/DC (CTSS+,LYZ+,CST3+)", "NK (KLRD1+,GNLY+)", 
                     "B cell (CD79A+, CD79B+)", "pDC (CCDC50+, IRF7+)", "Platelet (GNG11+, PPBP+)")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()




# Violin plot, Tracy's marker genes
library(Seurat)
library(Biobase)
library(BisqueRNA)
pbmc.eset <- SeuratToExpressionSet(pbmc, delimiter='-', position=2, version="v3")

gn.sel<-c("CD3D","CD27","IL7R","SELL","CCR7","IL32","GZMA",
          "GZMK","DUSP2","CD8A","GZMH","GZMB","CD79A","CD79B","CD86","CD14")


#' @title feature_vlnplot
#' @description This plot will visualize feature info in violin plot by outputing a ggplot object
#' @param input_eset Input expression set
#' @param feature character, which feature to visualize
#' @param target a character vector, the list of feature to visualize
#' @param stat a character, whether to plot median or mean as a black dot on violinplot
#' @param group_by character, which group info to visualize as x axis
#' @param color_by character, which group info to define color, if NULL, then violin plots will be colored by 'group_by'
#' @param colors character vector, default as NULL, will use ggplot default color palette
#' @param ylabel a character, title of y axis
#' @param boxplot logical, whether to plot boxplot on violinplot
#' @param title.size numerical, default as 5
#' @param ncol cordinates for y axis
#'
#' @export
feature_vlnplot <- function(input_eset,
                            target=NULL,feature="geneSymbol",
                            group_by="celltype",ylabel="Expression",
                            color_by=NULL,colors=NULL,
                            ncol=3,stat="median",
                            boxplot=FALSE,title.size=5){
  
  if(!group_by%in% colnames(pData(input_eset))) stop('Please check your group_by information!','\n')
  if(!feature%in% colnames(fData(input_eset))) stop('Please check your feature information!','\n')
  
  # extract input information
  input <- exprs(input_eset)
  indx<-which(fData(input_eset)[,feature]%in%target)
  gn<-fData(input_eset)[,feature][indx]
  
  if(length(indx)==0) stop('No target feature found in data!','\n')
  
  label <- as.factor(pData(input_eset)[,group_by])
  if (is.null(color_by)) color_by=group_by
  condition<- as.factor(pData(input_eset)[,color_by])
  
  # Gene expression visualized as columns
  if (length(target)!=1) {
    target_values <- t(as.matrix(input[indx,]))
    colnames(target_values)<-gn
    df <- data.frame(target_values,cluster=label,condition=condition)
  }else {
    target_values<-input[indx,]
    df <- data.frame(target_values,cluster=label,condition=condition)
    colnames(df)[1]<-gn
  }
  
  df_melt <- reshape2::melt(df, id.vars=c("cluster","condition"))
  
  p <- ggplot(df_melt, aes(x=cluster, y=value, fill=condition))+
    theme_classic()+
    geom_violin(trim=TRUE,scale="width",na.rm = TRUE,size=0.4)
  
  if(!is.null(stat)){
    if (stat=="median") p <- p + stat_summary(fun=median, geom="point", size=1.2, color="black",position=position_dodge(width=1))
    else if (stat=="mean") p <- p + stat_summary(fun=mean, geom="point", size=1.2, color="black",position=position_dodge(width=1))
    else cat("Stat not supported, please check your spelling.","\n")}
  
  if(boxplot) p <- p + geom_boxplot(fill="white",width=0.1,size=0.1,outlier.size = 0.001,show.legend = FALSE,na.rm = TRUE)
  
  p <- p + facet_wrap(~variable,scales = "free",ncol = ncol) +
    labs(x=group_by,y=ylabel)+
    theme(axis.text.x = element_text(size=10),
          plot.title = element_text(size = title.size, face = "bold"),
          strip.background = element_rect(fill="#FFFFFF"))
  
  if (!is.null(colors)) p <- p+ scale_fill_manual(values=colors)
  
  if (ylabel=="Activity") { p <- p + geom_boxplot(width=0.2,size=0.1,outlier.size = 0.001,show.legend = FALSE,na.rm = TRUE)}
  
  return(p)
}



gene.symbol <- rownames(exprs(pbmc.eset))
fData(pbmc.eset) <- data.frame(gene.symbol)

# Violin
p <- feature_vlnplot(input_eset=pbmc.eset, target=gn.sel, feature = "gene.symbol",
                     group_by = "cellType", ylabel = "log2Exp", ncol = 4)
p


# Heatmap
feature_heatmap(input_eset = pbmc.eset, feature="gene.symbol", target = gn.sel, group_name = "cellType",
                save_plot = TRUE, width = 10, height = 5,
                name = "log2Exp", 
                plot_name="/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_12k/Seurat/Heatmap_Tracy_markers.pdf")








#' Visualize marker score of different cell types on bubbleplot
#'
#' @title Generate visualization for marker scores via bubble plot
#' @description  Marker visualizatoin from known markers/signatures, requires knowledge-based marker list as input
#' @param ref reference dataframe, includes positive or negative markers for different cell types;
#' Specify first column as different cell types, second columns as markers, third columns as weight (postive or negative marker)
#' @param input_eset expressionSet/SparseExpressionSet object with clustering membership stored in pData
#' @param group_name a character, the variable containing clustering label in pData(eset); or any other group information stored in pData(eset)
#' @param save_plot logical, whether or not save your plot; if TRUE, plot will be saved as plot_name
#' @param width default as 8, inch as unit
#' @param height default as 5, inch as unit
#' @param plot_name plot name, please include plot type
#' @param feature feature type from second column of your reference , should be in colnames(fData(eset))
#' @return A ggplot object
#' @examples
#' \dontrun{
#' df.ref=data.frame(celltype="Cd4 T",markers=c("Cd8a","Cd4","Cd3g"),weight=c(-1,1,1))
#' draw.marker.bbp<-function(ref = df.ref,input_eset, feature='geneSymbol',group_name="ClusterRes", save_plot = FALSE, width=8, height=5)
#'
#'
#' }
#'
#' @export
draw.marker.bbp<-function(ref = NULL,input_eset,
                          feature='geneSymbol',group_name="ClusterRes",
                          save_plot = FALSE,
                          width=8, height=5,
                          plot_name="AnnotationBubbleplot.png"){
  
  #exp<-apply(exprs(eset),2,std)
  #filter reference marker sets
  
  if (!feature%in%colnames(fData(input_eset))) stop('Please check your feature!')
  colnames(ref)<-c("celltype","markers","weight")
  ref<-dplyr::filter(ref,markers%in%fData(input_eset)[,feature])
  indx<-which(fData(input_eset)[,feature]%in%ref$markers)
  if(length(indx)==0) stop("No genes from the reference list could be found in data!","\n")
  
  exp<-as.matrix(exprs(input_eset))[indx,]
  rownames(exp)<-fData(input_eset)[,feature][indx]
  
  celltypes<-unique(ref$celltype)
  # print(celltypes)
  
  ac<-matrix(NA,nrow=ncol(exp),ncol=length(celltypes),dimnames = list(colnames(exp),celltypes))
  for(i in 1:length(celltypes)){
    cat(i,"\n")
    ref.sel<-dplyr::filter(ref,celltype==celltypes[i])
    n <- length(unique(ref.sel$markers))
    
    if(n>1){
      mat<-t(exp[ref.sel$markers,])%*%as.vector(ref.sel$weight)
      ac[,i]<-mat[,1]/n
    }else if (n==1){
      ac[,i]<-exp[ref.sel$markers,]
    }
  }
  
  print(str(ac))
  
  ac_norm<-apply(ac,2,scale) #column normalization
  
  n_mtx<-(ac>0.5)
  df_n<-data.frame(label=pData(input_eset)[,group_name],n_mtx)
  df_n<-aggregate(.~label,df_n,mean)
  print(df_n)
  library(reshape2)
  df_n_melt<-melt(df_n,id.vars = "label")
  print(df_n_melt)
  print('hello1')
  
  df<-data.frame(label=pData(input_eset)[,group_name],ac_norm);
  print(str(df))
  df<-df[,colSums(is.na(df))<nrow(df)];   # remove NA columns
  df<-aggregate(.~label,df,mean)
  input<-t(apply(df[,-1],1,scale))        # row normalization
  input<-as.data.frame(cbind(df[,1],input))
  rownames(input)<-rownames(df)
  colnames(input)<-colnames(df)
  print('hello2')
  # print(input)
  df_melt<-melt(input,id.vars = "label")
  print(str(df_melt))
  print('hello3')
  print(df_melt[,c(1,2)])
  print(df_n_melt[,c(1,2)])
  
  print('hello4')
  if(all(df_melt[,c(1,2)]==df_n_melt[,c(1,2)])){
    d<-cbind(df_melt,df_n_melt[,3])
    print('hello5')
    print(d)
    colnames(d)<-c("Cluster", "CellType", "MarkerScore","ExpressionPercentage")
    d$Cluster<-as.factor(d$Cluster)
    
    p<-draw.bubblePlot2(df=d, xlab="Cluster", ylab="CellType",
                        clab="MarkerScore", slab="ExpressionPercentage",
                        plot.title="Cell type annotation for each cluster")
  }
  
  
  if(save_plot){ggsave(plot = p, filename = plot_name , units="in",
                       width = width,height = height,dpi = 300)}
  return(p)
}






# Visualize marker score of different cell types on bubbleplot
df.ref=data.frame(celltype=c("Naive T (CCR7+,LEF1+)", "CD4+ Memory T (AQP3+,IL7R+,IL32+)"),
                  markers=c("CCR7","LEF1"),
                  weight=c(1,1))

draw.marker.bbp(ref = df.ref, input_eset = pbmc.eset,
                          feature='gene.symbol',group_name="cellType",
                          save_plot = TRUE,
                          width=8, height=5,
                          plot_name="/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_12k/Seurat/AnnotationBubbleplot_Tracy_markers.png")










exp <- exprs(pbmc.eset)
celltypes=c("Naive T (CCR7+,LEF1+)", "CD4+ Memory T (AQP3+,IL7R+,IL32+)")
markers=c("CCR7","LEF1")
ac <- matrix(NA,nrow=ncol(exp),ncol=length(celltypes),dimnames = list(colnames(exp),celltypes))

ref=data.frame(celltype=c("Naive T (CCR7+,LEF1+)", "CD4+ Memory T (AQP3+,IL7R+,IL32+)"),
                  markers=c("CCR7","LEF1"),
                  weight=c(1,1))

for(i in 1:length(celltypes)){
  cat(i,"\n")
  ref.sel<-dplyr::filter(ref,celltype==celltypes[i])
  n <- length(unique(ref.sel$markers))
  
  if(n>1){
    mat<-t(exp[ref.sel$markers,])%*%as.vector(ref.sel$weight)
    ac[,i]<-mat[,1]/n
  }else if (n==1){
    ac[,i]<-exp[ref.sel$markers,]
  }
}



# Heatmap
feature_heatmap(input_eset = pbmc.eset, feature="gene.symbol", target = gn.sel, group_name = "cellType",
                save_plot = TRUE, width = 10, height = 5,
                name = "log2Exp", plot_name="/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_12k/Seurat/Heatmap_Tracy_markers.pdf")

