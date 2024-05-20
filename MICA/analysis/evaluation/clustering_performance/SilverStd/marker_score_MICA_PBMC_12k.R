# MICA clustering and marker score
library(scMINER)



d.12k <- readscRNAseqData(file = "/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_12k/PBMC12k_input",
                          is.10x = T, CreateSparseEset = F, add.meta = F)

eset.12k<-CreateSparseEset(data=d.12k$raw.data, feature.data = d.12k$feature.data, add.meta = T)


cutoffs <- draw.scRNAseq.QC(SparseEset=eset.12k, 
                            project.name = "PBMC12k",
                            plot.dir = "/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_12k/PBMC12k_input",
                            group = "group", # this indicate which meta data information will be use in x axis to group violin plots
                            output.cutoff = TRUE) # whether or not to output suggested cutoffs




d.12k.v3 <- read.table(file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_12k/PBMC12k_MICA_input_v3.txt', 
            sep = '\t', header = TRUE, row.names = 1)

d.12k.v3 <- as.matrix(d.12k.v3)

eset.12k.v3 <- eset.12k[fData(eset.12k)$ensembl %in% colnames(d.12k.v3), pData(eset.12k)$cellName %in% rownames(d.12k.v3)]






# do filtering
# eset.sel <- preMICA.filtering(SparseEset = eset.12k, cutoffs = cutoffs)

norm = 1e6
exp.norm <- sweep(exprs(eset.12k.v3), 2, norm/unname(Matrix::colSums(exprs(eset.12k.v3))), '*')

# log transformation
# Required for MICA
exp.log2 <- log(exp.norm+1, base=2)

# save as SparseEset
eset.log2.v3 <- CreateSparseEset(data=exp.log2, meta.data = pData(eset.12k.v3), feature.data = fData(eset.12k.v3), add.meta = F)




cluster.res <- read.table(file = '/Users/lding/Documents/MICA/outputs/PBMC_12k/MICA_GE/clustering_UMAP_euclidean_20_2.01.txt', 
                          sep = '\t', header = TRUE)

pData(eset.log2.v3)$cellType <- cluster.res$label






library(openxlsx)
ref<-read.xlsx("/Users/lding/Git/scMINER/docs/tests/Ref/Azimuth_PBMC_markers.xlsx",sheet="Markers_8_cell_types")
head(ref)

celltypes <- c("CD4Naive", "CD8Naive", "CD4TCM", "CD8TCM", "CD4TEM", "CD8TEM", "Monocyte", "NK", 
               "Bcell", "DC")

celltypes <- c("NaiveT", "CD4TCM", "CD8TCM", "TEM", "Monocyte", "NK", "Bcell", "DC")

rownames(exp.log2) <- fData(eset.12k.v3)$geneSymbol
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
df_n<-data.frame(label=pData(eset.log2.v3)[,group_name],n_mtx)
df_n<-aggregate(.~label,df_n,mean)
library(reshape2)
df_n_melt<-melt(df_n,id.vars = "label")


df <- data.frame(label=pData(eset.log2.v3)[,group_name],ac_norm);
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



pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_12k/MICA/Heatmap_marker_score_Azimuth_markers_8_cell_types.pdf", 
    height = 10, width=10)

h.mat <- t(input[,-1])
h.mat <- h.mat[,c(5,2,3,8,4,7,1,6)]

Heatmap(as.matrix(h.mat), col = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)), name = 'MarkerScore',
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,)
dev.off()





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
feature_vlnplot(input_eset=eset.log2.v3,target=gn.sel,feature = "geneSymbol",
                group_by = "cellType",ylabel = "log2Exp",ncol = 4)
dev.off()




# Heatmap
feature_heatmap(input_eset = eset.log2.v3, feature="geneSymbol", target = gn.sel, group_name = "cellType",
                save_plot = TRUE, width = 10, height = 5,
                name = "log2Exp", plot_name="/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_12k/MICA/Heatmap_21_markers_log2.pdf")
