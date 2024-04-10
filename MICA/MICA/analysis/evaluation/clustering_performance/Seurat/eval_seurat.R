#/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)


raw_counts <- read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_MICA_input_filter_14k.txt"), 
                       sep="\t", header=TRUE, row.names=1)

raw_counts_t <- transpose(raw_counts)
rownames(raw_counts_t) <- colnames(raw_counts)
colnames(raw_counts_t) <- rownames(raw_counts)




pbmc <- CreateSeuratObject(counts = raw_counts_t, project = "pbmc14k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 15280)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)


pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.27)
head(Idents(pbmc), 5)
write.csv(Idents(pbmc), file='/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/seurat_predicted_label_14k.csv')

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
write.csv(pbmc@reductions$umap@cell.embeddings, file='/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/seurat_umap.csv')

pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")

FeaturePlot(pbmc, features = c("FOXP3"))

saveRDS(pbmc, file = "/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/Filtered_DownSampled_SortedPBMC_data.rds")






meta.data<-pbmc@meta.data
feature.data<-data.frame(rownames(pbmc@assays$RNA@data))
colnames(feature.data) <- "geneSymbol"
rownames(feature.data) <- feature.data$geneSymbol
eset<-CreateSparseEset(data=pbmc@assays$RNA@data, meta.data = meta.data, feature.data = feature.data, add.meta = F)


true.label.20k <- read.table(file='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt', 
                             sep = '\t', header = T)
true.label.14k <- true.label.20k[true.label.20k$cell %in% rownames(pData(eset)),]
pData(eset)$trueLabel <- as.factor(true.label.14k$label)
pData(eset)$X <- pbmc@reductions$umap@cell.embeddings[,'UMAP_1']
pData(eset)$Y <- pbmc@reductions$umap@cell.embeddings[,'UMAP_2']



colors <- c('#C3C3C3', '#C3C3C3', '#F8766D', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3')
MICAplot(input_eset = eset,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "trueLabel", colors = colors, pct = 0.5)



colors <- c('#C3C3C3', '#C3C3C3', '#CD9600', '#C3C3C3', '#C3C3C3', '#C3C3C3', '#C3C3C3')
MICAplot(input_eset = eset,
         X="X", Y="Y", # which meta variable was treated as x or y coordinates
         color_by = "trueLabel", colors = colors, pct = 0.5)
