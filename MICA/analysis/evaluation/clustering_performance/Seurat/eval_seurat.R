#/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)


raw_counts<-read.table(file=paste0("/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/Filtered_DownSampled_SortedPBMC_data.csv"), 
                       sep=",", header=TRUE, row.names=1)

raw_counts_t <- transpose(raw_counts)
rownames(raw_counts_t) <- colnames(raw_counts)
colnames(raw_counts_t) <- rownames(raw_counts)




pbmc <- CreateSeuratObject(counts = raw_counts_t, project = "pbmc20k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)


pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
write.csv(Idents(pbmc), file='/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/seurat_predicted_label.csv')

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
write.csv(pbmc@reductions$umap@cell.embeddings, file='/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/seurat_umap.csv')

pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")

saveRDS(pbmc, file = "/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/Filtered_DownSampled_SortedPBMC_data.rds")
