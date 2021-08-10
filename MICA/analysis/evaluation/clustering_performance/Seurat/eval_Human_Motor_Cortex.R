#/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(pdfCluster)


data_path <- "/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/datasets/SilverStd"
raw_counts<-read.table(file=paste0(data_path, "/Human_Motor_Cortex/matrix.csv"), sep=",", header=TRUE, row.names=1)

raw_counts_t <- transpose(raw_counts)
rownames(raw_counts_t) <- colnames(raw_counts)
colnames(raw_counts_t) <- rownames(raw_counts)



cortex <- CreateSeuratObject(counts = raw_counts_t, project = "cortex", min.cells = 3, min.features = 200)
cortex[["percent.mt"]] <- PercentageFeatureSet(cortex, pattern = "^MT-")
# VlnPlot(cortex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cortex <- subset(cortex, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 5)
cortex <- NormalizeData(cortex, normalization.method = "LogNormalize", scale.factor = 10000)
cortex <- FindVariableFeatures(cortex, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cortex)
cortex <- ScaleData(cortex, features = all.genes)
cortex <- RunPCA(cortex, features = VariableFeatures(object = cortex))
cortex <- JackStraw(cortex, num.replicate = 100)
cortex <- ScoreJackStraw(cortex, dims = 1:20)
JackStrawPlot(cortex, dims = 1:15)
ElbowPlot(cortex)


cortex <- FindNeighbors(cortex, dims = 1:10)
cortex <- FindClusters(cortex, resolution = 1.0)
head(Idents(cortex), 5)

output_path <- '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/outputs/SilverStd/Human_Motor_Cortex/seurat'
write.csv(Idents(cortex), file=paste0(output_path, '/seurat_predicted_label.csv'))

cortex <- RunUMAP(cortex, dims = 1:10)

pdf(file=paste0(output_path, "/seurat_umap.pdf"))
DimPlot(cortex, reduction = "umap")
dev.off()

write.csv(cortex@reductions$umap@cell.embeddings, file=paste0(output_path, '/seurat_umap.csv'))

cortex <- RunTSNE(cortex, dims = 1:10)

pdf(file=paste0(output_path, "/seurat_tsne.pdf"))
DimPlot(cortex, reduction = "tsne")
dev.off()

saveRDS(cortex, file = paste0(output_path, "/Human_Motor_Cortex.rds"))



# Calculate ARI
true_label_file <- paste0(data_path, '/Human_Motor_Cortex/Human_Motor_Cortex_true_label.txt')
true_labels <- read.table(file=true_label_file, sep="\t", header=TRUE, row.names=1)
merged <- merge(true_labels, data.frame(cortex$seurat_clusters), by=0)
adj.rand.index(merged$class_label, merged$cortex.seurat_clusters)
