# Seurat clustering and singleR annotation

library(Seurat)
library(dplyr)
library(patchwork)
library(SingleR)
library(celldex)




# Seurat
pbmc.12k <- Read10X(data.dir = "/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_12k/PBMC12k_input")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.12k, project = "pbmc12k", min.cells = 3, min.features = 200)


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)


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



pbmc.markers <- FindAllMarkers(pbmc, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
markers_plus <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 40, order_by = avg_log2FC)

markers_minus <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_min(n = 40, order_by = avg_log2FC)



# Make maker list heatmap plot
# Seurat clustering
exp.count <- as.matrix(pbmc@assays$RNA@counts)

# CPM normalization
norm = 1e5
exp.norm <- sweep(exp.count, 2, norm/unname(Matrix::colSums(exp.count)), '*')

# log transformation
exp.log2 <- log(exp.norm+1, base=2)



# SingleR reference
ref <- DatabaseImmuneCellExpressionData()
pred.labels <- SingleR(test = exp.log2, ref = ref, assay.type.test=1,
                      labels = ref$label.main)


pdf(file = "/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3/PBMC_12k/Seurat/Heatmap_singleR_score.pdf", 
    height = 10, width=10)
Heatmap(pred.labels$scores, col = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)), name = 'SingleR Score',
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,)
dev.off()

