#!/usr/bin/env Rscript


library(Seurat)

mc_human = readRDS(file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/Motor_Cortex_human/allen_m1c_2019_ssv4.rds')


pbmc <- RunUMAP(mc_human, dims = 1:10)
DimPlot(mc_human, reduction = "umap")
