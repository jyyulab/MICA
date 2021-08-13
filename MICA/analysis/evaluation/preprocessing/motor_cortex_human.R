#!/usr/bin/env Rscript


library(Seurat)

pollen = readRDS(file = '/Users/lding/Documents/MICA/Datasets/HPC/GoldernStd/Pollen/pollen.rds')


pbmc <- RunUMAP(mc_human, dims = 1:10)
DimPlot(mc_human, reduction = "umap")
