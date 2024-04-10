#!/usr/bin/env Rscript

pollen = readRDS(file = '/Users/lding/Documents/MICA/Datasets/HPC/GoldernStd/Pollen/pollen.rds')


pollen <- readscRNAseqData(file="/Users/lding/Documents/MICA/Datasets/HPC/GoldenStd/Pollen/Pollen_MICA_input.txt",
                           is.10x = F, CreateSparseEset = F, add.meta=F)
pollen_t <- t(as.matrix(pollen))
eset.12k<-CreateSparseEset(data=as.matrix(pollen_t))


cutoffs <- draw.scRNAseq.QC(SparseEset=eset.12k, 
                            project.name = "Pollen",
                            plot.dir = "/Users/lding/Desktop",
                            # group = "group", # this indicate which meta data information will be use in x axis to group violin plots
                            output.cutoff = TRUE) #whether or not to output suggested cutoffs

cutoffs$nCell_cutoff = 3

eset.sel <- preMICA.filtering(SparseEset = eset.12k, cutoffs = cutoffs) 
