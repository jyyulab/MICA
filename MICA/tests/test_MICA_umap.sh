#!/usr/bin/env bash

./MICA/mica_umap.py -i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad -o ./test_data/outputs -d 10

# ./MICA/mica_f.py -i /Users/lding/Documents/MICA/Datasets/filtered_gene_bc_matrices/hg19/pbmc33k_preprocessed.h5ad -o /Users/lding/Documents/MICA/Datasets/filtered_gene_bc_matrices/hg19 -j 4 -m mi -d 40
