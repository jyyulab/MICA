#!/usr/bin/env bash

# ./MICA/mica_ge.py -i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad -o ./test_data/outputs -d 12 -e 1.0

./MICA/mica_u.py -i /Users/lding/Documents/MICA/Datasets/filtered_gene_bc_matrices/hg19/pbmc33k_preprocessed.h5ad -o /Users/lding/Documents/MICA/Datasets/filtered_gene_bc_matrices/hg19 -m euclidean -d 10

