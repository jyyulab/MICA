#!/usr/bin/env bash

./MICA/mica_f.py -i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad -o ./test_data/outputs -j 4 -m mi -d 500

#./MICA/mica_f.py -i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad -o ./test_data/outputs -j 4 -m euclidean -d 40
