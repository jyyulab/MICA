#!/usr/bin/env bash

mica local \
-i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad \
-p "cwl_local" \
-k 3 4 \
-o ./test_data/outputs/cwl_local/ \
--dist "mi" \
