#!/usr/bin/env bash

mica lsf \
-i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad \
-p "cwl_lsf" \
-k 3 4 \
-o ./test_data/outputs/cwl_lsf/ \
-j ./MICA/config/config_cwlexec.json
