#!/usr/bin/env bash

mica -pl lsf \
-i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad \
-o ./test_data/outputs/cwl_lsf/ \
-pn "cwl_lsf" \
-nc 3 4 \
-cj ./MICA/config/config_cwlexec.json
