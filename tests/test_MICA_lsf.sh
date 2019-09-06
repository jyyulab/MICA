#!/usr/bin/env bash

mica lsf \
-i ./test_data/inputs/PBMC_Demo_MICA_input_mini.txt \
-p "cwl_lsf" \
-k 3 4 \
-o ./test_data/outputs/cwl_lsf/ \
-j ./MICA/config/config_cwlexec.json
