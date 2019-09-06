#!/usr/bin/env bash

mica local \
-i ./test_data/inputs/PBMC_Demo_MICA_input_mini.txt \
-p "cwl_local" \
-k 3 4 \
-o ./test_data/outputs/cwl_local/ \
--dist "spearman" \
