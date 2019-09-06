#!/usr/bin/env bash

mica lsf \
-i /research/projects/yu3grp/scRNASeq/yu3grp/TracyQian/scMINER/dev/scMINER_dev/test_data/inputs/PBMC_Demo_MICA_input_mini.txt \
-p "cwl_lsf" \
-k 3 4 \
-o test_data/outputs/cwl_lsf/ \
-q standard
