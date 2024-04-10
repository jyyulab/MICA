#!/usr/bin/env bash

mica -pl local \
-i ./test_data/inputs/PBMC_Demo_MICA_input_mini.txt \
-o ./test_data/outputs/cwl_local/ \
-pn "cwl_local" \
-nc 3 4
