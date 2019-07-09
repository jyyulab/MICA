#!/usr/bin/env bash

export PATH=`pwd`/scMINER/bin:$PATH
./scMINER/scMINER.py lsf \
-i /home/cqian/PBMC12K/PBMC12k_MICA_input.txt \
-p "pbmc_12k" \
-k 8 9 10 \
-o test_data/outputs/pbmc_12k/ \
-q standard \