#!/usr/bin/env bash

mica lsf \
<<<<<<< HEAD
-i /research/projects/yu3grp/scRNASeq/yu3grp/TracyQian/scMINER/dev/scMINER_dev/test_data/inputs/PBMC_Demo_MICA_input_mini.txt \
=======
-i ./test_data/inputs/PBMC_Demo_MICA_input_mini.txt \
>>>>>>> 9b797370bb6cd79b4448be4440731461ef0e043a
-p "cwl_lsf" \
-k 3 4 \
-o ./test_data/outputs/cwl_lsf/ \
-q standard
