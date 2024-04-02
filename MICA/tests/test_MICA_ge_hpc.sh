#!/usr/bin/env bash

#BSUB -P MICA_GE
#BSUB -q superdome
#BSUB -oo MICA_GE.out -eo MICA_GE.err
#BSUB -n 22
#BSUB -R "span[hosts=1]"
#BSUB -M 3000
mica ge -i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad -o ./test_data/outputs -nw 20
