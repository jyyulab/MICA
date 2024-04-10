#!/usr/bin/env bash

#BSUB -P MICA_GE
#BSUB -q superdome
#BSUB -oo MICA_GE.out -eo MICA_GE.err
#BSUB -n 26
#BSUB -R "span[hosts=1]"
#BSUB -M 3000
mica ge -i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad -o ./test_data/outputs/dr_cell -da cell -dd 100 -nw 25
