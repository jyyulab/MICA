#!/bin/env bash
#BSUB -P DN
#BSUB -J eval_Human_Motor_Cortex
#BSUB -oo eval_Human_Motor_Cortex.out 
#BSUB -eo eval_Human_Motor_Cortex.err 
#BSUB -q superdome   
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -M 150000
/research/rgs01/home/clusterHome/lding/Git/MICA/MICA/analysis/evaluation/clustering_performance/Seurat/eval_Human_Motor_Cortex.R

