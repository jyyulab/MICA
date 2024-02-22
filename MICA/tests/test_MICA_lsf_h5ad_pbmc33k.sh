#!/usr/bin/env bash

mica lsf \
-i /research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/datasets/pbmc33k_preprocessed.h5ad \
-p "cwl_lsf" \
-k 3 4 \
-o /research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/outputs \
-j ./MICA/config/config_cwlexec.json
