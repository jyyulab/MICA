#!/usr/bin/env python3
from MICA.analysis.evaluation.clustering_performance import utils


root_dir = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'
level = 'SilverStd'
author = 'Human_Motor_Cortex'
input_file = '{}_MICA_input.h5ad'.format(author)
num_clusters = 20

#%%
max_resolution = 10.0
# utils.create_cmds_ge(root_dir, level, author, input_file, max_resolution)
utils.calc_ARIs_ge(root_dir, level, author, num_clusters)
