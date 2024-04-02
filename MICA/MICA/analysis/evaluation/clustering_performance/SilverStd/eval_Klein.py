#!/usr/bin/env python3
from MICA.analysis.evaluation.clustering_performance import utils


root_dir = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'
level = 'SilverStd'
author = 'Klein'
input_file = '{}_MICA_input.txt'.format(author)
num_clusters = 4
# utils.create_cmds_mds(root_dir, level, author, input_file, num_clusters)
# utils.calc_ARIs_mds(root_dir, level, author, num_clusters)

# max_resolution = 10.0
# utils.create_cmds_ge(root_dir, level, author, input_file, max_resolution)
# utils.calc_ARIs_ge(root_dir, level, author, num_clusters)

# utils.calc_AMIs_mds(root_dir, level, author, num_clusters)
utils.calc_AMIs_ge(root_dir, level, author, num_clusters)



#%%
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score

#%%
predict_label_file = '/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/after_tuning/Klein/clustering_UMAP_euclidean_20_4.95303.txt'

true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Klein/Klein_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
predict_num_clusters = len(set(predict_label['label']))
print(predict_num_clusters)
merged = true_label.merge(predict_label, left_on='cell', right_index=True)
ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
print(ari)
