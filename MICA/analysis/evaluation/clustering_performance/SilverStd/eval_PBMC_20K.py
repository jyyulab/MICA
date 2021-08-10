#!/usr/bin/env python3
from MICA.analysis.evaluation.clustering_performance import utils


root_dir = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'
level = 'SilverStd'
author = 'PBMC_20k'
input_file = '{}_MICA_input.h5ad'.format(author)
num_clusters = 10
# utils.create_cmds_mds(root_dir, level, author, input_file, num_clusters)
# utils.calc_ARIs_mds(root_dir, level, author, num_clusters)

# max_resolution = 10.0
# utils.create_cmds_ge(root_dir, level, author, input_file, max_resolution)
# utils.calc_ARIs_ge(root_dir, level, author, num_clusters)


# Check an old commit
# import pandas as pd
# from sklearn.metrics.cluster import adjusted_rand_score
# output_dir = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/outputs/SilverStd/PBMC_20k/ge_debug'
#
# true_label_file = '/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/datasets/' \
#                   'SilverStd/PBMC_20k/PBMC_20k_true_label.txt'
# true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
# predict_label_file = '{}/clustering_UMAP_euclidean_{}_{}.txt'.format(output_dir, 12, 2.2)
# predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
# predict_num_clusters = len(set(predict_label['label']))
# print(predict_num_clusters)
# merged = true_label.merge(predict_label, left_on='cell', right_index=True)
# ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
# print(ari)


#%%
from MICA.lib import preprocessing as pp
root_dir = '/Users/lding/Documents/MICA/Datasets/HPC/'
adata = pp.read_preprocessed_mat('{}/{}/{}/{}'.format(root_dir, level, author, input_file))
frame = adata.to_df()
frame.to_csv('{}/{}/{}/{}_MICA_input.txt'.format(root_dir, level, author, author), sep='\t')
