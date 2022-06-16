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


# utils.calc_AMIs_mds(root_dir, level, author, num_clusters)
# utils.calc_AMIs_ge(root_dir, level, author, num_clusters)


#%%
# Check an old commit
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score

#%%
predict_label_file = '/Users/lding/Documents/MICA/outputs/PBMC_20k/MICA/clustering_UMAP_euclidean_24_1.8.txt'

true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
predict_num_clusters = len(set(predict_label['label']))
print(predict_num_clusters)
merged = true_label.merge(predict_label, left_on='cell', right_index=True)
ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
print(ari)


#%%
predict_label_file = '/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3_1/PBMC_14k/MICA/filter_CTL/after_parameter_tuning/clustering_UMAP_euclidean_24_1.82212.txt'

true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
predict_num_clusters = len(set(predict_label['label']))
print(predict_num_clusters)
merged = true_label.merge(predict_label, left_on='cell', right_index=True)
ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
print(ari)


#%%
predict_label_file = '/Users/lding/Documents/scMINER/PBMC14k_input/ac_mat_metacell/clustering_UMAP_euclidean_24_2.01375.txt'

true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
predict_num_clusters = len(set(predict_label['label']))
print(predict_num_clusters)
merged = true_label.merge(predict_label, left_on='cell', right_index=True)
ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
print(ari)


#%%
predict_label_file = '/Users/lding/Documents/scMINER/PBMC14k_input/ac_mat_metacell/clustering_UMAP_euclidean_24_3.00417.txt'

true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
predict_num_clusters = len(set(predict_label['label']))
print(predict_num_clusters)
merged = true_label.merge(predict_label, left_on='cell', right_index=True)
ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
print(ari)



#%%
# predict_label_file = '/Users/lding/Documents/scMINER/PBMC14k_input/activity_6863_drivers/clustering_UMAP_euclidean_24_4.0552.txt'
# predict_label_file = '/Users/lding/Documents/scMINER/PBMC14k_input/activity_6863_drivers/clustering_UMAP_euclidean_24_4.95303.txt'
predict_label_file = '/Users/lding/Documents/scMINER/PBMC14k_input/activity_6863_drivers/clustering_UMAP_euclidean_24_2.01375.txt'

# predict_label_file = '/Users/lding/Documents/scMINER/PBMC14k_input/activity_6863_drivers/clustering_UMAP_euclidean_24_2.01375.txt'
predict_label_file = '/Users/lding/Documents/scMINER/PBMC14k_input/activity_6863_drivers/MDS/scatter2/clustering_umap_euclidean_19.txt'

true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
predict_num_clusters = len(set(predict_label['label']))
print(predict_num_clusters)
merged = true_label.merge(predict_label, left_on='cell', right_index=True)
ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
print(ari)


#%%
# from MICA.lib import preprocessing as pp
# # root_dir = '/Users/lding/Documents/MICA/Datasets/HPC/'
# adata = pp.read_preprocessed_mat('{}/{}/{}/{}'.format(root_dir, level, author, input_file))
#
# #%%
# frame = adata.to_df()
#
# #%%
# frame_T = frame.T
# # frame.to_csv('{}/{}/{}/{}_MICA_input.txt'.format(root_dir, level, author, author), sep='\t')
#
# #%%
# frame_T.insert(loc=0, column='geneSymbol', value=list(frame_T.index))
#
# #%%
# frame_T.to_csv('{}/{}/{}/{}_SJARACNe_input.txt'.format(root_dir, level, author, author), sep='\t')
