#!/usr/bin/env python3
from MICA.analysis.evaluation.clustering_performance import utils


root_dir = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'
level = 'SilverStd'
author = 'Tasic'
input_file = '{}_MICA_input.txt'.format(author)
num_clusters = 8
# utils.create_cmds_mds(root_dir, level, author, input_file, num_clusters)
utils.calc_ARIs_mds(root_dir, level, author, num_clusters)

# max_resolution = 10.0
# utils.create_cmds_ge(root_dir, level, author, input_file, max_resolution)
utils.calc_ARIr_ge(root_dir, level, author, num_clusters)


#%%
# from MICA.lib import preprocessing as pp
# root_dir = '/Users/lding/Documents/MICA/Datasets/HPC/'
# adata = pp.read_preprocessed_mat('{}/{}/{}/{}'.format(root_dir, level, author, input_file))
# frame = adata.to_df()
# frame.to_csv('{}/{}/{}/{}_MICA_input.txt'.format(root_dir, level, author, author), sep='\t')
