#!/usr/bin/env python3

import pandas as pd
from sklearn.metrics.cluster import adjusted_mutual_info_score

#%%
predict_label_file = '/Users/lding/Documents/MICA/Manuscript/Figures/Figure_S3_1/PBMC_14k/MICA/filter_CTL/after_parameter_tuning/clustering_UMAP_euclidean_24_1.82212.txt'

true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
predict_num_clusters = len(set(predict_label['label']))
print(predict_num_clusters)
merged = true_label.merge(predict_label, left_on='cell', right_index=True)
ami = adjusted_mutual_info_score(merged['label_x'], merged['label_y'])
print(ami)
